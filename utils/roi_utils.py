#!/usr/bin/env python3
"""
ROI utilities for mbrain mouse brain analysis project.

This module provides reusable functions for ROI manipulation, including:
- Voxel counting for NIfTI and MRtrix formats
- ROI size normalization utilities
- Batch processing for multiple subjects

Author: MBRAIN Project
Created: 2026-02-02
"""

import subprocess
from pathlib import Path
from typing import Union, Dict, List, Optional
import numpy as np

# Optional imports - will be loaded on demand
_nibabel_loaded = False
_nib = None


def _ensure_nibabel():
    """Lazy load nibabel to avoid import overhead when not needed."""
    global _nibabel_loaded, _nib
    if not _nibabel_loaded:
        import nibabel
        _nib = nibabel
        _nibabel_loaded = True
    return _nib


def count_voxels(roi_path: Union[str, Path], threshold: float = 0) -> int:
    """
    Count non-zero (or above threshold) voxels in an ROI file.
    
    Supports both NIfTI (.nii, .nii.gz) and MRtrix (.mif) formats.
    This is essential for ROI size normalization - larger ROIs generate more
    streamlines simply due to having more seed voxels.
    
    Args:
        roi_path: Path to ROI file (NIfTI or MIF format)
        threshold: Voxels with values > threshold are counted (default: 0)
    
    Returns:
        Number of non-zero voxels in the ROI
        
    Examples:
        >>> count_voxels('roi_thalamus_left.nii.gz')
        584
        >>> count_voxels('seg6_graft_left_anterior_dorsal_fod.mif')
        312
    """
    roi_path = Path(roi_path)
    
    if not roi_path.exists():
        return 0
    
    try:
        suffix = ''.join(roi_path.suffixes).lower()
        
        if '.mif' in suffix:
            # Use mrstats for MRtrix format (counts non-zero voxels)
            result = subprocess.run(
                ['mrstats', '-mask', str(roi_path), '-output', 'count', str(roi_path)],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                return int(float(result.stdout.strip()))
            return 0
        else:
            # NIfTI format - use nibabel
            nib = _ensure_nibabel()
            img = nib.load(roi_path)
            data = img.get_fdata()
            return int(np.sum(data > threshold))
            
    except Exception as e:
        print(f"  Warning: Could not count voxels in {roi_path}: {e}")
        return 0


def count_voxels_batch(roi_paths: List[Union[str, Path]], threshold: float = 0) -> Dict[str, int]:
    """
    Count voxels in multiple ROI files.
    
    Args:
        roi_paths: List of paths to ROI files
        threshold: Voxels with values > threshold are counted
    
    Returns:
        Dictionary mapping filename (stem) to voxel count
        
    Example:
        >>> roi_files = list(Path('segments/').glob('*.nii.gz'))
        >>> counts = count_voxels_batch(roi_files)
        >>> counts
        {'isocortex_left_anterior': 1234, 'isocortex_left_central': 987, ...}
    """
    results = {}
    for roi_path in roi_paths:
        path = Path(roi_path)
        # Use stem without all suffixes (e.g., 'roi_thalamus' from 'roi_thalamus.nii.gz')
        key = path.name.replace('.nii.gz', '').replace('.nii', '').replace('.mif', '')
        results[key] = count_voxels(path, threshold)
    return results


def get_roi_info(roi_path: Union[str, Path]) -> Dict:
    """
    Get comprehensive information about an ROI file.
    
    Args:
        roi_path: Path to ROI file
    
    Returns:
        Dictionary with voxel_count, volume_mm3, dimensions, voxel_size
    """
    roi_path = Path(roi_path)
    
    if not roi_path.exists():
        return {'error': f'File not found: {roi_path}'}
    
    try:
        suffix = ''.join(roi_path.suffixes).lower()
        
        if '.mif' in suffix:
            # Use mrinfo for MRtrix format
            voxels = count_voxels(roi_path)
            
            # Get voxel size
            result = subprocess.run(
                ['mrinfo', '-spacing', str(roi_path)],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                spacing = [float(x) for x in result.stdout.strip().split()]
                voxel_volume = np.prod(spacing)
            else:
                spacing = [0, 0, 0]
                voxel_volume = 0
            
            return {
                'voxel_count': voxels,
                'volume_mm3': voxels * voxel_volume,
                'voxel_size': spacing
            }
        else:
            # NIfTI format
            nib = _ensure_nibabel()
            img = nib.load(roi_path)
            data = img.get_fdata()
            voxels = int(np.sum(data > 0))
            
            # Get voxel dimensions from header
            voxel_size = img.header.get_zooms()[:3]
            voxel_volume = np.prod(voxel_size)
            
            return {
                'voxel_count': voxels,
                'volume_mm3': voxels * voxel_volume,
                'voxel_size': list(voxel_size),
                'dimensions': data.shape
            }
            
    except Exception as e:
        return {'error': str(e)}


def calculate_connectivity_density(streamline_count: int, roi_voxels: int, 
                                   scale_factor: float = 1000) -> float:
    """
    Calculate connectivity density normalized by ROI size.
    
    Formula: Connectivity Density = streamlines / roi_voxels × scale_factor
    
    This controls for ROI size bias - larger ROIs naturally generate more
    streamlines simply due to having more seed voxels.
    
    Args:
        streamline_count: Number of streamlines (or SIFT2 sum)
        roi_voxels: Number of voxels in the seed ROI
        scale_factor: Multiplier for readability (default 1000)
    
    Returns:
        Connectivity density value
        
    Example:
        >>> calculate_connectivity_density(5000, 500)
        10000.0  # 5000 streamlines from 500 voxels × 1000
    """
    if roi_voxels <= 0:
        return 0.0
    return (streamline_count / roi_voxels) * scale_factor


def count_streamlines(tck_path: Union[str, Path]) -> int:
    """
    Count streamlines in a .tck tractography file using tckinfo.
    
    Args:
        tck_path: Path to .tck file
    
    Returns:
        Number of streamlines
    """
    tck_path = Path(tck_path)
    
    if not tck_path.exists():
        return 0
    
    try:
        result = subprocess.run(
            ['tckinfo', '-count', str(tck_path)],
            capture_output=True, text=True
        )
        # Parse output - format is typically "count: N" or just "N"
        for line in result.stdout.strip().split('\n'):
            if 'count' in line.lower():
                return int(line.split(':')[1].strip())
            elif line.strip().isdigit():
                return int(line.strip())
    except Exception as e:
        print(f"  Warning: Could not count streamlines in {tck_path}: {e}")
    
    return 0


def get_subject_seed_voxels(subject_path: Path, group: str, 
                            segment_style: str = 'seg6') -> Dict[str, int]:
    """
    Get voxel counts for all seed segments of a subject.
    
    This is the main function for collecting seed ROI sizes for normalization.
    
    Args:
        subject_path: Path to subject's session directory
        group: 'CTRL', 'AP', or 'APTR'
        segment_style: 'seg3' (3 dorsal segments), 'seg6' (6 segments), or 'seg9' (9 segments)
    
    Returns:
        Dictionary mapping segment name to voxel count
        
    Example:
        >>> path = Path('data/derivatives/mrtrix_preproc_flip1/sub-129-921-R1F/ses-260115121231')
        >>> voxels = get_subject_seed_voxels(path, 'CTRL')
        >>> voxels
        {'isocortex_left_anterior': 1234, 'isocortex_left_central': 987, ...}
    """
    voxels = {}
    
    if group == 'CTRL':
        # CTRL uses isocortex segments (5 segments)
        seg_dir = subject_path / 'dwi' / 'isocortex_segments_v2'
        segments = ['isocortex_left_anterior', 'isocortex_left_anterior_central',
                    'isocortex_left_central', 'isocortex_left_posterior_central',
                    'isocortex_left_posterior']
        
        for seg in segments:
            seg_file = seg_dir / f'{seg}.nii.gz'
            if seg_file.exists():
                voxels[seg] = count_voxels(seg_file)
                
    elif group in ['AP', 'APTR']:
        # AP/APTR use graft segments
        # Try multiple possible locations
        possible_dirs = [
            subject_path / 'dwi' / 'graft_segments_fod',  # Current location
            subject_path / 'dwi' / 'graft_segments_v2' / 'fod_space',  # Alternative
        ]
        
        seg_dir = None
        for d in possible_dirs:
            if d.exists():
                seg_dir = d
                break
        
        if seg_dir is None:
            return voxels
        
        # Define segments based on style
        if segment_style == 'seg3':
            # Only dorsal segments (for cortico-thalamic connectivity)
            segments = ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']
        elif segment_style == 'seg6':
            segments = ['anterior_dorsal', 'anterior_ventral',
                        'central_dorsal', 'central_ventral',
                        'posterior_dorsal', 'posterior_ventral']
        else:  # seg9
            segments = ['anterior_dorsal', 'anterior_medial', 'anterior_ventral',
                        'central_dorsal', 'central_medial', 'central_ventral',
                        'posterior_dorsal', 'posterior_medial', 'posterior_ventral']
        
        for seg in segments:
            # Try different naming patterns used in the project
            for pattern in [
                f'seg9_graft_left_{seg}_fod.mif',  # Most common format
                f'seg6_graft_left_{seg}_fod.mif',
                f'graft_left_{seg}.mif',
                f'graft_left_{seg}_fod.mif',
            ]:
                seg_file = seg_dir / pattern
                if seg_file.exists():
                    voxels[f'graft_left_{seg}'] = count_voxels(seg_file)
                    break
    
    return voxels


def get_all_subjects_seed_voxels(config: dict, groups: Optional[List[str]] = None,
                                  segment_style: str = 'seg3') -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Get seed voxel counts for all subjects across specified groups.
    
    Args:
        config: Configuration dictionary from config_loader
        groups: List of groups to process (default: ['CTRL', 'AP', 'APTR'])
        segment_style: 'seg3', 'seg6', or 'seg9'
    
    Returns:
        Nested dictionary: {group: {subject_id: {segment: voxels}}}
        
    Example:
        >>> from functions.config_loader import load_config
        >>> config = load_config(pipeline='flip1')
        >>> voxels = get_all_subjects_seed_voxels(config)
        >>> voxels['CTRL']['sub-129-921-R1F']['isocortex_left_anterior']
        1234
    """
    if groups is None:
        groups = ['CTRL', 'AP', 'APTR']
    
    base_path = Path(config['paths']['derivatives_flip1'])
    results = {}
    
    for group in groups:
        if group not in config['groups']:
            continue
            
        results[group] = {}
        
        for subj_info in config['groups'][group]:
            subject_id = subj_info['id']
            session_id = subj_info['session']
            
            subject_path = base_path / subject_id / session_id
            
            if subject_path.exists():
                results[group][subject_id] = get_subject_seed_voxels(
                    subject_path, group, segment_style
                )
    
    return results


if __name__ == '__main__':
    # Demo/test code
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from functions.config_loader import load_config
    
    config = load_config(pipeline='flip1')
    
    print("=" * 80)
    print("ROI VOXEL COUNTS - ALL SUBJECTS")
    print("=" * 80)
    
    all_voxels = get_all_subjects_seed_voxels(config, segment_style='seg3')
    
    for group in ['CTRL', 'AP', 'APTR']:
        if group not in all_voxels:
            continue
        print(f"\n{group}:")
        print("-" * 60)
        
        for subject_id, segments in all_voxels[group].items():
            total = sum(segments.values())
            print(f"\n  {subject_id}: (total: {total})")
            for seg, count in sorted(segments.items()):
                print(f"    {seg}: {count}")
