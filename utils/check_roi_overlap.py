#!/usr/bin/env python
"""
Check ROI overlap between segmented regions.

This utility verifies that ROI segments (e.g., isocortex segments for CTRL,
graft segments for AP/APTR) do not overlap. Overlapping ROIs would cause
double-counting in tractography and connectivity analysis.

Usage:
    python utils/check_roi_overlap.py --subject sub-129-921-R1F --session ses-260115121231 --group CTRL
    python utils/check_roi_overlap.py --subject sub-AP72-862-901-NmF --session ses-251210093342 --group APTR
    python utils/check_roi_overlap.py --all  # Check all subjects

Author: MBRAIN Project
Date: 2026-01-30
"""

import argparse
import sys
import os
from pathlib import Path

import nibabel as nib
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config

# Use flip1 pipeline (config1.yaml)
os.environ['MBRAIN_PIPELINE'] = 'flip1'


def check_roi_overlap(roi_files: list, roi_names: list = None) -> dict:
    """
    Check for overlap between multiple ROI files.
    
    Parameters
    ----------
    roi_files : list
        List of paths to ROI NIfTI files
    roi_names : list, optional
        List of names for each ROI (for reporting). If None, uses filenames.
    
    Returns
    -------
    dict
        Dictionary with:
        - 'voxel_counts': dict of ROI name -> voxel count
        - 'overlaps': list of (roi1, roi2, overlap_count) tuples
        - 'total_unique': total unique voxels across all ROIs
        - 'total_sum': sum of all ROI voxel counts
        - 'has_overlap': bool, True if any overlap exists
    """
    if roi_names is None:
        roi_names = [Path(f).stem for f in roi_files]
    
    # Load all ROIs
    rois = {}
    voxel_counts = {}
    for name, path in zip(roi_names, roi_files):
        if not os.path.exists(path):
            print(f"  WARNING: {path} not found, skipping")
            continue
        img = nib.load(path)
        data = img.get_fdata() > 0
        rois[name] = data
        voxel_counts[name] = int(np.sum(data))
    
    # Check pairwise overlaps
    overlaps = []
    names = list(rois.keys())
    for i, name1 in enumerate(names):
        for name2 in names[i+1:]:
            overlap = int(np.sum(rois[name1] & rois[name2]))
            if overlap > 0:
                overlaps.append((name1, name2, overlap))
    
    # Calculate totals
    if rois:
        combined = np.any(list(rois.values()), axis=0)
        total_unique = int(np.sum(combined))
    else:
        total_unique = 0
    total_sum = sum(voxel_counts.values())
    
    return {
        'voxel_counts': voxel_counts,
        'overlaps': overlaps,
        'total_unique': total_unique,
        'total_sum': total_sum,
        'has_overlap': len(overlaps) > 0
    }


def get_segment_files(subject: str, session: str, group: str, hemisphere: str = 'left') -> tuple:
    """
    Get segment ROI file paths for a subject based on group.
    
    Parameters
    ----------
    subject : str
        Subject ID (e.g., 'sub-129-921-R1F')
    session : str
        Session ID (e.g., 'ses-260115121231')
    group : str
        Group name: 'CTRL', 'AP', or 'APTR'
    hemisphere : str
        'left' or 'right'
    
    Returns
    -------
    tuple
        (roi_files, roi_names, space_description)
    """
    config = load_config()
    # Use derivatives_flip1 for flip1 pipeline
    base = Path(config['paths']['derivatives_flip1']) / subject / session
    
    if group == 'CTRL':
        # CTRL uses isocortex segments (5 A-P divisions)
        segment_dir = base / 'dwi' / 'isocortex_segments_v2'
        segments = ['anterior', 'anterior_central', 'central', 'posterior_central', 'posterior']
        roi_files = [segment_dir / f'isocortex_{hemisphere}_{seg}.nii.gz' for seg in segments]
        space = 'FOD space (dwi/isocortex_segments_v2/)'
    else:
        # AP/APTR uses graft segments (3 dorsal segments for APTR comparison)
        segment_dir = base / 'dwi' / 'graft_segments_v2' / 'fod_space'
        segments = ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']
        roi_files = [segment_dir / f'graft_{hemisphere}_{seg}.nii.gz' for seg in segments]
        space = 'FOD space (dwi/graft_segments_v2/fod_space/)'
    
    return [str(f) for f in roi_files], segments, space


def check_subject(subject: str, session: str, group: str, verbose: bool = True) -> dict:
    """
    Check ROI overlap for a single subject.
    
    Parameters
    ----------
    subject : str
        Subject ID
    session : str
        Session ID
    group : str
        Group name
    verbose : bool
        Print results to stdout
    
    Returns
    -------
    dict
        Results from check_roi_overlap()
    """
    results = {}
    
    for hemi in ['left', 'right']:
        roi_files, roi_names, space = get_segment_files(subject, session, group, hemi)
        
        if verbose:
            print(f"\n=== {subject} / {session} ({group}) - {hemi.upper()} hemisphere ===")
            print(f"Space: {space}")
        
        result = check_roi_overlap(roi_files, roi_names)
        results[hemi] = result
        
        if verbose:
            print("\nVoxel counts:")
            for name, count in result['voxel_counts'].items():
                print(f"  {name}: {count} voxels")
            
            print(f"\nTotal unique voxels: {result['total_unique']}")
            print(f"Sum of all segments: {result['total_sum']}")
            print(f"Difference (overlap): {result['total_sum'] - result['total_unique']}")
            
            if result['has_overlap']:
                print("\n⚠️  OVERLAPS DETECTED:")
                for roi1, roi2, overlap in result['overlaps']:
                    print(f"  {roi1} & {roi2}: {overlap} voxels")
            else:
                print("\n✅ No overlap detected")
    
    return results


def check_all_subjects(verbose: bool = True) -> dict:
    """
    Check ROI overlap for all subjects in config.
    
    Returns
    -------
    dict
        Dictionary of subject -> results
    """
    config = load_config()
    all_results = {}
    issues = []
    
    for group_name, group_info in config['groups'].items():
        if group_name not in ['CTRL', 'AP', 'APTR']:
            continue
        
        subjects = group_info.get('subjects', [])
        for subj_info in subjects:
            subject = subj_info['id']
            session = subj_info['session']
            
            results = check_subject(subject, session, group_name, verbose=verbose)
            all_results[subject] = results
            
            # Check for issues
            for hemi, result in results.items():
                if result['has_overlap']:
                    issues.append(f"{subject} ({hemi}): overlap detected")
    
    if verbose:
        print("\n" + "="*60)
        if issues:
            print("⚠️  ISSUES FOUND:")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print("✅ All subjects passed overlap check")
    
    return all_results


def main():
    parser = argparse.ArgumentParser(
        description='Check ROI overlap between segmented regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python utils/check_roi_overlap.py --subject sub-129-921-R1F --session ses-260115121231 --group CTRL
  python utils/check_roi_overlap.py --all
        """
    )
    
    parser.add_argument('--subject', type=str, help='Subject ID')
    parser.add_argument('--session', type=str, help='Session ID')
    parser.add_argument('--group', type=str, choices=['CTRL', 'AP', 'APTR'], help='Group name')
    parser.add_argument('--all', action='store_true', help='Check all subjects')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    if args.all:
        check_all_subjects(verbose=not args.quiet)
    elif args.subject and args.session and args.group:
        check_subject(args.subject, args.session, args.group, verbose=not args.quiet)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
