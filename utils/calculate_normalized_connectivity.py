#!/usr/bin/env python3
"""
Calculate normalized connectivity (connectivity density) to control for ROI size bias.

Literature basis:
- Smith et al., 2015 (SIFT2): Recommends SIFT2 weights but acknowledges ROI size effects remain
- Yeh et al., 2016: Recommend normalizing by ROI volume for cross-region comparisons
- van den Heuvel & Sporns, 2011: Use connection density for connectome analysis

Formula: Connectivity Density = SIFT2_sum / seed_voxels (× 1000 for readability)
Units: SIFT2 weight per 1000 seed voxels
"""

import os
import sys
import numpy as np
import nibabel as nib
from pathlib import Path

# Add parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config

# Use config1.yaml
os.environ['MBRAIN_CONFIG'] = str(Path(__file__).parent.parent / 'config1.yaml')

def count_roi_voxels(roi_path):
    """Count non-zero voxels in an ROI."""
    img = nib.load(roi_path)
    data = img.get_fdata()
    return int(np.sum(data > 0))

def load_sift2_connectivity(sift2_csv_path):
    """Load SIFT2 connectivity from CSV file."""
    connectivity = {}
    with open(sift2_csv_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    name = parts[0]
                    value = float(parts[1])
                    connectivity[name] = value
    return connectivity

def main():
    config = load_config(pipeline='flip1')
    base_path = Path(config['paths']['derivatives_flip1'])
    
    segments = ['anterior', 'anterior_central', 'central', 'posterior_central', 'posterior']
    
    print("=" * 100)
    print("ROI SIZE NORMALIZATION ANALYSIS")
    print("Formula: Connectivity Density = SIFT2_sum / seed_voxels × 1000")
    print("=" * 100)
    
    # Process CTRL subjects
    print("\n" + "=" * 100)
    print("CTRL GROUP (Isocortex segments)")
    print("=" * 100)
    
    ctrl_subjects = config['groups']['CTRL']
    
    for subj_info in ctrl_subjects:
        subject_id = subj_info['id']
        session_id = subj_info['session']
        
        subj_path = base_path / subject_id / session_id
        seg_path = subj_path / 'dwi' / 'isocortex_segments_v2'  # Segments are directly here, not in fod_space
        tracts_path = subj_path / 'dwi' / 'segment_tracts'
        
        print(f"\n{subject_id}:")
        print(f"  {'Segment':25s} {'Raw SIFT2':>12s} {'Voxels':>8s} {'Density':>12s} {'Rank(raw)':>10s} {'Rank(norm)':>10s}")
        print(f"  " + "-" * 80)
        
        raw_values = {}
        norm_values = {}
        
        for seg in segments:
            # Load left hemisphere segment
            seg_file = seg_path / f'isocortex_left_{seg}.nii.gz'
            
            if not seg_file.exists():
                print(f"  {seg:25s} {'N/A':>12s}")
                continue
            
            voxels = count_roi_voxels(seg_file)
            
            # Load SIFT2 connectivity from weights file
            tract_name = f'isocortex_left_{seg}_to_thalamus_left'
            weights_file = tracts_path / f'tracks_{tract_name}_cleaned_sift2_weights.txt'
            
            if weights_file.exists():
                weights = np.loadtxt(weights_file)
                raw = float(np.sum(weights))
            else:
                print(f"  {seg:25s} {'No weights file':>12s}")
                raw = 0
                continue
            
            norm = (raw / voxels * 1000) if voxels > 0 else 0
            raw_values[seg] = raw
            norm_values[seg] = norm
        
        # Calculate ranks
        if raw_values:
            raw_sorted = sorted(raw_values.items(), key=lambda x: x[1], reverse=True)
            norm_sorted = sorted(norm_values.items(), key=lambda x: x[1], reverse=True)
            
            raw_ranks = {seg: i+1 for i, (seg, _) in enumerate(raw_sorted)}
            norm_ranks = {seg: i+1 for i, (seg, _) in enumerate(norm_sorted)}
            
            for seg in segments:
                if seg in raw_values:
                    raw = raw_values[seg]
                    voxels = count_roi_voxels(seg_path / f'isocortex_left_{seg}.nii.gz')
                    norm = norm_values[seg]
                    rr = raw_ranks[seg]
                    nr = norm_ranks[seg]
                    rank_change = "↑" if nr < rr else ("↓" if nr > rr else "=")
                    print(f"  {seg:25s} {raw:12.2f} {voxels:8d} {norm:12.4f} {rr:10d} {nr:>9d}{rank_change}")
    
    # Process APTR subjects
    print("\n" + "=" * 100)
    print("APTR GROUP (Graft dorsal segments)")
    print("=" * 100)
    
    aptr_subjects = config['groups']['APTR']
    aptr_segments = ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']
    
    for subj_info in aptr_subjects:
        subject_id = subj_info['id']
        session_id = subj_info['session']
        
        subj_path = base_path / subject_id / session_id
        seg_path = subj_path / 'dwi' / 'graft_segments_fod'
        tracts_path = subj_path / 'dwi' / 'segment_tracts'
        
        print(f"\n{subject_id}:")
        print(f"  {'Segment':25s} {'Raw SIFT2':>12s} {'Voxels':>8s} {'Density':>12s} {'Rank(raw)':>10s} {'Rank(norm)':>10s}")
        print(f"  " + "-" * 80)
        
        raw_values = {}
        norm_values = {}
        voxel_counts = {}
        
        for seg in aptr_segments:
            # APTR segment files are .mif format: seg9_graft_left_{seg}_fod.mif
            seg_file = seg_path / f'seg9_graft_left_{seg}_fod.mif'
            
            if not seg_file.exists():
                continue
            
            # Use mrstats to count non-zero voxels in .mif file
            import subprocess
            result = subprocess.run(
                ['mrstats', '-mask', str(seg_file), '-output', 'count', str(seg_file)],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                voxels = int(float(result.stdout.strip()))
            else:
                voxels = 0
            
            voxel_counts[seg] = voxels
            
            # APTR weights files are named: {seg}_sift2.txt
            weights_file = tracts_path / f'{seg}_sift2.txt'
            
            if weights_file.exists():
                weights = np.loadtxt(weights_file)
                raw = float(np.sum(weights))
            else:
                raw = 0
            
            norm = (raw / voxels * 1000) if voxels > 0 else 0
            raw_values[seg] = raw
            norm_values[seg] = norm
        
        # Calculate ranks
        if raw_values:
            raw_sorted = sorted(raw_values.items(), key=lambda x: x[1], reverse=True)
            norm_sorted = sorted(norm_values.items(), key=lambda x: x[1], reverse=True)
            
            raw_ranks = {seg: i+1 for i, (seg, _) in enumerate(raw_sorted)}
            norm_ranks = {seg: i+1 for i, (seg, _) in enumerate(norm_sorted)}
            
            for seg in aptr_segments:
                if seg in raw_values:
                    raw = raw_values[seg]
                    voxels = voxel_counts[seg]
                    norm = norm_values[seg]
                    rr = raw_ranks[seg]
                    nr = norm_ranks[seg]
                    rank_change = "↑" if nr < rr else ("↓" if nr > rr else "=")
                    print(f"  {seg:25s} {raw:12.2f} {voxels:8d} {norm:12.4f} {rr:10d} {nr:>9d}{rank_change}")
    
    print("\n" + "=" * 100)
    print("INTERPRETATION")
    print("=" * 100)
    print("""
    Connectivity Density = SIFT2 weight per 1000 voxels
    
    ↑ = Rank improved after normalization (smaller ROI, high connectivity per voxel)
    ↓ = Rank decreased after normalization (larger ROI, lower connectivity per voxel)
    = = Rank unchanged
    
    If ranks change significantly, ROI size was biasing the raw connectivity values.
    Use Connectivity Density for fair comparisons between regions of different sizes.
    """)

if __name__ == '__main__':
    main()
