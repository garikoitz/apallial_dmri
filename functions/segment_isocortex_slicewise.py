#!/usr/bin/env python3
"""
Segment isocortex ROI into anterior/central/posterior AND dorsal/ventral.

For each slice in the medial-lateral (Z) direction:
- Divide each slice into dorsal (top half) and ventral (bottom half)
- Accumulate across slices within each anterior/central/posterior section

Output: 6 segments per hemisphere:
- {hemi}_anterior_dorsal, {hemi}_anterior_ventral
- {hemi}_central_dorsal, {hemi}_central_ventral  
- {hemi}_posterior_dorsal, {hemi}_posterior_ventral

Usage: 
    python segment_isocortex_slicewise.py --subject sub-129-921-R1F
    python segment_isocortex_slicewise.py --group CTRL  # All CTRL subjects
"""

import nibabel as nib
import numpy as np
import os
import argparse
from pathlib import Path

from config_loader import load_config


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session) tuples based on filters."""
    cfg = load_config()
    groups = cfg.get('groups', {})
    subjects = []
    
    if subject:
        # Find specific subject in any group
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session']))
                    return subjects
        raise ValueError(f"Subject {subject} not found in any group")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found. Available: {list(groups.keys())}")
        for s in groups[group]:
            subjects.append((s['id'], s['session']))
    else:
        # Default to CTRL group for isocortex segmentation
        for s in groups.get('CTRL', []):
            subjects.append((s['id'], s['session']))
    
    return subjects


def segment_isocortex(subject, session):
    """Segment isocortex ROI for a single subject."""
    cfg = load_config()
    paths = cfg.get('paths', {})
    
    # Use flip2 derivatives folder
    base_dir = Path(paths['derivatives_dir']) / subject / session

    # Output directory
    output_dir = base_dir / "dwi" / "isocortex_segments_v2"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Segmenting isocortex: {subject}/{session}")
    print(f"{'='*60}")

    # Process both hemispheres
    for hemi in ["left", "right"]:
        print(f"\n--- Processing {hemi} hemisphere ---")
        
        # Input ROI (from atlas coregistration)
        roi_path = base_dir / "dwi" / "atlas" / "rois" / f"roi_isocortex_{hemi}.nii.gz"
        
        if not roi_path.exists():
            print(f"  WARNING: ROI not found: {roi_path}")
            continue
        
        # Load the ROI
        img = nib.load(str(roi_path))
        data = img.get_fdata()
        affine = img.affine
        
        # Find voxels in the ROI
        coords = np.where(data > 0)
        if len(coords[0]) == 0:
            print(f"  WARNING: No voxels found in {hemi} isocortex ROI")
            continue
        
        x_coords = coords[0]  # Anterior-Posterior (axis 0)
        y_coords = coords[1]  # Dorsal-Ventral (axis 1)
        z_coords = coords[2]  # Medial-Lateral left-right (axis 2)
        
        print(f"  Total voxels: {len(x_coords)}")
        print(f"  X range (ant-post): {x_coords.min()} - {x_coords.max()}")
        print(f"  Y range (dor-vent): {y_coords.min()} - {y_coords.max()}")
        print(f"  Z range (med-lat): {z_coords.min()} - {z_coords.max()}")
        
        # Divide X (anterior-posterior) into thirds
        x_min, x_max = x_coords.min(), x_coords.max()
        x_range = x_max - x_min + 1
        x_third = x_range / 3
        
        x_ant_max = x_min + x_third
        x_cent_max = x_min + 2 * x_third
        
        print(f"  Anterior: X >= {x_cent_max:.1f}")
        print(f"  Central: {x_ant_max:.1f} <= X < {x_cent_max:.1f}")
        print(f"  Posterior: X < {x_ant_max:.1f}")
        
        # Initialize output arrays
        segments = {
            "anterior_dorsal": np.zeros_like(data),
            "anterior_ventral": np.zeros_like(data),
            "central_dorsal": np.zeros_like(data),
            "central_ventral": np.zeros_like(data),
            "posterior_dorsal": np.zeros_like(data),
            "posterior_ventral": np.zeros_like(data),
        }
        
        # Process slice by slice in Z (medial-lateral direction)
        unique_z = np.unique(z_coords)
        
        for z in unique_z:
            # Get voxels in this Z slice
            z_mask = (z_coords == z)
            slice_x = x_coords[z_mask]
            slice_y = y_coords[z_mask]
            
            if len(slice_y) == 0:
                continue
            
            # Find dorsal/ventral midpoint for this slice
            y_mid = (slice_y.min() + slice_y.max()) / 2
            
            # Assign each voxel to appropriate segment
            for i, (x, y) in enumerate(zip(slice_x, slice_y)):
                # Determine anterior/central/posterior
                # In flip2: Low X = posterior, High X = anterior
                if x < x_ant_max:
                    ap_section = "posterior"  # Low X = posterior
                elif x < x_cent_max:
                    ap_section = "central"
                else:
                    ap_section = "anterior"  # High X = anterior
                
                # Determine dorsal/ventral (Y axis: higher = more ventral)
                if y >= y_mid:
                    dv_section = "ventral"  # High Y = ventral
                else:
                    dv_section = "dorsal"  # Low Y = dorsal
                
                segment_name = f"{ap_section}_{dv_section}"
                segments[segment_name][int(x), int(y), int(z)] = 1
        
        # Save each segment
        print(f"\n  Segment voxel counts for {hemi}:")
        for seg_name, seg_data in segments.items():
            count = int(np.sum(seg_data))
            print(f"    {seg_name}: {count} voxels")
            
            if count > 0:
                out_path = output_dir / f"isocortex_{hemi}_{seg_name}.nii.gz"
                out_img = nib.Nifti1Image(seg_data.astype(np.float32), affine)
                nib.save(out_img, str(out_path))
    
    print(f"\n  ✓ Done! Segments saved to: {output_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Segment isocortex ROI for CTRL subjects',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, default='CTRL', 
                        help='Process specific group (default: CTRL)')
    parser.add_argument('--subject', type=str, 
                        help='Process specific subject')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    
    print(f"Processing {len(subjects)} subjects...")
    
    results = []
    for subject_id, session in subjects:
        success = segment_isocortex(subject_id, session)
        results.append((subject_id, success))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for subj, success in results:
        status = "✓" if success else "✗"
        print(f"  {status} {subj}")


if __name__ == "__main__":
    main()
