#!/usr/bin/env python3
"""
Generic ROI splitting script.

Splits ROIs into left/right hemispheres based on midline.
Supports all subjects and ROIs defined in config.

Usage:
    python split_rois.py                                    # All ROIs, all subjects
    python split_rois.py --group CTRL                       # CTRL group only
    python split_rois.py --subject sub-129-921-R1F          # Specific subject
    python split_rois.py --roi isocortex                    # Specific ROI
    python split_rois.py --resample-to-fod                  # Also resample to FOD space
    python split_rois.py --dry-run                          # Preview commands
"""

import argparse
import subprocess
import numpy as np
import nibabel as nib
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups, get_roi_definitions
)


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session) tuples based on filters."""
    groups = get_subject_groups()
    subjects = []
    
    if subject:
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session']))
                    return subjects
        raise ValueError(f"Subject {subject} not found")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found")
        for s in groups[group]:
            subjects.append((s['id'], s['session']))
    else:
        for g, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session']))
    
    return subjects


def get_rois(roi_name=None):
    """Get ROI definitions, optionally filtered."""
    roi_defs = get_roi_definitions()
    if roi_name:
        if roi_name in roi_defs:
            return {roi_name: roi_defs[roi_name]}
        raise ValueError(f"ROI {roi_name} not found")
    return roi_defs


def split_roi_hemispheres(input_file, output_left, output_right):
    """Split a ROI NIfTI into left and right hemispheres."""
    img = nib.load(input_file)
    data = img.get_fdata()
    affine = img.affine
    
    # Get midline (assume x is left-right axis)
    midline_x = data.shape[0] // 2
    
    # Create hemisphere masks
    data_left = data.copy()
    data_left[midline_x:, :, :] = 0
    
    data_right = data.copy()
    data_right[:midline_x, :, :] = 0
    
    # Save
    nib.save(nib.Nifti1Image(data_left.astype(np.uint8), affine), output_left)
    nib.save(nib.Nifti1Image(data_right.astype(np.uint8), affine), output_right)
    
    return data_left.sum() > 0, data_right.sum() > 0


def resample_to_fod_space(input_file, reference_file, output_file, mask_file=None):
    """Resample ROI to FOD space using mrtransform."""
    cmd = [
        "mrtransform",
        str(input_file),
        str(output_file),
        "-template", str(reference_file),
        "-interp", "nearest",
        "-force"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        return False
    
    # Apply mask if provided
    if mask_file and Path(mask_file).exists():
        masked_output = output_file.replace('.nii.gz', '_masked.nii.gz')
        cmd = ["mrcalc", str(output_file), str(mask_file), "-mult", str(masked_output), "-force"]
        subprocess.run(cmd, check=True)
    
    return True


def process_subject(subject_id, session, rois, resample_to_fod, dry_run):
    """Process ROI splitting for a single subject."""
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    
    roi_dir = dwi_dir / "atlas" / "rois"
    fod_space_dir = roi_dir / "fod_space"
    fod_file = dwi_dir / "fod_wm.mif"
    mask_file = dwi_dir / "mask.mif"
    
    print(f"\n  Subject: {subject_id}")
    
    if not roi_dir.exists():
        print(f"    SKIP: ROI directory not found")
        return []
    
    if dry_run:
        print(f"    [DRY RUN] Would split {len(rois)} ROIs")
        return [(roi, True) for roi in rois]
    
    results = []
    for roi_name in rois:
        input_file = roi_dir / f"roi_{roi_name}.nii.gz"
        output_left = roi_dir / f"roi_{roi_name}_left.nii.gz"
        output_right = roi_dir / f"roi_{roi_name}_right.nii.gz"
        
        if not input_file.exists():
            print(f"    {roi_name}: SKIP (not found)")
            results.append((roi_name, None))
            continue
        
        has_left, has_right = split_roi_hemispheres(input_file, output_left, output_right)
        print(f"    {roi_name}: L={has_left}, R={has_right}")
        
        # Resample to FOD space if requested
        if resample_to_fod and fod_file.exists():
            fod_space_dir.mkdir(parents=True, exist_ok=True)
            
            for hemi, src in [('left', output_left), ('right', output_right)]:
                if src.exists():
                    dst = fod_space_dir / f"roi_{roi_name}_{hemi}.nii.gz"
                    resample_to_fod_space(src, fod_file, dst, mask_file)
        
        results.append((roi_name, True))
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Split ROIs into hemispheres',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--roi', type=str, help='Process specific ROI')
    parser.add_argument('--resample-to-fod', action='store_true', 
                        help='Also resample ROIs to FOD space')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    rois = get_rois(roi_name=args.roi)
    
    print(f"Processing {len(subjects)} subjects, {len(rois)} ROIs")
    print(f"Resample to FOD space: {args.resample_to_fod}")
    
    all_results = []
    for subject_id, session in subjects:
        results = process_subject(
            subject_id, session, rois, args.resample_to_fod, args.dry_run
        )
        all_results.extend([(subject_id, r[0], r[1]) for r in results])
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success_count = sum(1 for r in all_results if r[2])
    print(f"  Processed: {success_count}/{len(all_results)}")


if __name__ == "__main__":
    main()
