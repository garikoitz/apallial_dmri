#!/usr/bin/env python3
"""
Demo: Convert ROIs and tracts from subject space to group template space.

This script demonstrates the workflow for:
1. Converting ROIs from subject FOD space to template space
2. Converting tracts from subject space to template space

Usage:
    python demo_template_conversion.py --group CTRL --demo

Note: Requires completed population_template with warp files.
"""

import os
import subprocess
import argparse
from pathlib import Path

# Subject info for each group (randomly selected for ROI visualization)
DEMO_SUBJECTS = {
    'CTRL': {
        'id': 'sub-129-1014-L1M',
        'session': 'ses-260115113814',
        'roi_type': 'atlas',  # CTRL uses atlas ROIs
    },
    'AP': {
        'id': 'sub-596-C-OO-NmM',
        'session': 'ses-240911103349',
        'roi_type': 'manual',  # AP uses manual graft ROIs
    },
    'APTR': {
        'id': 'sub-Apallial17NmF',
        'session': 'ses-240816091633',
        'roi_type': 'manual',  # APTR uses manual graft ROIs
    }
}

# All subjects per group
ALL_SUBJECTS = {
    'CTRL': [
        ('sub-129-921-R1F', 'ses-260115121231'),
        ('sub-129-921-R1L1F', 'ses-260116085433'),
        ('sub-129-1014-L1M', 'ses-260115113814'),
        ('sub-138-1010-L1F', 'ses-260115101441'),
    ],
    'AP': [
        ('sub-596-C-OO-NmM', 'ses-240911103349'),
        ('sub-769-C-O5-R1M', 'ses-240911115733'),
        ('sub-769-C-OO-NmM', 'ses-240911111853'),
        ('sub-875-C-OO-NM', 'ses-240911100253'),
    ],
    'APTR': [
        ('sub-AP72-862-901-NmF', 'ses-251210093342'),
        ('sub-AP72-862-901-R1M', 'ses-251210100043'),
        ('sub-Apallial17NmF', 'ses-240816091633'),
        ('sub-Apallial20R1F', 'ses-240909094232'),
    ]
}


def run_cmd(cmd, capture=True, verbose=True):
    """Run shell command."""
    if verbose:
        print(f"  CMD: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=capture, text=True)
    if result.returncode != 0 and capture:
        print(f"  ERROR: {result.stderr}")
    return result.returncode == 0


def get_subject_warp_file(group, subject_id, deriv_dir):
    """Get the warp file for a subject from population_template output."""
    # population_template saves warps with subject names
    # e.g., warps/29-1014-L1M.mif (uses truncated names)
    warp_dir = deriv_dir / "templates" / group / "warps"
    
    # Try different naming patterns
    patterns = [
        f"{subject_id}.mif",
        f"{subject_id.replace('sub-', '')}.mif",
        # population_template truncates names, try partial match
    ]
    
    for pattern in patterns:
        warp_file = warp_dir / pattern
        if warp_file.exists():
            return warp_file
    
    # List available warps and find match
    if warp_dir.exists():
        for f in warp_dir.iterdir():
            if f.suffix == '.mif':
                # Check if subject_id is contained in the filename
                if subject_id.replace('sub-', '').replace('-', '') in f.stem.replace('-', ''):
                    return f
    
    return None


def convert_roi_to_template(roi_path, warp_file, template_path, output_path):
    """
    Convert an ROI from subject space to template space using the warp field.
    
    The warp from population_template goes from TEMPLATE to SUBJECT space.
    To warp ROIs (subject -> template), we need to INVERT the warp.
    
    Steps:
    1. Invert the warp field
    2. Apply inverted warp to ROI
    """
    # Invert warp (template->subject becomes subject->template)
    warp_inv = str(output_path).replace('.nii.gz', '_warp_inv.mif').replace('.mif', '_warp_inv.mif')
    
    cmd = f"warpinvert {warp_file} {warp_inv} -force"
    if not run_cmd(cmd):
        return False
    
    # Apply inverted warp to transform ROI to template space
    cmd = f"mrtransform {roi_path} {output_path} -warp {warp_inv} -template {template_path} -interp nearest -force"
    if not run_cmd(cmd):
        return False
    
    # Clean up temp warp
    if os.path.exists(warp_inv):
        os.remove(warp_inv)
    
    return True


def convert_tract_to_template(tract_path, warp_file, output_path):
    """
    Convert a tract from subject space to template space.
    
    For tracts, we use tcktransform with the inverted warp.
    """
    # Invert warp
    warp_inv = str(output_path).replace('.tck', '_warp_inv.mif')
    
    cmd = f"warpinvert {warp_file} {warp_inv} -force"
    if not run_cmd(cmd):
        return False
    
    # Apply to tract
    cmd = f"tcktransform {tract_path} {warp_inv} {output_path} -force"
    if not run_cmd(cmd):
        return False
    
    # Clean up
    if os.path.exists(warp_inv):
        os.remove(warp_inv)
    
    return True


def demo_conversion(group, deriv_dir):
    """Demonstrate ROI and tract conversion for one subject."""
    
    print(f"\n{'='*60}")
    print(f"DEMO: Converting ROIs and tracts for {group} group")
    print(f"{'='*60}")
    
    # Get demo subject info
    subj_info = DEMO_SUBJECTS[group]
    subj_id = subj_info['id']
    session = subj_info['session']
    
    print(f"\nDemo subject: {subj_id}")
    print(f"Session: {session}")
    
    # Paths
    template_dir = deriv_dir / "templates" / group
    template_fod = template_dir / "wm_fod_template.mif"
    warp_dir = template_dir / "warps"
    
    subj_dir = deriv_dir / subj_id / session / "dwi"
    
    # Check if template exists
    if not template_fod.exists():
        print(f"\n[ERROR] Template not found: {template_fod}")
        print("  Please wait for population_template to complete.")
        return False
    
    # Check if warps exist
    if not warp_dir.exists():
        print(f"\n[ERROR] Warp directory not found: {warp_dir}")
        return False
    
    # Find warp file for this subject
    warp_file = get_subject_warp_file(group, subj_id, deriv_dir)
    if not warp_file:
        print(f"\n[ERROR] No warp file found for {subj_id}")
        print(f"  Available warps in {warp_dir}:")
        for f in warp_dir.iterdir():
            print(f"    {f.name}")
        return False
    
    print(f"\nWarp file: {warp_file}")
    
    # Create output directory
    out_dir = template_dir / "template_space"
    out_dir.mkdir(exist_ok=True)
    
    # Demo 1: Convert ROIs
    print(f"\n--- Converting ROIs to template space ---")
    
    if subj_info['roi_type'] == 'atlas':
        # CTRL: use atlas ROIs
        roi_dir = subj_dir / "atlas" / "rois" / "fod_space"
        rois_to_convert = [
            "roi_isocortex_left.nii.gz",
            "roi_isocortex_right.nii.gz",
            "roi_thalamus_left.nii.gz",
            "roi_thalamus_right.nii.gz",
            "roi_striatum_left.nii.gz",
            "roi_striatum_right.nii.gz",
        ]
    else:
        # AP/APTR: use manual graft ROIs
        roi_dir = subj_dir / "manual_roi" / "fod_space"
        rois_to_convert = [
            "Graft_left.nii.gz",
            "Graft_right.nii.gz",
        ]
    
    for roi_name in rois_to_convert:
        roi_path = roi_dir / roi_name
        if roi_path.exists():
            out_roi = out_dir / f"{subj_id}_{roi_name}"
            print(f"\nConverting: {roi_name}")
            convert_roi_to_template(str(roi_path), str(warp_file), str(template_fod), str(out_roi))
        else:
            print(f"\n[SKIP] ROI not found: {roi_path}")
    
    # Demo 2: Convert one tract
    print(f"\n--- Converting tract to template space ---")
    
    tract_dir = subj_dir / "segment_tracts"
    
    # Find a cleaned tract
    if subj_info['roi_type'] == 'atlas':
        tract_pattern = "tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck"
    else:
        tract_pattern = "tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck"
    
    tract_path = tract_dir / tract_pattern
    
    if tract_path.exists():
        out_tract = out_dir / f"{subj_id}_{tract_pattern}"
        print(f"\nConverting: {tract_pattern}")
        convert_tract_to_template(str(tract_path), str(warp_file), str(out_tract))
    else:
        print(f"\n[SKIP] Tract not found: {tract_path}")
        # Try to find any cleaned tract
        if tract_dir.exists():
            for f in tract_dir.iterdir():
                if 'cleaned' in f.name and f.suffix == '.tck':
                    out_tract = out_dir / f"{subj_id}_{f.name}"
                    print(f"\nConverting alternative: {f.name}")
                    convert_tract_to_template(str(f), str(warp_file), str(out_tract))
                    break
    
    # Visualization command
    print(f"\n--- Visualization command ---")
    print(f"\nmrview {template_fod} \\")
    for roi_name in rois_to_convert[:2]:  # First 2 ROIs
        out_roi = out_dir / f"{subj_id}_{roi_name}"
        if out_roi.exists():
            print(f"  -overlay.load {out_roi} -overlay.opacity 0.5 \\")
    print(f"  -tractography.load {out_dir}/{subj_id}_*_cleaned.tck &")
    
    print(f"\n{'='*60}")
    print("Demo complete!")
    print(f"{'='*60}")
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Convert ROIs and tracts to template space')
    parser.add_argument('--group', type=str, choices=['CTRL', 'AP', 'APTR'], required=True,
                        help='Group to process')
    parser.add_argument('--demo', action='store_true',
                        help='Run demo conversion for one subject')
    parser.add_argument('--all', action='store_true',
                        help='Convert all subjects and tracts')
    
    args = parser.parse_args()
    
    deriv_dir = Path("/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2")
    
    if args.demo:
        demo_conversion(args.group, deriv_dir)
    elif args.all:
        print("Full conversion not yet implemented. Use --demo first.")


if __name__ == '__main__':
    main()
