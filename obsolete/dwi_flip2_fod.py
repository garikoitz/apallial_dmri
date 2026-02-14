#!/usr/bin/env python3
"""
Apply flip2 transformation to preprocessed DWI and generate FOD.

Takes the preprocessed DWI from old pipeline (mrtrix_preproc) and applies
flip2 transformation (permutation + flip 1 + flip 2) to match T2w orientation.
Then generates response functions and FOD for tractography.

Usage:
    python dwi_flip2_fod.py --subject sub-129-921-R1F
    python dwi_flip2_fod.py --group CTRL  # All CTRL subjects
"""

import argparse
import subprocess
import os
from pathlib import Path

from config_loader import load_config


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session) tuples based on filters."""
    cfg = load_config()
    groups = cfg.get('groups', {})
    subjects = []
    
    if subject:
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
    
    return subjects


def run_cmd(cmd, description=None):
    """Run shell command with optional description."""
    if description:
        print(f"  {description}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:500]}")
        return False
    return True


def process_subject(subject, session):
    """Apply flip2 to DWI and generate FOD for a single subject."""
    cfg = load_config()
    paths = cfg.get('paths', {})
    
    deriv_old = Path(paths['derivatives_dir_old'])
    deriv_new = Path(paths['derivatives_dir'])
    
    old_dir = deriv_old / subject / session / "dwi"
    new_dir = deriv_new / subject / session / "dwi"
    
    print(f"\n{'='*60}")
    print(f"Processing: {subject}/{session}")
    print(f"{'='*60}")
    
    # Check for source file
    old_dwi = old_dir / "dwi_flip_gradcheck_denoise_bias.mif"
    if not old_dwi.exists():
        print(f"  ERROR: Source DWI not found: {old_dwi}")
        return False
    
    # Check if already done
    fod_file = new_dir / "fod_wm_norm.mif"
    if fod_file.exists():
        print(f"  Skip: FOD already exists")
        return True
    
    new_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Permute axes (2,1,0,3 for 4D)
    dwi_perm = new_dir / "dwi_perm.mif"
    if not run_cmd(f"mrconvert '{old_dwi}' '{dwi_perm}' -axes 2,1,0,3 -strides 1,2,3,4 -force -quiet",
                   "Step 1: Permuting axes..."):
        return False
    
    # Step 2: Flip axis 1
    dwi_flip1 = new_dir / "dwi_flip1.mif"
    if not run_cmd(f"mrtransform '{dwi_perm}' '{dwi_flip1}' -flip 1 -reorient_fod no -force -quiet",
                   "Step 2: Flipping axis 1..."):
        return False
    
    # Step 3: Flip axis 2
    dwi_flip2 = new_dir / "dwi_flip2.mif"
    if not run_cmd(f"mrtransform '{dwi_flip1}' '{dwi_flip2}' -flip 2 -reorient_fod no -force -quiet",
                   "Step 3: Flipping axis 2..."):
        return False
    
    # Clean up intermediates
    dwi_perm.unlink(missing_ok=True)
    dwi_flip1.unlink(missing_ok=True)
    
    # Step 4: Create mask
    mask = new_dir / "mask_flip2.mif"
    if not run_cmd(f"dwi2mask '{dwi_flip2}' '{mask}' -force -quiet",
                   "Step 4: Creating brain mask..."):
        return False
    
    # Step 5: Response functions
    resp_wm = new_dir / "response_wm.txt"
    resp_gm = new_dir / "response_gm.txt"
    resp_csf = new_dir / "response_csf.txt"
    if not run_cmd(f"dwi2response dhollander '{dwi_flip2}' '{resp_wm}' '{resp_gm}' '{resp_csf}' -mask '{mask}' -force -quiet",
                   "Step 5: Estimating response functions..."):
        return False
    
    # Step 6: FOD
    fod_wm = new_dir / "fod_wm.mif"
    fod_gm = new_dir / "fod_gm.mif"
    fod_csf = new_dir / "fod_csf.mif"
    if not run_cmd(f"dwi2fod msmt_csd '{dwi_flip2}' '{resp_wm}' '{fod_wm}' '{resp_gm}' '{fod_gm}' '{resp_csf}' '{fod_csf}' -mask '{mask}' -force -quiet",
                   "Step 6: Computing FOD..."):
        return False
    
    # Step 7: Skip mtnormalise (often fails for mouse brain), just copy
    fod_wm_norm = new_dir / "fod_wm_norm.mif"
    if not run_cmd(f"cp '{fod_wm}' '{fod_wm_norm}'",
                   "Step 7: Using unnormalized WM FOD (mtnormalise often fails for mouse)..."):
        return False
    
    # Step 8: Extract b0
    b0 = new_dir / "b0_flip2.mif"
    if not run_cmd(f"dwiextract '{dwi_flip2}' -bzero - | mrmath - mean '{b0}' -axis 3 -force -quiet",
                   "Step 8: Extracting mean b0..."):
        return False
    
    print(f"  ✓ Done!")
    return True


def main():
    parser = argparse.ArgumentParser(description='Apply flip2 to DWI and generate FOD')
    parser.add_argument('--group', type=str, default='CTRL',
                        help='Process specific group (default: CTRL)')
    parser.add_argument('--subject', type=str,
                        help='Process specific subject')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    
    print(f"Processing {len(subjects)} subjects...")
    
    results = []
    for subject_id, session in subjects:
        success = process_subject(subject_id, session)
        results.append((subject_id, success))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for subj, success in results:
        status = "✓" if success else "✗"
        print(f"  {status} {subj}")


if __name__ == "__main__":
    main()
