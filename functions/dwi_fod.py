#!/usr/bin/env python3
"""
Generate FOD from preprocessed DWI for tractography (flip1 pipeline).

Reads configuration from config1.yaml. Assumes DWI has already been preprocessed
and registered to T2w space.

Usage:
    python dwi_fod.py --subject sub-Apallial17NmF
    python dwi_fod.py --group APTR
    python dwi_fod.py --group APTR --force  # Force regeneration

Expected input:  dwi/dwi_flip_gradcheck_denoise_bias_anat.mif
Expected output: dwi/wm_fod_norm.mif (for tractography)
"""

import argparse
import subprocess
import os
from pathlib import Path
from typing import Optional, Tuple, List

# Force flip1 config
os.environ['MBRAIN_PIPELINE'] = 'flip1'
from config_loader import load_config


def get_subjects(cfg: dict, group: Optional[str] = None, subject: Optional[str] = None) -> List[Tuple[str, str]]:
    """Get list of (subject_id, session) tuples based on filters."""
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


def run_cmd(cmd: str, description: Optional[str] = None) -> bool:
    """Run shell command with optional description."""
    if description:
        print(f"  {description}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:500]}")
        return False
    return True


def process_subject(subject: str, session: str, cfg: dict, force: bool = False) -> bool:
    """
    Generate FOD for flip1 pipeline.
    
    Steps:
      1. Create mask (if not exists)
      2. Response function estimation (dhollander)
      3. FOD estimation (msmt_csd)
      4. FOD normalization (WM + CSF only)
    """
    paths = cfg.get('paths', {})
    deriv = Path(paths['derivatives_dir'])
    
    dwi_dir = deriv / subject / session / "dwi"
    
    print(f"\n{'='*60}")
    print(f"Processing: {subject}/{session}")
    print(f"{'='*60}")
    
    # Input file - final preprocessed DWI
    dwi_input = dwi_dir / "dwi_flip_gradcheck_denoise_bias_anat.mif"
    if not dwi_input.exists():
        print(f"  ERROR: Preprocessed DWI not found: {dwi_input}")
        return False
    
    # Output files
    mask = dwi_dir / "mask.mif"
    resp_wm = dwi_dir / "response_wm.txt"
    resp_gm = dwi_dir / "response_gm.txt"
    resp_csf = dwi_dir / "response_csf.txt"
    resp_voxels = dwi_dir / "response_voxels.mif"
    wm_fod = dwi_dir / "wm_fod.mif"
    gm_fod = dwi_dir / "gm_fod.mif"
    csf_fod = dwi_dir / "csf_fod.mif"
    wm_fod_norm = dwi_dir / "wm_fod_norm.mif"
    csf_fod_norm = dwi_dir / "csf_fod_norm.mif"
    
    # Check if already done
    if wm_fod_norm.exists() and not force:
        print(f"  Skip: FOD already exists (use --force to regenerate)")
        return True
    
    # Step 1: Create mask (if not exists)
    if not mask.exists() or force:
        if not run_cmd(f"dwi2mask '{dwi_input}' '{mask}' -force -quiet",
                       "Step 1: Creating brain mask..."):
            return False
    
    # Step 2: Response functions (dhollander)
    if not run_cmd(f"dwi2response dhollander '{dwi_input}' '{resp_wm}' '{resp_gm}' '{resp_csf}' "
                   f"-voxels '{resp_voxels}' -mask '{mask}' -force -quiet",
                   "Step 2: Estimating response functions (dhollander)..."):
        return False
    
    # Step 3: FOD estimation (msmt_csd)
    if not run_cmd(f"dwi2fod msmt_csd '{dwi_input}' "
                   f"'{resp_wm}' '{wm_fod}' "
                   f"'{resp_gm}' '{gm_fod}' "
                   f"'{resp_csf}' '{csf_fod}' "
                   f"-mask '{mask}' -force -quiet",
                   "Step 3: Computing FOD (msmt_csd)..."):
        return False
    
    # Step 4: FOD normalization (WM + CSF only, GM often fails for mouse brain)
    result = subprocess.run(
        f"mtnormalise '{wm_fod}' '{wm_fod_norm}' '{csf_fod}' '{csf_fod_norm}' -mask '{mask}' -force -quiet",
        shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  WARNING: mtnormalise failed (common for mouse brain): {result.stderr[:200]}")
        print(f"  Using unnormalized WM FOD...")
        run_cmd(f"cp '{wm_fod}' '{wm_fod_norm}'", "Copying WM FOD as normalized...")
    else:
        print(f"  Step 4: FOD normalization (WM + CSF only)... done")
    
    print(f"  ✓ Done! Output: {wm_fod_norm}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Generate FOD from preprocessed DWI (flip1 pipeline)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python dwi_fod.py --subject sub-Apallial17NmF
    python dwi_fod.py --group APTR
    python dwi_fod.py --group APTR --force
        """
    )
    parser.add_argument('--group', type=str,
                        help='Process all subjects in group (CTRL, AP, APTR)')
    parser.add_argument('--subject', type=str,
                        help='Process specific subject')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration even if output exists')
    
    args = parser.parse_args()
    
    # Load config1.yaml
    cfg = load_config()
    
    print(f"Pipeline: flip1")
    print(f"Derivatives: {cfg['paths']['derivatives_dir']}")
    
    if not args.group and not args.subject:
        parser.error("Must specify --group or --subject")
    
    subjects = get_subjects(cfg, group=args.group, subject=args.subject)
    
    print(f"Processing {len(subjects)} subjects...")
    
    results = []
    for subject_id, session in subjects:
        success = process_subject(subject_id, session, cfg, force=args.force)
        results.append((subject_id, success))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for subj, success in results:
        status = "✓" if success else "✗"
        print(f"  {status} {subj}")


if __name__ == "__main__":
    main()
