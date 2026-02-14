#!/usr/bin/env python3
"""
Generic preprocessing script for DWI data.

Runs preprocessing steps on subjects from any group.
All parameters are read from config.yaml.

Usage:
    python preprocess_subjects.py                              # All subjects
    python preprocess_subjects.py --group CTRL                 # CTRL group only
    python preprocess_subjects.py --subject sub-129-921-R1F    # Specific subject
    python preprocess_subjects.py --steps denoise unring mask  # Specific steps
    python preprocess_subjects.py --dry-run                    # Preview commands
"""

import argparse
import subprocess
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups
)


PREPROCESSING_STEPS = [
    'convert',      # Convert NIfTI to MIF
    'gradcheck',    # Gradient orientation check
    'denoise',      # Denoising
    'unring',       # Gibbs ringing removal
    'mask',         # Brain mask
    'bias',         # Bias field correction
]


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session) tuples."""
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


def run_step(step, dwi_dir, subject_id, session, dry_run=False):
    """Run a single preprocessing step."""
    cfg = load_config()
    
    # Define input/output patterns for each step
    steps_config = {
        'convert': {
            'input': f"{subject_id}_{session}_dwi.nii.gz",
            'output': f"{subject_id}_{session}_dwi.mif",
            'cmd': lambda i, o: ["mrconvert", str(i), str(o), "-force"]
        },
        'gradcheck': {
            'input': f"{subject_id}_{session}_dwi_flip.mif",
            'output': f"{subject_id}_{session}_dwi_gradchecked.mif",
            'cmd': lambda i, o: ["dwigradcheck", str(i), "-export_grad_mrtrix", 
                                  str(dwi_dir / "raw.b"), "-force"]
        },
        'denoise': {
            'input': f"{subject_id}_{session}_dwi_gradchecked.mif",
            'output': f"{subject_id}_{session}_dwi_denoised.mif",
            'cmd': lambda i, o: ["dwidenoise", str(i), str(o), 
                                  "-noise", str(dwi_dir / "noise.mif"), "-force"]
        },
        'unring': {
            'input': f"{subject_id}_{session}_dwi_denoised.mif",
            'output': f"{subject_id}_{session}_dwi_unringed.mif",
            'cmd': lambda i, o: ["mrdegibbs", str(i), str(o), "-force"]
        },
        'mask': {
            'input': f"{subject_id}_{session}_dwi_unringed.mif",
            'output': "mask.mif",
            'cmd': lambda i, o: ["dwi2mask", str(i), str(o), "-force"]
        },
        'bias': {
            'input': f"{subject_id}_{session}_dwi_unringed.mif",
            'output': f"{subject_id}_{session}_dwi_preprocessed.mif",
            'cmd': lambda i, o: ["dwibiascorrect", "ants", str(i), str(o),
                                  "-mask", str(dwi_dir / "mask.mif"), "-force"]
        },
    }
    
    if step not in steps_config:
        print(f"    Unknown step: {step}")
        return False
    
    step_cfg = steps_config[step]
    input_file = dwi_dir / step_cfg['input']
    output_file = dwi_dir / step_cfg['output']
    
    print(f"  Step: {step}")
    
    if not input_file.exists():
        print(f"    SKIP: Input not found: {input_file.name}")
        return False
    
    if output_file.exists():
        print(f"    SKIP: Output exists: {output_file.name}")
        return True
    
    cmd = step_cfg['cmd'](input_file, output_file)
    
    if dry_run:
        print(f"    [DRY RUN] {' '.join(cmd)}")
        return True
    
    print(f"    Running: {' '.join(cmd[:3])}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[:200]}")
        return False
    
    print(f"    ✓ Complete")
    return True


def process_subject(subject_id, session, steps, dry_run):
    """Process a single subject."""
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    
    print(f"\n{'='*60}")
    print(f"Subject: {subject_id}/{session}")
    print(f"{'='*60}")
    
    if not dwi_dir.exists():
        print(f"  ERROR: DWI directory not found")
        return []
    
    results = []
    for step in steps:
        success = run_step(step, dwi_dir, subject_id, session, dry_run)
        results.append((step, success))
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Preprocess DWI data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Available steps: {', '.join(PREPROCESSING_STEPS)}"
    )
    parser.add_argument('--group', type=str, help='Process specific group')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--steps', nargs='+', choices=PREPROCESSING_STEPS,
                        help='Specific steps to run')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    steps = args.steps or PREPROCESSING_STEPS
    
    print(f"Processing {len(subjects)} subjects")
    print(f"Steps: {', '.join(steps)}")
    
    all_results = []
    for subject_id, session in subjects:
        results = process_subject(subject_id, session, steps, args.dry_run)
        all_results.extend([(subject_id, r[0], r[1]) for r in results])
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success_count = sum(1 for r in all_results if r[2])
    print(f"  Completed: {success_count}/{len(all_results)}")


if __name__ == "__main__":
    main()
