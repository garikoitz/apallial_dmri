#!/usr/bin/env python3
"""
Generate whole-brain tractograms for mouse brain DWI data.

This script generates dense whole-brain tractograms using tckgen with the iFOD2
algorithm. Unlike ROI-to-ROI tractography, this generates streamlines throughout
the entire brain mask without anatomical priors (no ACT/5TT).

The whole-brain approach:
- Seeds streamlines throughout the brain mask
- No anatomical constraints (mask-only, no 5TT tissue segmentation)
- Generates 1 million streamlines per subject
- Works consistently across all groups (CTRL, AP, APTR)

Usage:
    python functions/generate_wholebrain_tractogram.py --subject sub-129-921-R1F
    python functions/generate_wholebrain_tractogram.py --group CTRL --n-subjects 2
    python functions/generate_wholebrain_tractogram.py --all

Author: MBRAIN Project
Date: 2026-01-31
"""

import argparse
import subprocess
import sys
from pathlib import Path
from config_loader import load_config


def run_command(cmd, description=None):
    """Run a shell command and handle errors."""
    if description:
        print(f"\n{'='*60}")
        print(f"  {description}")
        print(f"{'='*60}")
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: Command failed with code {result.returncode}")
        print(f"STDERR: {result.stderr}")
        return False
    
    if result.stdout:
        print(result.stdout)
    return True


def get_subject_session(config, subject_id):
    """Get the session for a subject from config."""
    for group_name in ['CTRL', 'AP', 'APTR']:
        group = config.get('groups', {}).get(group_name, [])
        for subj in group:
            if subj['id'] == subject_id:
                return subj['session'], group_name
    return None, None


def get_subjects_for_group(config, group_name, n_subjects=None):
    """Get subjects for a specific group."""
    group = config.get('groups', {}).get(group_name, [])
    subjects = [(s['id'], s['session']) for s in group]
    if n_subjects:
        subjects = subjects[:n_subjects]
    return subjects


def generate_wholebrain_tractogram(
    subject_id,
    session,
    base_path,
    config,
    num_streamlines=1000000,
    force=False
):
    """
    Generate whole-brain tractogram for a single subject.
    
    Parameters
    ----------
    subject_id : str
        Subject ID (e.g., 'sub-129-921-R1F')
    session : str
        Session ID (e.g., 'ses-260115121231')
    base_path : Path
        Base path to derivatives/mrtrix_preproc_flip1
    config : dict
        Configuration dictionary
    num_streamlines : int
        Number of streamlines to generate (default: 1 million)
    force : bool
        Overwrite existing files
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    # Build paths
    subj_path = base_path / subject_id / session / "dwi"
    
    # Input files - use group-calibrated FOD for comparable FBA analysis
    fod_file = subj_path / "wm_fod_norm_group.mif"
    mask_file = subj_path / "mask_anat.mif"
    
    # Output files
    output_dir = subj_path / "wholebrain"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    tractogram_file = output_dir / f"wholebrain_{num_streamlines//1000}k.tck"
    
    # Check if output exists
    if tractogram_file.exists() and not force:
        print(f"Tractogram already exists: {tractogram_file}")
        print("Use --force to regenerate")
        return True
    
    # Verify input files exist
    if not fod_file.exists():
        print(f"ERROR: FOD file not found: {fod_file}")
        return False
    if not mask_file.exists():
        print(f"ERROR: Mask file not found: {mask_file}")
        return False
    
    # Get tractography parameters from config
    tract_params = config.get('tractography', {})
    step_size = tract_params.get('step_size', 0.02)
    angle = tract_params.get('angle', 60)
    cutoff = tract_params.get('cutoff', 0.1)
    min_length = tract_params.get('min_length', 2)
    max_length = tract_params.get('max_length', 15)
    
    print(f"\n{'#'*70}")
    print(f"# WHOLE-BRAIN TRACTOGRAPHY: {subject_id}")
    print(f"{'#'*70}")
    print(f"  Session: {session}")
    print(f"  FOD: {fod_file}")
    print(f"  Mask: {mask_file}")
    print(f"  Output: {tractogram_file}")
    print(f"  Streamlines: {num_streamlines:,}")
    print(f"  Parameters:")
    print(f"    - step_size: {step_size} mm")
    print(f"    - angle: {angle}°")
    print(f"    - cutoff: {cutoff}")
    print(f"    - min_length: {min_length} mm")
    print(f"    - max_length: {max_length} mm")
    
    # Build tckgen command
    # Note: No -seed_image needed when using -select with mask
    # The -mask option constrains where streamlines can exist
    tckgen_cmd = [
        "tckgen",
        str(fod_file),
        str(tractogram_file),
        "-algorithm", "iFOD2",
        "-select", str(num_streamlines),
        "-step", str(step_size),
        "-angle", str(angle),
        "-cutoff", str(cutoff),
        "-minlength", str(min_length),
        "-maxlength", str(max_length),
        "-mask", str(mask_file),
        "-seed_image", str(mask_file),  # Seed throughout the mask
        "-nthreads", "8",
        "-force" if force else ""
    ]
    
    # Remove empty strings
    tckgen_cmd = [c for c in tckgen_cmd if c]
    
    success = run_command(
        tckgen_cmd,
        f"Generating {num_streamlines:,} streamlines with tckgen iFOD2"
    )
    
    if not success:
        return False
    
    # Get tractogram info
    tckinfo_cmd = ["tckinfo", str(tractogram_file)]
    run_command(tckinfo_cmd, "Tractogram statistics")
    
    # Get file size
    file_size_mb = tractogram_file.stat().st_size / (1024 * 1024)
    print(f"\nTractogram file size: {file_size_mb:.1f} MB")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate whole-brain tractograms for mouse brain DWI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single subject
  python functions/generate_wholebrain_tractogram.py --subject sub-129-921-R1F

  # First 2 subjects from each group
  python functions/generate_wholebrain_tractogram.py --all --n-subjects 2

  # Specific group
  python functions/generate_wholebrain_tractogram.py --group CTRL --n-subjects 2

  # Custom streamline count
  python functions/generate_wholebrain_tractogram.py --subject sub-129-921-R1F --streamlines 500000
        """
    )
    
    parser.add_argument(
        "--subject", "-s",
        help="Subject ID (e.g., sub-129-921-R1F)"
    )
    parser.add_argument(
        "--group", "-g",
        choices=["CTRL", "AP", "APTR"],
        help="Process subjects from a specific group"
    )
    parser.add_argument(
        "--all", "-a",
        action="store_true",
        help="Process all subjects (use with --n-subjects to limit)"
    )
    parser.add_argument(
        "--n-subjects", "-n",
        type=int,
        help="Number of subjects to process per group"
    )
    parser.add_argument(
        "--streamlines", "-k",
        type=int,
        default=1000000,
        help="Number of streamlines to generate (default: 1000000)"
    )
    parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Overwrite existing tractograms"
    )
    parser.add_argument(
        "--config", "-c",
        default="config1.yaml",
        help="Path to config file (default: config1.yaml)"
    )
    parser.add_argument(
        "--derivatives-dir",
        help="Override derivatives directory (default: from config or MBRAIN_DERIVATIVES_DIR env var)"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Get base path - CLI flag > env var > config
    if args.derivatives_dir:
        base_path = Path(args.derivatives_dir)
    else:
        import os
        env_deriv = os.environ.get('MBRAIN_DERIVATIVES_DIR')
        if env_deriv:
            base_path = Path(env_deriv)
        else:
            paths = config.get('paths', {})
            base_path = Path(paths.get('derivatives_flip1') or paths.get('derivatives_dir'))
            if not base_path.exists():
                repo_dir = Path(paths.get('repo_dir', '/Users/glerma/soft/mbrain'))
                base_path = repo_dir / "data" / "derivatives" / "mrtrix_preproc_flip1"
    
    # Determine which subjects to process
    subjects_to_process = []
    
    if args.subject:
        session, group = get_subject_session(config, args.subject)
        if session:
            subjects_to_process.append((args.subject, session, group))
        else:
            print(f"ERROR: Subject {args.subject} not found in config")
            sys.exit(1)
    
    elif args.group:
        subjects = get_subjects_for_group(config, args.group, args.n_subjects)
        for subj_id, session in subjects:
            subjects_to_process.append((subj_id, session, args.group))
    
    elif args.all:
        for group_name in ['CTRL', 'AP', 'APTR']:
            subjects = get_subjects_for_group(config, group_name, args.n_subjects)
            for subj_id, session in subjects:
                subjects_to_process.append((subj_id, session, group_name))
    
    else:
        parser.print_help()
        print("\nERROR: Please specify --subject, --group, or --all")
        sys.exit(1)
    
    # Process subjects
    print(f"\n{'='*70}")
    print(f"WHOLE-BRAIN TRACTOGRAPHY PIPELINE")
    print(f"{'='*70}")
    print(f"Subjects to process: {len(subjects_to_process)}")
    print(f"Streamlines per subject: {args.streamlines:,}")
    print(f"Base path: {base_path}")
    
    results = []
    for subj_id, session, group in subjects_to_process:
        print(f"\n[{len(results)+1}/{len(subjects_to_process)}] Processing {subj_id} ({group})")
        
        success = generate_wholebrain_tractogram(
            subject_id=subj_id,
            session=session,
            base_path=base_path,
            config=config,
            num_streamlines=args.streamlines,
            force=args.force
        )
        
        results.append((subj_id, group, success))
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    successful = [r for r in results if r[2]]
    failed = [r for r in results if not r[2]]
    
    print(f"Successful: {len(successful)}/{len(results)}")
    for subj_id, group, _ in successful:
        print(f"  ✓ {subj_id} ({group})")
    
    if failed:
        print(f"\nFailed: {len(failed)}/{len(results)}")
        for subj_id, group, _ in failed:
            print(f"  ✗ {subj_id} ({group})")
    
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
