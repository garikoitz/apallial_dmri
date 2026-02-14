#!/usr/bin/env python3
"""
Generic tract extraction script.

Extracts top N tracts by SIFT2 weight from tract files.
Supports all subjects, tracts, and configurable N.

Usage:
    python extract_top_tracts.py                                    # All tracts, all subjects
    python extract_top_tracts.py --group CTRL                       # CTRL group only
    python extract_top_tracts.py --subject sub-129-921-R1F          # Specific subject
    python extract_top_tracts.py --tract Isocortex_Striatum         # Specific tract
    python extract_top_tracts.py --top-n 100                        # Extract top 100 (overrides config)
    python extract_top_tracts.py --cleaned                          # Use cleaned tracks
    python extract_top_tracts.py --dry-run                          # Preview commands
"""

import argparse
import subprocess
import numpy as np
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups, get_tract_definitions
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


def get_tracts(tract_name=None):
    """Get tract definitions, optionally filtered."""
    tract_defs = get_tract_definitions()
    if tract_name:
        for t in tract_defs:
            if t['name'] == tract_name:
                return [t]
        raise ValueError(f"Tract {tract_name} not found")
    return tract_defs


def run_sift2(tracks_file, fod_file, output_weights):
    """Run tcksift2 to compute track weights."""
    cmd = [
        "tcksift2",
        str(tracks_file),
        str(fod_file),
        str(output_weights),
        "-force"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def extract_top_n(tracks_file, weights_file, output_file, top_n):
    """Extract top N tracks by SIFT2 weight."""
    # Load weights
    weights = np.loadtxt(weights_file)
    
    if len(weights) == 0:
        print(f"    ERROR: No weights found")
        return False, 0
    
    # Get threshold for top N
    n_select = min(top_n, len(weights))
    sorted_weights = np.sort(weights)[::-1]
    threshold = sorted_weights[n_select - 1]
    
    # Use tckedit with weight threshold
    cmd = [
        "tckedit",
        str(tracks_file),
        str(output_file),
        "-tck_weights_in", str(weights_file),
        "-minweight", str(threshold),
        "-force"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    ERROR: tckedit failed: {result.stderr[:200]}")
        return False, 0
    
    # Get actual count
    info_result = subprocess.run(
        ["tckinfo", str(output_file)],
        capture_output=True, text=True
    )
    count = 0
    for line in info_result.stdout.split('\n'):
        if 'count:' in line and 'total_count' not in line:
            count = int(line.split(':')[1].strip())
            break
    
    return True, count


def process_tract(subject_id, session, tract_name, hemisphere, top_n, use_cleaned, dry_run):
    """Process a single tract for a subject."""
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    
    # Build filenames based on whether using cleaned tracks
    suffix = "_cleaned" if use_cleaned else ""
    base_name = f"tracks_{tract_name}_{hemisphere}{suffix}"
    
    tracks_file = dwi_dir / f"{base_name}.tck"
    weights_file = dwi_dir / f"{base_name}_sift2_weights.txt"
    output_file = dwi_dir / f"{base_name}_top{top_n}.tck"
    fod_file = dwi_dir / "fod_wm.mif"
    
    print(f"  {tract_name}_{hemisphere}{suffix}:")
    
    if not tracks_file.exists():
        print(f"    SKIP: Track file not found: {tracks_file.name}")
        return None
    
    if dry_run:
        print(f"    [DRY RUN] Would extract top {top_n} from {tracks_file.name}")
        return True
    
    # Run SIFT2 if weights don't exist
    if not weights_file.exists():
        print(f"    Running SIFT2...")
        if not fod_file.exists():
            print(f"    ERROR: FOD file not found")
            return False
        if not run_sift2(tracks_file, fod_file, weights_file):
            print(f"    ERROR: SIFT2 failed")
            return False
    
    # Extract top N
    success, count = extract_top_n(tracks_file, weights_file, output_file, top_n)
    
    if success:
        print(f"    ✓ Extracted {count} tracks -> {output_file.name}")
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description='Extract top N tracts by SIFT2 weight',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--tract', type=str, help='Process specific tract type')
    parser.add_argument('--hemisphere', type=str, choices=['left', 'right', 'both'],
                        default='both', help='Hemisphere to process')
    parser.add_argument('--top-n', type=int, help='Number of top tracks to extract (overrides config)')
    parser.add_argument('--cleaned', action='store_true', help='Use cleaned track files')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    # Get top_n from config if not specified
    cfg = load_config()
    top_n = args.top_n or cfg.get('top_tracts', {}).get('count', 200)
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    tracts = get_tracts(tract_name=args.tract)
    hemispheres = ['left', 'right'] if args.hemisphere == 'both' else [args.hemisphere]
    
    print(f"Extracting top {top_n} tracts")
    print(f"Subjects: {len(subjects)}, Tracts: {len(tracts)}, Hemispheres: {len(hemispheres)}")
    print(f"Using cleaned tracks: {args.cleaned}")
    
    results = []
    for subject_id, session in subjects:
        print(f"\n{'='*60}")
        print(f"Subject: {subject_id}")
        print(f"{'='*60}")
        
        for tract in tracts:
            for hemi in hemispheres:
                result = process_tract(
                    subject_id, session, tract['name'], hemi,
                    top_n, args.cleaned, args.dry_run
                )
                if result is not None:
                    results.append((subject_id, tract['name'], hemi, result))
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success_count = sum(1 for r in results if r[3])
    print(f"  Processed: {success_count}/{len(results)}")


if __name__ == "__main__":
    main()
