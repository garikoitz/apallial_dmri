#!/usr/bin/env python3
"""
Generic tract cleaning script using MATLAB AFQ_removeFiberOutliers.

Cleans fiber tracts by removing outliers based on distance and length.
Parameters are read from config.yaml.

Usage:
    python clean_tracts.py                                    # All tracts, all subjects
    python clean_tracts.py --group CTRL                       # CTRL group only
    python clean_tracts.py --subject sub-129-921-R1F          # Specific subject
    python clean_tracts.py --tract Isocortex_Striatum         # Specific tract
    python clean_tracts.py --hemisphere left                  # Only left hemisphere
    python clean_tracts.py --dry-run                          # Preview commands
"""

import argparse
import subprocess
import os
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups,
    get_tract_definitions, get_matlab_cleaning_params
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


def get_tracts(tract_name=None, segment=None):
    """Get tract definitions, optionally filtered.
    
    For segment-based tracts, returns a synthetic tract definition.
    """
    # For segment-based tracts, create a synthetic tract definition
    if segment and tract_name:
        return [{'name': tract_name}]
    
    tract_defs = get_tract_definitions()
    if tract_name:
        for t in tract_defs:
            if t['name'] == tract_name:
                return [t]
        raise ValueError(f"Tract {tract_name} not found")
    return tract_defs


def run_matlab_cleaning(tracks_file, output_file, dry_run=False):
    """Run MATLAB AFQ cleaning on a tract file."""
    paths = get_paths()
    params = get_matlab_cleaning_params()
    
    matlab_bin = paths['matlab_bin']
    afq_path = paths['afq_path']
    
    max_dist = params['max_dist']
    max_len = params['max_len']
    num_nodes = params['num_nodes']
    max_iter = params.get('max_iter', 10)
    
    dwidir = str(tracks_file.parent)
    tckfilename = tracks_file.name
    output_name = output_file.name
    
    # Use -batch flag for true non-interactive mode (faster and more reliable than -r)
    matlab_cmd = (
        f"addpath(genpath('{afq_path}')); "
        f"dwidir='{dwidir}'; "
        f"tckfilename='{tckfilename}'; "
        f"fg=fgRead(fullfile(dwidir,tckfilename)); "
        f"[fg_clean,~,SF]=AFQ_removeFiberOutliers(fg,{max_dist},{max_len},{num_nodes},'mean',0,{max_iter},0); "
        f"AFQ_fgWrite(fg_clean,fullfile(dwidir,'{output_name}'),'tck'); "
        f"fprintf('Cleaned: %d -> %d fibers\\n',length(fg.fibers),length(fg_clean.fibers));"
    )
    
    cmd = [matlab_bin, "-nodisplay", "-nosplash", "-batch", matlab_cmd]
    
    if dry_run:
        print(f"    [DRY RUN] Would run MATLAB cleaning")
        print(f"    Params: maxDist={max_dist}, maxLen={max_len}, numNodes={num_nodes}")
        return True
    
    print(f"    Running MATLAB AFQ cleaning (maxDist={max_dist}, maxLen={max_len})...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    
    if result.returncode != 0:
        print(f"    ERROR: MATLAB failed")
        return False
    
    # Check if output was created
    if output_file.exists():
        # Get fiber count
        info_result = subprocess.run(
            ["tckinfo", str(output_file)],
            capture_output=True, text=True
        )
        for line in info_result.stdout.split('\n'):
            if 'count:' in line and 'total_count' not in line:
                print(f"    ✓ {line.strip()}")
                break
        return True
    else:
        print(f"    ERROR: Output file not created")
        return False


def process_tract(subject_id, session, tract_name, hemisphere, dry_run, segment=None, target_hemisphere=None):
    """Process a single tract for cleaning.
    
    Args:
        subject_id: Subject ID
        session: Session ID
        tract_name: Tract name (e.g., 'Isocortex_Striatum' or 'graft' or 'isocortex')
        hemisphere: Source hemisphere ('left' or 'right')
        dry_run: If True, show commands without executing
        segment: Optional segment name (e.g., 'anterior_dorsal') for segment-based tracts
        target_hemisphere: Target hemisphere for segment tracts (defaults to opposite if None)
    """
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    
    # Determine file path based on whether this is a segment tract
    if segment:
        # Segment-based tract: tracks_{source}_{hemi}_{segment}_to_thalamus_{target_hemi}.tck
        segment_dir = dwi_dir / "segment_tracts"
        target_hemi = target_hemisphere if target_hemisphere else ('left' if hemisphere == 'right' else 'right')
        input_file = segment_dir / f"tracks_{tract_name}_{hemisphere}_{segment}_to_thalamus_{target_hemi}.tck"
        output_file = segment_dir / f"tracks_{tract_name}_{hemisphere}_{segment}_to_thalamus_{target_hemi}_cleaned.tck"
        display_name = f"{tract_name}_{hemisphere}_{segment}_to_thalamus_{target_hemi}"
    else:
        # Standard tract: tracks_{tract_name}_{hemisphere}.tck
        input_file = dwi_dir / f"tracks_{tract_name}_{hemisphere}.tck"
        output_file = dwi_dir / f"tracks_{tract_name}_{hemisphere}_cleaned.tck"
        display_name = f"{tract_name}_{hemisphere}"
    
    print(f"  {display_name}:")
    
    if not input_file.exists():
        print(f"    SKIP: Track file not found: {input_file.name}")
        return None
    
    if output_file.exists() and not dry_run:
        print(f"    SKIP: Already cleaned: {output_file.name}")
        return True
    
    return run_matlab_cleaning(input_file, output_file, dry_run)


def main():
    parser = argparse.ArgumentParser(
        description='Clean fiber tracts using MATLAB AFQ',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--tract', type=str, help='Process specific tract type')
    parser.add_argument('--hemisphere', type=str, choices=['left', 'right', 'both'],
                        default='both', help='Hemisphere to process')
    parser.add_argument('--segment', type=str, help='Segment name for segment-based tracts (e.g., anterior_dorsal)')
    parser.add_argument('--target-hemisphere', type=str, choices=['left', 'right'],
                        help='Target hemisphere for segment tracts (defaults to contralateral)')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    tracts = get_tracts(tract_name=args.tract, segment=args.segment)
    hemispheres = ['left', 'right'] if args.hemisphere == 'both' else [args.hemisphere]
    
    # Show cleaning parameters
    params = get_matlab_cleaning_params()
    print(f"Cleaning parameters from config:")
    print(f"  max_dist: {params['max_dist']}")
    print(f"  max_len: {params['max_len']}")
    print(f"  num_nodes: {params['num_nodes']}")
    print(f"\nProcessing {len(subjects)} subjects, {len(tracts)} tracts")
    
    results = []
    for subject_id, session in subjects:
        print(f"\n{'='*60}")
        print(f"Subject: {subject_id}")
        print(f"{'='*60}")
        
        for tract in tracts:
            for hemi in hemispheres:
                result = process_tract(
                    subject_id, session, tract['name'], hemi, args.dry_run,
                    segment=args.segment, target_hemisphere=args.target_hemisphere
                )
                if result is not None:
                    results.append((subject_id, tract['name'], hemi, result))
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success_count = sum(1 for r in results if r[3])
    print(f"  Cleaned: {success_count}/{len(results)}")


if __name__ == "__main__":
    main()
