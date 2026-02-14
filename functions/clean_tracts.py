#!/usr/bin/env python3
"""
Generic tract cleaning script using MATLAB AFQ_removeFiberOutliers.

Cleans fiber tracts by removing outliers based on distance and length.
Parameters are read from config1.yaml.

Usage:
    python clean_tracts.py                                    # All tracts, all subjects
    python clean_tracts.py --group CTRL                       # CTRL group only
    python clean_tracts.py --subject sub-129-921-R1F          # Specific subject
    python clean_tracts.py --tract Isocortex_Striatum         # Specific tract
    python clean_tracts.py --hemisphere left                  # Only left hemisphere
    python clean_tracts.py --dry-run                          # Preview commands
    
    # Passthrough tract cleaning (batch MATLAB per subject):
    python clean_tracts.py --passthrough                      # Clean all passthrough tracts
    python clean_tracts.py --passthrough --shrink             # Clean shrink passthrough (AP/APTR)
    python clean_tracts.py --passthrough --subject sub-129-921-R1F
    python clean_tracts.py --passthrough --force              # Re-clean even if _cleaned exists
"""

import argparse
import subprocess
import os
import time
from pathlib import Path
from datetime import datetime

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


def run_batch_matlab_cleaning(tract_dir, tck_files, params, paths, dry_run=False, force=False):
    """Run MATLAB AFQ cleaning on multiple tracts in a single MATLAB session.
    
    This is much more efficient than one MATLAB invocation per tract,
    as MATLAB startup takes ~10s each time.
    
    Args:
        tract_dir: Path to directory containing .tck files
        tck_files: List of .tck filenames to clean
        params: MATLAB cleaning parameters dict
        paths: Config paths dict
        dry_run: Preview mode
        force: Re-clean even if _cleaned exists
    
    Returns:
        dict mapping filename -> (success, n_orig, n_cleaned)
    """
    matlab_bin = paths['matlab_bin']
    afq_path = paths['afq_path']
    max_dist = params['max_dist']
    max_len = params['max_len']
    num_nodes = params['num_nodes']
    max_iter = params.get('max_iter', 5)
    min_streamlines = params.get('min_streamlines', 20)
    
    # Filter out already-cleaned files unless --force
    to_clean = []
    skipped = []
    for tck in tck_files:
        cleaned_name = tck.replace('.tck', '_cleaned.tck')
        cleaned_path = tract_dir / cleaned_name
        if cleaned_path.exists() and not force:
            skipped.append(tck)
        else:
            to_clean.append(tck)
    
    if skipped:
        print(f"    Skipping {len(skipped)} already-cleaned tracts")
    
    if not to_clean:
        print(f"    Nothing to clean (all done)")
        return {t: (True, 0, 0) for t in tck_files}
    
    if dry_run:
        print(f"    [DRY RUN] Would clean {len(to_clean)} tracts in single MATLAB session:")
        for t in to_clean:
            print(f"      - {t}")
        return {}
    
    # Build a single MATLAB command that processes all tracts in one session
    # Use a cell array of filenames and loop through them
    file_list = "', '".join(to_clean)
    
    matlab_cmd = (
        f"addpath(genpath('{afq_path}')); "
        f"dwidir = '{tract_dir}'; "
        f"files = {{'{file_list}'}}; "
        f"for i = 1:length(files), "
        f"  tckfile = files{{i}}; "
        f"  outfile = strrep(tckfile, '.tck', '_cleaned.tck'); "
        f"  sffile = strrep(tckfile, '.tck', '_SF.tck'); "
        f"  try, "
        f"    fg = fgRead(fullfile(dwidir, tckfile)); "
        f"    nOrig = length(fg.fibers); "
        f"    if nOrig < {min_streamlines}, "
        f"      fprintf('SKIP_LOW: %s has %d fibers (min={min_streamlines})\\n', tckfile, nOrig); "
        f"      continue; "
        f"    end; "
        f"    lens = cellfun(@(f) size(f,2), fg.fibers); "
        f"    fg.fibers = fg.fibers(lens > 3); "
        f"    nAfterFilter = length(fg.fibers); "
        f"    if nAfterFilter < {min_streamlines}, "
        f"      fprintf('SKIP_LOW: %s has %d fibers after short-filter (min={min_streamlines})\\n', tckfile, nAfterFilter); "
        f"      continue; "
        f"    end; "
        f"    [fg_clean, ~, SF] = AFQ_removeFiberOutliers(fg, {max_dist}, {max_len}, {num_nodes}, 'mean', 0, {max_iter}, 0); "
        f"    nClean = length(fg_clean.fibers); "
        f"    fg_SF = dtiNewFiberGroup([fg.name '_SF'], [0 0 155], [], [], SF.fibers); "
        f"    AFQ_fgWrite(fg_clean, fullfile(dwidir, outfile), 'tck'); "
        f"    AFQ_fgWrite(fg_SF, fullfile(dwidir, sffile), 'tck'); "
        f"    fprintf('CLEANED: %s %d -> %d fibers\\n', tckfile, nOrig, nClean); "
        f"  catch ME, "
        f"    fprintf('ERROR: %s: %s\\n', tckfile, ME.message); "
        f"  end; "
        f"end; "
        f"fprintf('BATCH_DONE\\n'); "
        f"exit;"
    )
    
    cmd = [matlab_bin, "-nodisplay", "-nosplash", "-batch", matlab_cmd]
    
    print(f"    Running MATLAB batch cleaning for {len(to_clean)} tracts...")
    t0 = time.time()
    
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=600 * len(to_clean)  # 10min per tract max
        )
    except subprocess.TimeoutExpired:
        print(f"    ERROR: MATLAB timed out after {600 * len(to_clean)}s")
        return {}
    
    elapsed = time.time() - t0
    
    # Parse MATLAB output
    results = {}
    for line in result.stdout.split('\n'):
        line = line.strip()
        if line.startswith('CLEANED:'):
            parts = line.split()
            fname = parts[1]
            n_orig = int(parts[2])
            n_clean = int(parts[4])
            results[fname] = (True, n_orig, n_clean)
            print(f"      ✓ {fname}: {n_orig} → {n_clean} fibers")
        elif line.startswith('SKIP_LOW:'):
            fname = line.split()[1]
            results[fname] = (False, 0, 0)
            print(f"      ⊘ {line.replace('SKIP_LOW: ', '')}")
        elif line.startswith('ERROR:'):
            fname = line.split()[1].rstrip(':')
            results[fname] = (False, 0, 0)
            print(f"      ✗ {line}")
    
    if result.returncode != 0 and 'BATCH_DONE' not in result.stdout:
        print(f"    WARNING: MATLAB exited with code {result.returncode}")
        if result.stderr:
            # Only print last few lines of stderr
            stderr_lines = result.stderr.strip().split('\n')
            for l in stderr_lines[-5:]:
                print(f"      stderr: {l}")
    
    print(f"    Batch complete in {elapsed:.1f}s ({len(results)}/{len(to_clean)} processed)")
    
    # Add skipped files as success
    for t in skipped:
        results[t] = (True, 0, 0)
    
    return results


def run_passthrough_cleaning(group=None, subject=None, dry_run=False,
                              use_shrink=False, force=False):
    """Clean all passthrough tracts for specified subjects.
    
    Uses batch MATLAB processing (one session per subject directory)
    for efficiency.
    
    Args:
        group: Optional group filter ('CTRL', 'AP', 'APTR')
        subject: Optional subject filter
        dry_run: Preview mode
        use_shrink: Use passthrough_tracts_shrink/ for AP/APTR
        force: Re-clean even if _cleaned exists
    """
    cfg = load_config()
    paths_cfg = get_paths()
    params = get_matlab_cleaning_params()
    deriv_dir = Path(paths_cfg['derivatives_dir'])
    groups = get_subject_groups()
    
    print(f"{'='*70}")
    print(f"PASSTHROUGH TRACT CLEANING")
    print(f"{'='*70}")
    print(f"Parameters: max_dist={params['max_dist']}, max_len={params['max_len']}, "
          f"num_nodes={params['num_nodes']}, max_iter={params.get('max_iter', 5)}")
    print(f"min_streamlines={params.get('min_streamlines', 20)}")
    print(f"Shrink mode: {use_shrink}")
    print(f"Force re-clean: {force}")
    print()
    
    # Determine which subjects to process
    subjects_to_process = []
    for grp_name, subj_list in groups.items():
        if group and grp_name != group:
            continue
        for s in subj_list:
            if subject and s['id'] != subject:
                continue
            subjects_to_process.append((grp_name, s['id'], s['session']))
    
    total_cleaned = 0
    total_skipped_low = 0
    total_errors = 0
    total_tracts = 0
    t_global = time.time()
    
    for grp_name, subj_id, session in subjects_to_process:
        print(f"\n{'='*60}")
        print(f"[{grp_name}] {subj_id}")
        print(f"{'='*60}")
        
        wholebrain_dir = deriv_dir / subj_id / session / "dwi" / "wholebrain"
        
        # Determine which passthrough dirs to clean
        dirs_to_clean = []
        
        # Always clean regular passthrough_tracts
        pt_dir = wholebrain_dir / "passthrough_tracts"
        if pt_dir.exists():
            dirs_to_clean.append(("passthrough_tracts", pt_dir))
        else:
            print(f"  WARNING: {pt_dir} does not exist, skipping")
        
        # For AP/APTR with --shrink, also clean passthrough_tracts_shrink
        if use_shrink and grp_name in ('AP', 'APTR'):
            pt_shrink_dir = wholebrain_dir / "passthrough_tracts_shrink"
            if pt_shrink_dir.exists():
                dirs_to_clean.append(("passthrough_tracts_shrink", pt_shrink_dir))
            else:
                print(f"  WARNING: {pt_shrink_dir} does not exist, skipping")
        
        for dir_label, pt_path in dirs_to_clean:
            # Find all raw .tck files (exclude _cleaned, _SF, afd_samples)
            tck_files = sorted([
                f.name for f in pt_path.glob("*.tck")
                if not f.name.endswith('_cleaned.tck')
                and not f.name.endswith('_SF.tck')
            ])
            
            if not tck_files:
                print(f"  [{dir_label}] No .tck files found")
                continue
            
            print(f"  [{dir_label}] Found {len(tck_files)} tracts to process")
            total_tracts += len(tck_files)
            
            results = run_batch_matlab_cleaning(
                pt_path, tck_files, params, paths_cfg,
                dry_run=dry_run, force=force
            )
            
            for fname, (success, n_orig, n_clean) in results.items():
                if success and n_orig > 0:
                    total_cleaned += 1
                elif not success:
                    if n_orig == 0:
                        total_skipped_low += 1
                    else:
                        total_errors += 1
    
    elapsed_total = time.time() - t_global
    
    print(f"\n{'='*70}")
    print(f"CLEANING SUMMARY")
    print(f"{'='*70}")
    print(f"Subjects processed: {len(subjects_to_process)}")
    print(f"Total tracts found: {total_tracts}")
    print(f"Successfully cleaned: {total_cleaned}")
    print(f"Skipped (too few fibers): {total_skipped_low}")
    print(f"Errors: {total_errors}")
    print(f"Total time: {elapsed_total:.1f}s ({elapsed_total/60:.1f} min)")


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
    parser.add_argument('--passthrough', action='store_true',
                        help='Clean passthrough tracts (batch MATLAB per subject)')
    parser.add_argument('--shrink', action='store_true',
                        help='Use passthrough_tracts_shrink/ (AP/APTR only, with --passthrough)')
    parser.add_argument('--force', action='store_true',
                        help='Re-clean even if _cleaned files already exist')
    
    args = parser.parse_args()
    
    if args.passthrough:
        run_passthrough_cleaning(
            group=args.group,
            subject=args.subject,
            dry_run=args.dry_run,
            use_shrink=args.shrink,
            force=args.force
        )
        return
    
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
