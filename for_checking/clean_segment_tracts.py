#!/usr/bin/env python3
"""
Clean segment-based tractography files using MATLAB AFQ_removeFiberOutliers.

Cleans fiber tracts by removing outliers based on distance and length.
Parameters are read from config.yaml.

Usage:
    python clean_segment_tracts.py                            # All segments, all subjects
    python clean_segment_tracts.py --group CTRL               # CTRL group only
    python clean_segment_tracts.py --subject sub-129-921-R1F  # Specific subject
    python clean_segment_tracts.py --dry-run                  # Preview commands
"""

import argparse
import subprocess
import os
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups,
    get_matlab_cleaning_params
)

# Segment names
SEGMENTS = [
    "anterior_dorsal", "anterior_ventral",
    "central_dorsal", "central_ventral", 
    "posterior_dorsal", "posterior_ventral"
]


def get_subjects_from_config():
    """Get subjects organized by group from config.yaml."""
    groups = get_subject_groups()
    subjects = {"APTR": [], "AP": [], "CTRL": []}
    for group, subj_list in groups.items():
        if group in subjects:
            for s in subj_list:
                subjects[group].append((s['id'], s['session']))
    return subjects


def get_streamline_count(tck_file):
    """Get streamline count from tck file."""
    result = subprocess.run(["tckinfo", str(tck_file)], capture_output=True, text=True)
    for line in result.stdout.split('\n'):
        if 'count' in line.lower() and ':' in line:
            try:
                return int(line.split(':')[1].strip())
            except:
                pass
    return 0


def run_matlab_cleaning(tracks_file, output_file, params, paths, dry_run=False):
    """Run MATLAB AFQ cleaning on a tract file."""
    # Check streamline count first
    min_streamlines = params.get('min_streamlines', 10)
    count = get_streamline_count(tracks_file)
    
    if count == 0:
        # Copy empty file (or skip)
        return False, "SKIP (0 streamlines)"
    
    if count < min_streamlines:
        # Too few streamlines - just copy the file as-is
        import shutil
        if not dry_run:
            shutil.copy(str(tracks_file), str(output_file))
        return True, f"COPY ({count} < {min_streamlines} min)"
    
    matlab_bin = paths['matlab_bin']
    afq_path = paths['afq_path']
    
    max_dist = params['max_dist']
    max_len = params['max_len']
    num_nodes = params['num_nodes']
    max_iter = params.get('max_iter', 10)
    
    dwidir = str(tracks_file.parent)
    tckfilename = tracks_file.name
    output_name = output_file.name
    
    matlab_cmd = f"""
addpath(genpath('{afq_path}'));
dwidir='{dwidir}';
tckfilename='{tckfilename}';
fg=fgRead(fullfile(dwidir,tckfilename));
[fg_clean,~,SF]=AFQ_removeFiberOutliers(fg,{max_dist},{max_len},{num_nodes},'mean',0,{max_iter},0);
AFQ_fgWrite(fg_clean,fullfile(dwidir,'{output_name}'),'tck');
fprintf('Cleaned: %d -> %d fibers\\n',length(fg.fibers),length(fg_clean.fibers));
exit;
"""
    
    if dry_run:
        print(f"      [DRY RUN] Would clean {tracks_file.name} -> {output_file.name}")
        return True, "DRY RUN"
    
    cmd = [matlab_bin, "-nodisplay", "-nosplash", "-r", matlab_cmd.replace('\n', ' ')]
    
    # Timeout after 120 seconds (2 minutes) per tract
    timeout_seconds = params.get('timeout', 120)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=dwidir, timeout=timeout_seconds)
        
        # Parse output for cleaned count
        for line in result.stdout.split('\n'):
            if 'Cleaned:' in line:
                return True, line.strip()
        
        if result.returncode != 0:
            return False, result.stderr[:100]
        
        return True, "Done"
    
    except subprocess.TimeoutExpired:
        return False, f"TIMEOUT ({timeout_seconds}s)"


def process_subject(sub, ses, group, paths, params, dry_run=False):
    """Process all segment tracts for a subject."""
    base_path = paths['derivatives_dir']
    session_dir = Path(base_path) / sub / ses / "dwi"
    tract_dir = session_dir / "segment_tracts"
    
    if not tract_dir.exists():
        print(f"    SKIP: No segment_tracts folder")
        return 0, 0
    
    # Determine seed type based on group
    if group in ["AP", "APTR"]:
        seed_type = "graft"
    else:
        seed_type = "isocortex"
    
    cleaned = 0
    failed = 0
    
    # Process both thalamus and striatum targets
    targets = ["thalamus", "striatum"]
    
    for target in targets:
        for hemi in ["left", "right"]:
            for seg in SEGMENTS:
                # Track file naming: tracks_graft_left_anterior_dorsal_to_thalamus_left.tck
                track_name = f"tracks_{seed_type}_{hemi}_{seg}_to_{target}_{hemi}.tck"
                track_file = tract_dir / track_name
                
                if not track_file.exists():
                    continue
                
                output_name = track_name.replace(".tck", "_cleaned.tck")
                output_file = tract_dir / output_name
                
                print(f"    {hemi} {seg} → {target}...", end=" ", flush=True)
                
                success, msg = run_matlab_cleaning(track_file, output_file, params, paths, dry_run)
                
                if success:
                    print(msg)
                    cleaned += 1
                else:
                    print(f"FAILED: {msg}")
                    failed += 1
    
    return cleaned, failed


def main():
    parser = argparse.ArgumentParser(description="Clean segment-based tractography")
    parser.add_argument("--group", choices=["APTR", "AP", "CTRL", "all"], default="all")
    parser.add_argument("--subject", help="Process specific subject (sub-xxx)")
    parser.add_argument("--segment", help="Specific segment (e.g., anterior_dorsal)")
    parser.add_argument("--dry-run", action="store_true", help="Preview commands")
    args = parser.parse_args()
    
    paths = get_paths()
    params = get_matlab_cleaning_params()
    subjects = get_subjects_from_config()
    
    # Filter segments if specified
    global SEGMENTS
    if args.segment:
        SEGMENTS = [args.segment]
    
    print("="*60)
    print("SEGMENT TRACT CLEANING")
    print("="*60)
    print(f"Cleaning parameters from config:")
    print(f"  max_dist:  {params['max_dist']}")
    print(f"  max_len:   {params['max_len']}")
    print(f"  num_nodes: {params['num_nodes']}")
    print(f"Base path:   {paths['derivatives_dir']}")
    print(f"Segments:    {SEGMENTS}")
    
    if args.dry_run:
        print("\n*** DRY RUN MODE ***\n")
    
    total_cleaned = 0
    total_failed = 0
    
    groups_to_process = ["APTR", "AP", "CTRL"] if args.group == "all" else [args.group]
    
    for group in groups_to_process:
        if group not in subjects:
            continue
            
        print(f"\n{'#'*60}")
        print(f"# {group} GROUP")
        print(f"{'#'*60}")
        
        for sub, ses in subjects[group]:
            if args.subject and sub != args.subject:
                continue
            
            print(f"\n  Subject: {sub}/{ses}")
            cleaned, failed = process_subject(sub, ses, group, paths, params, args.dry_run)
            total_cleaned += cleaned
            total_failed += failed
    
    print(f"\n{'='*60}")
    print(f"CLEANING COMPLETE: {total_cleaned} cleaned, {total_failed} failed")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
