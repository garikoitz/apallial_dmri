#!/usr/bin/env python3
"""
Master pipeline script for passthrough tract analysis.

This script runs the complete pipeline:
1. Extract passthrough streamlines from whole-brain tractograms
2. Run MATLAB AFQ cleaning on extracted tracts
3. Apply SIFT2 weighting
4. Calculate connectivity summaries

Usage:
    python functions/run_passthrough_pipeline.py                    # All steps, all subjects
    python functions/run_passthrough_pipeline.py --group CTRL       # CTRL only
    python functions/run_passthrough_pipeline.py --step extract     # Extraction only
    python functions/run_passthrough_pipeline.py --step clean       # Cleaning only
    python functions/run_passthrough_pipeline.py --step sift2       # SIFT2 only
    python functions/run_passthrough_pipeline.py --dry-run          # Preview
"""

import argparse
import subprocess
import os
import sys
from pathlib import Path
from datetime import datetime
import json

# Set pipeline before importing config_loader
os.environ['MBRAIN_PIPELINE'] = 'flip1'

from config_loader import (
    load_config, get_paths, get_subject_groups, 
    get_matlab_cleaning_params
)


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session, group) tuples based on filters."""
    groups = get_subject_groups()
    subjects = []
    
    if subject:
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session'], g))
                    return subjects
        raise ValueError(f"Subject {subject} not found")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found")
        for s in groups[group]:
            subjects.append((s['id'], s['session'], group))
    else:
        for g, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session'], g))
    
    return subjects


def run_extraction(subjects, deriv_dir, hemi="left", dry_run=False, use_shrink=False):
    """Run passthrough extraction for all subjects."""
    shrink_label = " (SHRUNK SEGMENTS)" if use_shrink else ""
    print("\n" + "=" * 70)
    print(f"STEP 1: PASSTHROUGH EXTRACTION ({hemi.upper()} hemisphere){shrink_label}")
    print("=" * 70)
    
    # Import and run extraction
    sys.path.insert(0, str(Path(__file__).parent))
    from extract_passthrough_tracts import process_ctrl_subject, process_apaptr_subject
    
    all_results = {}
    for subject_id, session, group in subjects:
        print(f"\n--- {subject_id} ({group}) ---")
        
        if group == "CTRL":
            results = process_ctrl_subject(subject_id, session, deriv_dir, hemi, dry_run)
        else:
            results = process_apaptr_subject(subject_id, session, group, deriv_dir, hemi, dry_run, use_shrink)
        
        all_results[subject_id] = results
    
    return all_results


def run_matlab_cleaning_single(track_file, output_file, dry_run=False):
    """Run MATLAB AFQ cleaning on a single tract file."""
    paths = get_paths()
    params = get_matlab_cleaning_params()
    
    matlab_bin = paths['matlab_path']
    afq_path = paths['rtp2_pipeline_path']
    
    max_dist = params['max_dist']
    max_len = params['max_len']
    num_nodes = params['num_nodes']
    max_iter = params.get('max_iter', 3)
    min_streamlines = params.get('min_streamlines', 20)
    timeout = params.get('timeout', 300)
    
    # Check streamline count
    info_result = subprocess.run(
        ["tckinfo", str(track_file), "-count"],
        capture_output=True, text=True
    )
    count = 0
    for line in info_result.stdout.split('\n'):
        if 'count:' in line and 'total_count' not in line:
            count = int(line.split()[-1])
            break
    
    if count < min_streamlines:
        print(f"        SKIP: Only {count} streamlines (min: {min_streamlines})")
        return -1
    
    dwidir = str(track_file.parent)
    tckfilename = track_file.name
    output_name = output_file.name
    sf_name = output_name.replace('_cleaned.tck', '_SF.tck')
    
    matlab_cmd = (
        f"addpath(genpath('{afq_path}')); "
        f"dwidir='{dwidir}'; "
        f"tckfilename='{tckfilename}'; "
        f"fg=fgRead(fullfile(dwidir,tckfilename)); "
        f"[fg_clean,~,SF]=AFQ_removeFiberOutliers(fg,{max_dist},{max_len},{num_nodes},'mean',0,{max_iter},0); "
        f"fg_SF=dtiNewFiberGroup([fg.name '_SF'],[0 0 155],[],[],SF.fibers); "
        f"AFQ_fgWrite(fg_clean,fullfile(dwidir,'{output_name}'),'tck'); "
        f"AFQ_fgWrite(fg_SF,fullfile(dwidir,'{sf_name}'),'tck'); "
        f"fprintf('Cleaned: %d -> %d fibers\\n',length(fg.fibers),length(fg_clean.fibers));"
    )
    
    cmd = [matlab_bin, "-nodisplay", "-nosplash", "-batch", matlab_cmd]
    
    if dry_run:
        print(f"        [DRY RUN] Would run MATLAB cleaning")
        return 0
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        
        if result.returncode != 0:
            print(f"        ERROR: MATLAB failed")
            return -1
        
        if output_file.exists():
            # Get cleaned count
            info_result = subprocess.run(
                ["tckinfo", str(output_file), "-count"],
                capture_output=True, text=True
            )
            for line in info_result.stdout.split('\n'):
                if 'count:' in line and 'total_count' not in line:
                    cleaned_count = int(line.split()[-1])
                    return cleaned_count
        return 0
        
    except subprocess.TimeoutExpired:
        print(f"        ERROR: MATLAB timed out (>{timeout}s)")
        return -1
    except Exception as e:
        print(f"        ERROR: {e}")
        return -1


def run_cleaning(subjects, deriv_dir, hemi="left", dry_run=False, use_shrink=False):
    """Run MATLAB cleaning on all extracted passthrough tracts."""
    shrink_label = " (SHRUNK SEGMENTS)" if use_shrink else ""
    print("\n" + "=" * 70)
    print(f"STEP 2: MATLAB AFQ CLEANING ({hemi.upper()} hemisphere){shrink_label}")
    print("=" * 70)
    
    params = get_matlab_cleaning_params()
    print(f"Parameters: maxDist={params['max_dist']}, maxLen={params['max_len']}, "
          f"numNodes={params['num_nodes']}, minStreamlines={params.get('min_streamlines', 20)}")
    
    all_results = {}
    
    # Determine which folder to use
    tract_folder = "passthrough_tracts_shrink" if use_shrink else "passthrough_tracts"
    
    for subject_id, session, group in subjects:
        print(f"\n--- {subject_id} ({group}) ---")
        
        base_dir = deriv_dir / subject_id / session
        passthrough_dir = base_dir / "dwi" / "wholebrain" / tract_folder
        
        if not passthrough_dir.exists():
            print(f"    SKIP: No {tract_folder} folder")
            continue
        
        # Find all .tck files for this hemisphere (excluding _cleaned and _SF)
        tck_files = list(passthrough_dir.glob(f"*_{hemi}_*.tck"))
        tck_files = [f for f in tck_files if '_cleaned' not in f.name and '_SF' not in f.name]
        
        results = {}
        for i, tck_file in enumerate(tck_files, 1):
            output_file = tck_file.with_name(tck_file.stem + '_cleaned.tck')
            
            if output_file.exists() and not dry_run:
                # Check if cleaned file is newer than original
                if output_file.stat().st_mtime > tck_file.stat().st_mtime:
                    print(f"    [{i}/{len(tck_files)}] {tck_file.stem}: SKIP (already cleaned)")
                    results[tck_file.stem] = "skipped"
                    continue
            
            print(f"    [{i}/{len(tck_files)}] {tck_file.stem}...", end=" ", flush=True)
            count = run_matlab_cleaning_single(tck_file, output_file, dry_run)
            
            if count >= 0:
                print(f"-> {count} streamlines")
                results[tck_file.stem] = count
            else:
                results[tck_file.stem] = "failed"
        
        all_results[subject_id] = results
    
    return all_results


def run_sift2(subjects, deriv_dir, dry_run=False):
    """Apply SIFT2 to whole-brain tractograms."""
    print("\n" + "=" * 70)
    print("STEP 3: SIFT2 WEIGHTING")
    print("=" * 70)
    
    all_results = {}
    
    for subject_id, session, group in subjects:
        print(f"\n--- {subject_id} ({group}) ---")
        
        base_dir = deriv_dir / subject_id / session
        dwi_dir = base_dir / "dwi"
        wholebrain_dir = dwi_dir / "wholebrain"
        
        wholebrain_tck = wholebrain_dir / "wholebrain_1000k.tck"
        fod_file = dwi_dir / "wm_fod_norm_group.mif"
        weights_file = wholebrain_dir / "sift2_weights.txt"
        mu_file = wholebrain_dir / "sift2_mu.txt"
        
        if not wholebrain_tck.exists():
            print(f"    SKIP: No wholebrain tractogram")
            all_results[subject_id] = "no_tractogram"
            continue
        
        if not fod_file.exists():
            print(f"    SKIP: No FOD file")
            all_results[subject_id] = "no_fod"
            continue
        
        if weights_file.exists() and not dry_run:
            print(f"    SKIP: SIFT2 weights already exist")
            all_results[subject_id] = "exists"
            continue
        
        cmd = [
            "tcksift2",
            str(wholebrain_tck),
            str(fod_file),
            str(weights_file),
            "-out_mu", str(mu_file),
            "-force"
        ]
        
        if dry_run:
            print(f"    [DRY RUN] {' '.join(cmd)}")
            all_results[subject_id] = "dry_run"
            continue
        
        print(f"    Running tcksift2...", end=" ", flush=True)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode == 0 and weights_file.exists():
                # Count weights
                with open(weights_file, 'r') as f:
                    num_weights = sum(1 for _ in f)
                print(f"Done ({num_weights} weights)")
                all_results[subject_id] = num_weights
            else:
                print(f"FAILED")
                all_results[subject_id] = "failed"
                
        except subprocess.TimeoutExpired:
            print(f"TIMEOUT")
            all_results[subject_id] = "timeout"
        except Exception as e:
            print(f"ERROR: {e}")
            all_results[subject_id] = "error"
    
    return all_results


def generate_summary(subjects, deriv_dir):
    """Generate connectivity summary for all subjects."""
    print("\n" + "=" * 70)
    print("STEP 4: CONNECTIVITY SUMMARY")
    print("=" * 70)
    
    summary = {}
    
    for subject_id, session, group in subjects:
        base_dir = deriv_dir / subject_id / session
        passthrough_dir = base_dir / "dwi" / "wholebrain" / "passthrough_tracts"
        
        if not passthrough_dir.exists():
            continue
        
        subject_summary = {
            'group': group,
            'session': session,
            'tracts': {}
        }
        
        # Get all cleaned tracts
        cleaned_files = list(passthrough_dir.glob("*_cleaned.tck"))
        
        for cleaned_file in cleaned_files:
            tract_name = cleaned_file.stem.replace('_cleaned', '')
            
            # Get streamline count
            info_result = subprocess.run(
                ["tckinfo", str(cleaned_file), "-count"],
                capture_output=True, text=True
            )
            count = 0
            for line in info_result.stdout.split('\n'):
                if 'count:' in line and 'total_count' not in line:
                    count = int(line.split()[-1])
                    break
            
            subject_summary['tracts'][tract_name] = {
                'cleaned_streamlines': count
            }
        
        summary[subject_id] = subject_summary
    
    # Save summary
    summary_file = deriv_dir / "passthrough_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nSummary saved to: {summary_file}")
    
    # Print summary
    print("\n--- CONNECTIVITY SUMMARY ---")
    for group in ['CTRL', 'AP', 'APTR']:
        group_subjects = [s for s in summary.keys() if summary[s]['group'] == group]
        if group_subjects:
            print(f"\n{group}:")
            for subj in group_subjects:
                total = sum(t.get('cleaned_streamlines', 0) for t in summary[subj]['tracts'].values())
                n_tracts = len(summary[subj]['tracts'])
                print(f"  {subj}: {n_tracts} tracts, {total} total streamlines")
    
    return summary


def main():
    parser = argparse.ArgumentParser(
        description='Run complete passthrough tract analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--group', type=str, choices=['CTRL', 'AP', 'APTR'],
                        help='Process specific group only')
    parser.add_argument('--subject', type=str, help='Process specific subject only')
    parser.add_argument('--hemi', type=str, choices=['left', 'right'], default='left',
                        help='Hemisphere to process (default: left)')
    parser.add_argument('--step', type=str, choices=['extract', 'clean', 'sift2', 'summary', 'all'],
                        default='all', help='Run specific step only')
    parser.add_argument('--dry-run', action='store_true', help='Preview commands without executing')
    parser.add_argument('--shrink', action='store_true',
                        help='Use shrunk graft segments (eroded by 1 voxel) - AP/APTR only')
    parser.add_argument('--derivatives-dir', type=str,
                        help='Override derivatives directory (default: from config1.yaml)')
    
    args = parser.parse_args()
    
    # Override derivatives directory if specified
    if args.derivatives_dir:
        os.environ['MBRAIN_DERIVATIVES_DIR'] = args.derivatives_dir
    
    # Shrink mode only applies to AP/APTR
    if args.shrink and args.group == 'CTRL':
        print("ERROR: --shrink option only applies to AP/APTR groups (CTRL has no graft)")
        sys.exit(1)
    
    # Load configuration
    cfg = load_config()
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_flip1'])
    
    shrink_label = " (SHRUNK SEGMENTS)" if args.shrink else ""
    print("=" * 70)
    print(f"PASSTHROUGH TRACT ANALYSIS PIPELINE{shrink_label}")
    print(f"Hemisphere: {args.hemi.upper()}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Pipeline: flip1")
    print(f"Derivatives: {deriv_dir}")
    print("=" * 70)
    
    # Get subjects
    subjects = get_subjects(group=args.group, subject=args.subject)
    
    # If shrink mode with no group specified, filter to AP/APTR only
    if args.shrink and not args.group:
        subjects = [(s, ses, g) for (s, ses, g) in subjects if g in ['AP', 'APTR']]
    
    print(f"\nSubjects to process: {len(subjects)} ({args.hemi} hemisphere)")
    for s in subjects:
        print(f"  - {s[0]} ({s[2]})")
    
    # Run pipeline steps
    if args.step in ['extract', 'all']:
        run_extraction(subjects, deriv_dir, args.hemi, args.dry_run, args.shrink)
    
    if args.step in ['clean', 'all']:
        run_cleaning(subjects, deriv_dir, args.hemi, args.dry_run, args.shrink)
    
    if args.step in ['sift2', 'all']:
        run_sift2(subjects, deriv_dir, args.dry_run)
    
    if args.step in ['summary', 'all'] and not args.dry_run:
        generate_summary(subjects, deriv_dir)
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)


if __name__ == "__main__":
    main()
