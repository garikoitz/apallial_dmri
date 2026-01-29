#!/usr/bin/env python3
"""
Tractography generation script for flip1 pipeline.

Generates streamlines from graft segments (AP/APTR) or isocortex segments (CTRL)
to target ROIs (thalamus, striatum, etc.) using tckgen.

All parameters are read from config1.yaml via config_loader.

Usage:
    # Process single subject
    python functions/generate_tracts.py --subject sub-Apallial17NmF
    
    # Process all subjects in a group
    python functions/generate_tracts.py --group APTR
    
    # Process specific target and hemisphere
    python functions/generate_tracts.py --subject sub-Apallial17NmF --target thalamus --hemisphere left
    
    # Use 9-segment mode (default is 6)
    python functions/generate_tracts.py --subject sub-Apallial17NmF --seg-mode 9
    
    # Dry run to preview commands
    python functions/generate_tracts.py --group AP --dry-run

Author: MBRAIN Project
Pipeline: flip1 (config1.yaml)
"""

import os
import subprocess
import argparse
from pathlib import Path

# Set pipeline to flip1
os.environ['MBRAIN_PIPELINE'] = 'flip1'

from config_loader import (
    load_config, get_paths, get_subject_groups,
    get_tractography_params, get_matlab_cleaning_params
)


def get_config():
    """Load and return config data."""
    config = load_config()
    paths = get_paths()
    tractography_params = get_tractography_params()
    matlab_params = get_matlab_cleaning_params()
    
    return {
        'config': config,
        'base_path': paths['derivatives_dir'],
        'algorithm': tractography_params.get('algorithm', 'iFOD2'),
        'streamlines': tractography_params.get('select', 2000),
        'cutoff': tractography_params.get('cutoff', 0.1),
        'maxlength': tractography_params.get('max_length', 15),
        'minlength': tractography_params.get('min_length', 2),
        'step_size': tractography_params.get('step_size', 0.02),
        'angle': tractography_params.get('angle', 60),
        # MATLAB cleaning parameters
        'matlab_bin': matlab_params.get('matlab_bin', '/Applications/MATLAB_R2025b.app/bin/matlab'),
        'afq_path': matlab_params.get('afq_path', '/Users/glerma/soft/gitlab/rtp2-pipeline'),
        'max_dist': matlab_params.get('max_dist', 3),
        'max_len': matlab_params.get('max_len', 3),
        'num_nodes': matlab_params.get('num_nodes', 100),
        'max_iter': matlab_params.get('max_iter', 3),
        'min_streamlines': matlab_params.get('min_streamlines', 20),
        'matlab_timeout': matlab_params.get('timeout', 300),
    }


def get_subjects_from_config():
    """Load subject definitions from config1.yaml."""
    groups = get_subject_groups()
    subjects = {}
    for group_name, subject_list in groups.items():
        subjects[group_name] = [(s['id'], s['session']) for s in subject_list]
    return subjects


# Segment names for 6-segment mode (default)
SEGMENTS_6 = [
    "anterior_dorsal",
    "anterior_ventral",
    "central_dorsal",
    "central_ventral",
    "posterior_dorsal",
    "posterior_ventral",
]

# Segment names for 9-segment mode
SEGMENTS_9 = [
    "anterior_dorsal",
    "anterior_medial",
    "anterior_ventral",
    "central_dorsal",
    "central_medial",
    "central_ventral",
    "posterior_dorsal",
    "posterior_medial",
    "posterior_ventral",
]


def run_tckgen(cfg, fod_path, seed_roi, include_roi, output_tck, mask_path=None, dry_run=False):
    """Run tckgen from seed to include ROI using config parameters.
    
    Key options:
        -stop: Forces algorithm to terminate streamline as soon as it enters an include region.
               This prevents overshooting past the target ROI.
        -seed_unidirectional: Prevents "backward" growth from seed point that often causes
                              streamlines to overshoot in the wrong direction.
    """
    cmd = [
        "tckgen", fod_path, output_tck,
        "-algorithm", cfg['algorithm'],
        "-seed_image", seed_roi,
        "-include", include_roi,
        "-select", str(cfg['streamlines']),
        "-cutoff", str(cfg['cutoff']),
        "-maxlength", str(cfg['maxlength']),
        "-minlength", str(cfg['minlength']),
        "-step", str(cfg['step_size']),
        "-angle", str(cfg['angle']),
        "-stop",                  # Terminate when entering include region (no overshooting)
        "-seed_unidirectional",   # Prevent backward growth from seed (cleaner tracts)
        "-force"
    ]
    if mask_path and os.path.exists(mask_path):
        cmd.extend(["-mask", mask_path])
    
    if dry_run:
        print(f"      [DRY RUN] tckgen {os.path.basename(seed_roi)} -> {os.path.basename(output_tck)}")
        return True
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"      ERROR: {result.stderr}")
        return False
    return True


def get_track_count(tck_path):
    """Get streamline count from tck file."""
    result = subprocess.run(["tckinfo", tck_path], capture_output=True, text=True)
    for line in result.stdout.split('\n'):
        if 'count' in line.lower() and ':' in line:
            try:
                return int(line.split(':')[1].strip())
            except:
                pass
    return 0


def run_matlab_cleaning(cfg, input_tck, output_tck, dry_run=False):
    """Run MATLAB AFQ cleaning on a tract file.
    
    Uses AFQ_removeFiberOutliers from vistasoft to remove outlier streamlines.
    All parameters are read from config1.yaml matlab_cleaning section.
    
    Args:
        cfg: Config dictionary with MATLAB parameters
        input_tck: Input tract file
        output_tck: Output cleaned tract file
        dry_run: If True, show commands without executing
    """
    # Check streamline count first
    count = get_track_count(input_tck)
    if count < cfg['min_streamlines']:
        print(f"      SKIP cleaning (only {count} fibers, need ≥{cfg['min_streamlines']})")
        return False
    
    if dry_run:
        print(f"      [DRY RUN] Would run MATLAB AFQ cleaning")
        return True
    
    dwidir = os.path.dirname(input_tck)
    tckfilename = os.path.basename(input_tck)
    output_name = os.path.basename(output_tck)
    
    # MATLAB command using AFQ_removeFiberOutliers
    matlab_cmd = (
        f"addpath(genpath('{cfg['afq_path']}')); "
        f"dwidir='{dwidir}'; "
        f"tckfilename='{tckfilename}'; "
        f"fg=fgRead(fullfile(dwidir,tckfilename)); "
        f"try, "
        f"[fg_clean,~,SF]=AFQ_removeFiberOutliers(fg,{cfg['max_dist']},{cfg['max_len']},{cfg['num_nodes']},'mean',0,{cfg['max_iter']},0); "
        f"AFQ_fgWrite(fg_clean,fullfile(dwidir,'{output_name}'),'tck'); "
        f"fprintf('Cleaned: %d -> %d fibers\\n',length(fg.fibers),length(fg_clean.fibers)); "
        f"catch ME, "
        f"fprintf('Cleaning failed: %s. Saving original.\\n', ME.message); "
        f"AFQ_fgWrite(fg,fullfile(dwidir,'{output_name}'),'tck'); "
        f"end"
    )
    
    cmd = [cfg['matlab_bin'], "-nodisplay", "-nosplash", "-batch", matlab_cmd]
    
    print(f"      Cleaning with MATLAB AFQ...", end=" ", flush=True)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=cfg['matlab_timeout'])
        
        if result.returncode != 0:
            print(f"FAILED")
            if result.stderr:
                print(f"      STDERR: {result.stderr[:300]}")
            return False
        
        if os.path.exists(output_tck):
            cleaned_count = get_track_count(output_tck)
            print(f"→ {cleaned_count} streamlines (was {count})")
            return True
        else:
            print(f"FAILED (no output)")
            return False
    except subprocess.TimeoutExpired:
        print(f"TIMEOUT after {cfg['matlab_timeout']}s")
        return False


def process_graft_subject(cfg, sub, ses, group, targets=None, hemispheres=None, 
                          seg_filter=None, seg_mode=6, clean=False, dry_run=False):
    """Process AP/APTR subject - tractography from graft segments to target ROIs.
    
    Paths for flip1 pipeline:
        FOD: dwi/wm_fod_norm.mif
        Mask: dwi/mask.mif
        Segments: dwi/graft_segments_fod/seg{6,9}_graft_{hemi}_{segment}_fod.mif
        Target ROIs: manual_roi/fod_space/roi_{Target}_{hemi}_flip_fod.mif
    """
    print(f"\n{'='*60}")
    print(f"Tractography {group}: {sub}/{ses} (seg{seg_mode} mode)")
    print(f"{'='*60}")
    
    if targets is None:
        targets = ["Thalamus", "Striatum"]
    if hemispheres is None:
        hemispheres = ["left", "right"]
    
    # Choose segment list based on mode
    segments = SEGMENTS_9 if seg_mode == 9 else SEGMENTS_6
    prefix = f"seg{seg_mode}"
    
    session_dir = os.path.join(cfg['base_path'], sub, ses)
    dwi_dir = os.path.join(session_dir, "dwi")
    
    # FOD image (flip1 pipeline uses wm_fod_norm.mif)
    fod_path = os.path.join(dwi_dir, "wm_fod_norm.mif")
    if not os.path.exists(fod_path):
        print(f"  ERROR: FOD not found: {fod_path}")
        return
    
    # Mask for tractography (flip1 uses mask.mif)
    mask_path = os.path.join(dwi_dir, "mask.mif")
    if not os.path.exists(mask_path):
        mask_path = None
        print(f"  WARNING: mask.mif not found, proceeding without mask")
    
    # Output directory
    output_dir = os.path.join(dwi_dir, "segment_tracts")
    os.makedirs(output_dir, exist_ok=True)
    
    # Graft segments directory (flip1)
    seg_dir = os.path.join(dwi_dir, "graft_segments_fod")
    if not os.path.exists(seg_dir):
        print(f"  ERROR: Graft segments not found: {seg_dir}")
        return
    
    total_tracts = 0
    total_streamlines = 0
    
    # Process each target
    for target in targets:
        # Process each hemisphere
        for seed_hemi in hemispheres:
            target_hemi = seed_hemi  # Ipsilateral
            
            # Target ROI (flip1 uses roi_{Target}_{hemi}_flip_fod.mif)
            target_roi = os.path.join(session_dir, "manual_roi", "fod_space",
                                      f"roi_{target}_{target_hemi}_flip_fod.mif")
            
            if not os.path.exists(target_roi):
                print(f"  WARNING: {target} ROI not found for {target_hemi}: {target_roi}")
                continue
            
            print(f"\n  {seed_hemi.upper()} graft → {target_hemi.upper()} {target}:")
            
            # Process segments (filter if specified)
            segs_to_process = [seg_filter] if seg_filter else segments
            
            for seg in segs_to_process:
                # Flip1 naming: seg6_graft_{hemi}_{segment}_fod.mif
                seed_roi = os.path.join(seg_dir, f"{prefix}_graft_{seed_hemi}_{seg}_fod.mif")
                
                if not os.path.exists(seed_roi):
                    print(f"    {seg}: SKIP (no segment)")
                    continue
                
                # Output filename
                output_tck = os.path.join(output_dir,
                    f"tracks_{prefix}_graft_{seed_hemi}_{seg}_to_{target}_{target_hemi}.tck")
                
                print(f"    {seg}...", end=" ", flush=True)
                
                if run_tckgen(cfg, fod_path, seed_roi, target_roi, output_tck,
                              mask_path=mask_path, dry_run=dry_run):
                    if not dry_run:
                        count = get_track_count(output_tck)
                        print(f"{count} streamlines")
                        total_tracts += 1
                        total_streamlines += count
                        
                        # Run MATLAB cleaning if requested
                        if clean and count > 0:
                            output_cleaned = output_tck.replace('.tck', '_cleaned.tck')
                            run_matlab_cleaning(cfg, output_tck, output_cleaned, dry_run=dry_run)
                    else:
                        print("OK")
                        if clean:
                            print(f"      [DRY RUN] Would run MATLAB cleaning")
                else:
                    print("FAILED")
    
    if not dry_run:
        print(f"\n  SUMMARY: {total_tracts} tract files, {total_streamlines} total streamlines")


def process_isocortex_subject(cfg, sub, ses, targets=None, hemispheres=None, 
                              seg_filter=None, clean=False, dry_run=False):
    """Process CTRL subject - tractography from isocortex segments to target ROIs.
    
    Paths for flip1 pipeline:
        FOD: dwi/wm_fod_norm.mif
        Mask: dwi/mask.mif
        Segments: dwi/isocortex_segments_fod/isocortex_{hemi}_{segment}_fod.mif
        Target ROIs: dwi/atlas/rois/fod_space/roi_{target}_{hemi}_fod.mif
    """
    print(f"\n{'='*60}")
    print(f"Tractography CTRL: {sub}/{ses}")
    print(f"{'='*60}")
    
    if targets is None:
        targets = ["thalamus", "striatum"]
    if hemispheres is None:
        hemispheres = ["left", "right"]
    
    # Use 6-segment mode for isocortex
    segments = SEGMENTS_6
    
    session_dir = os.path.join(cfg['base_path'], sub, ses)
    dwi_dir = os.path.join(session_dir, "dwi")
    
    # FOD image (flip1 pipeline)
    fod_path = os.path.join(dwi_dir, "wm_fod_norm.mif")
    if not os.path.exists(fod_path):
        print(f"  ERROR: FOD not found: {fod_path}")
        return
    
    # Mask for tractography
    mask_path = os.path.join(dwi_dir, "mask.mif")
    if not os.path.exists(mask_path):
        mask_path = None
    
    # Output directory
    output_dir = os.path.join(dwi_dir, "segment_tracts")
    os.makedirs(output_dir, exist_ok=True)
    
    # Isocortex segments directory
    seg_dir = os.path.join(dwi_dir, "isocortex_segments_fod")
    if not os.path.exists(seg_dir):
        print(f"  ERROR: Isocortex segments not found: {seg_dir}")
        return
    
    total_tracts = 0
    total_streamlines = 0
    
    # Process each target
    for target in targets:
        # Process each hemisphere
        for seed_hemi in hemispheres:
            target_hemi = seed_hemi  # Ipsilateral
            
            # Target ROI
            target_roi = os.path.join(dwi_dir, "atlas", "rois", "fod_space",
                                      f"roi_{target}_{target_hemi}_fod.mif")
            
            if not os.path.exists(target_roi):
                print(f"  WARNING: {target} ROI not found for {target_hemi}")
                continue
            
            print(f"\n  {seed_hemi.upper()} isocortex → {target_hemi.upper()} {target}:")
            
            # Process segments (filter if specified)
            segs_to_process = [seg_filter] if seg_filter else segments
            
            for seg in segs_to_process:
                seed_roi = os.path.join(seg_dir, f"isocortex_{seed_hemi}_{seg}_fod.mif")
                
                if not os.path.exists(seed_roi):
                    print(f"    {seg}: SKIP (no segment)")
                    continue
                
                output_tck = os.path.join(output_dir,
                    f"tracks_isocortex_{seed_hemi}_{seg}_to_{target}_{target_hemi}.tck")
                
                print(f"    {seg}...", end=" ", flush=True)
                
                if run_tckgen(cfg, fod_path, seed_roi, target_roi, output_tck,
                              mask_path=mask_path, dry_run=dry_run):
                    if not dry_run:
                        count = get_track_count(output_tck)
                        print(f"{count} streamlines")
                        total_tracts += 1
                        total_streamlines += count
                        
                        # Run MATLAB cleaning if requested
                        if clean and count > 0:
                            output_cleaned = output_tck.replace('.tck', '_cleaned.tck')
                            run_matlab_cleaning(cfg, output_tck, output_cleaned, dry_run=dry_run)
                    else:
                        print("OK")
                        if clean:
                            print(f"      [DRY RUN] Would run MATLAB cleaning")
                else:
                    print("FAILED")
    
    if not dry_run:
        print(f"\n  SUMMARY: {total_tracts} tract files, {total_streamlines} total streamlines")


def main():
    parser = argparse.ArgumentParser(
        description="Tractography generation for flip1 pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single subject (APTR):
  python functions/generate_tracts.py --subject sub-Apallial17NmF
  
  # Process all subjects in a group:
  python functions/generate_tracts.py --group APTR
  
  # Specific target and hemisphere:
  python functions/generate_tracts.py --subject sub-Apallial17NmF --target Thalamus --hemisphere left
  
  # Use 9-segment mode:
  python functions/generate_tracts.py --subject sub-Apallial17NmF --seg-mode 9
  
  # Dry run to preview:
  python functions/generate_tracts.py --group AP --dry-run
  
  # Run with MATLAB cleaning:
  python functions/generate_tracts.py --subject sub-Apallial17NmF --segment anterior_dorsal --clean
"""
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--group", choices=["APTR", "AP", "CTRL", "all"], default=None,
                       help="Process all subjects in group")
    parser.add_argument("--subject", help="Process specific subject (sub-xxx)")
    parser.add_argument("--target", default="all",
                       help="Target ROI (Thalamus, Striatum, or all)")
    parser.add_argument("--hemisphere", choices=["left", "right", "both"], default="both",
                       help="Hemisphere to process (default: both)")
    parser.add_argument("--segment", help="Specific segment (e.g., anterior_dorsal)")
    parser.add_argument("--seg-mode", type=int, choices=[6, 9], default=6,
                       help="Segmentation mode: 6 (3x2) or 9 (3x3) segments (default: 6)")
    parser.add_argument("--clean", action="store_true",
                       help="Run MATLAB AFQ cleaning after tractography")
    args = parser.parse_args()
    
    # Require either --subject or --group
    if not args.subject and not args.group:
        parser.error("Either --subject or --group is required")
    
    # Load config
    cfg = get_config()
    subjects = get_subjects_from_config()
    
    # Parse hemispheres
    hemispheres = ["left", "right"] if args.hemisphere == "both" else [args.hemisphere]
    
    # Parse targets (capitalize for AP/APTR manual ROI naming)
    if args.target == "all":
        targets_aptr = ["Thalamus", "Striatum"]
        targets_ctrl = ["thalamus", "striatum"]
    else:
        targets_aptr = [args.target.capitalize()]
        targets_ctrl = [args.target.lower()]
    
    print("="*60)
    print("TRACTOGRAPHY: flip1 Pipeline")
    print("="*60)
    print(f"Config: config1.yaml")
    print(f"Parameters:")
    print(f"  algorithm:  {cfg['algorithm']}")
    print(f"  select:     {cfg['streamlines']} streamlines")
    print(f"  step_size:  {cfg['step_size']}mm")
    print(f"  angle:      {cfg['angle']}°")
    print(f"  cutoff:     {cfg['cutoff']}")
    print(f"  minlength:  {cfg['minlength']}mm")
    print(f"  maxlength:  {cfg['maxlength']}mm")
    print(f"Base path:    {cfg['base_path']}")
    print(f"Hemispheres:  {hemispheres}")
    print(f"Targets:      {targets_aptr}")
    print(f"Seg mode:     {args.seg_mode} segments")
    print(f"Clean:        {args.clean}")
    
    if args.dry_run:
        print("\n*** DRY RUN MODE - No actual tractography will run ***\n")
    
    # Build list of subjects to process
    subjects_to_process = []
    
    if args.subject:
        # Find subject in any group
        for group_name in ["APTR", "AP", "CTRL"]:
            for sub, ses in subjects.get(group_name, []):
                if sub == args.subject:
                    subjects_to_process.append((sub, ses, group_name))
                    break
        if not subjects_to_process:
            print(f"ERROR: Subject {args.subject} not found in config1.yaml")
            return
    else:
        # Process all subjects in specified groups
        groups_to_process = ["APTR", "AP", "CTRL"] if args.group == "all" else [args.group]
        for group_name in groups_to_process:
            for sub, ses in subjects.get(group_name, []):
                subjects_to_process.append((sub, ses, group_name))
    
    # Process subjects
    for sub, ses, group_name in subjects_to_process:
        if group_name in ["APTR", "AP"]:
            process_graft_subject(cfg, sub, ses, group_name, targets=targets_aptr,
                                  hemispheres=hemispheres, seg_filter=args.segment,
                                  seg_mode=args.seg_mode, clean=args.clean, dry_run=args.dry_run)
        else:  # CTRL
            process_isocortex_subject(cfg, sub, ses, targets=targets_ctrl,
                                      hemispheres=hemispheres, seg_filter=args.segment,
                                      clean=args.clean, dry_run=args.dry_run)
    
    print("\n" + "="*60)
    print("TRACTOGRAPHY COMPLETE!")
    print("="*60)


if __name__ == "__main__":
    main()
