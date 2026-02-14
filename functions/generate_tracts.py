#!/usr/bin/env python3
"""
Tractography generation script for flip1 pipeline.

Generates streamlines from seed segments (Graft or Isocortex)
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
# Set pipeline to flip1 *before* any other imports that might load config
os.environ['MBRAIN_PIPELINE'] = 'flip1'

import subprocess
import argparse
import sys
from pathlib import Path
import logging

# Add parent dir to path to import utils
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import setup_logger

from functions.config_loader import (
    load_config, get_paths, get_subject_groups,
    get_tractography_params, get_matlab_cleaning_params
)

# Global logger initialization
logger = logging.getLogger(__name__)


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


# Segment definitions
SEGMENTS_6 = [
    "anterior_dorsal", "anterior_ventral",
    "central_dorsal", "central_ventral",
    "posterior_dorsal", "posterior_ventral",
]

SEGMENTS_9 = [
    "anterior_dorsal", "anterior_medial", "anterior_ventral",
    "central_dorsal", "central_medial", "central_ventral",
    "posterior_dorsal", "posterior_medial", "posterior_ventral",
]

SEGMENTS_5_AP = [
    "anterior", "anterior_central", "central", "posterior_central", "posterior",
]


class TractographyStrategy:
    """Base strategy for subject processing."""
    def __init__(self, session_dir, args):
        self.session_dir = session_dir
        self.dwi_dir = os.path.join(session_dir, "dwi")
        self.args = args

    def get_fod_path(self):
        return os.path.join(self.dwi_dir, "wm_fod_norm.mif")

    def get_mask_path(self):
        mask = os.path.join(self.dwi_dir, "mask_anat.mif")
        return mask if os.path.exists(mask) else None

    def get_output_dir(self):
        return os.path.join(self.dwi_dir, "segment_tracts")

    def get_segments(self):
        raise NotImplementedError

    def get_target_roi(self, target, hemi):
        raise NotImplementedError

    def get_seed_roi(self, hemi, segment):
        raise NotImplementedError

    def get_output_filename(self, seed_hemi, segment, target, target_hemi):
        raise NotImplementedError


class GraftStrategy(TractographyStrategy):
    """Strategy for AP/APTR subjects with manual graft segments."""
    def __init__(self, session_dir, args):
        super().__init__(session_dir, args)
        self.seg_dir = os.path.join(self.dwi_dir, "graft_segments_fod")
        self.prefix = f"seg{args.seg_mode}"

    def get_segments(self):
        if self.args.segment:
            return [self.args.segment]
        return SEGMENTS_9 if self.args.seg_mode == 9 else SEGMENTS_6

    def get_target_roi(self, target, hemi):
        # Graft targets are capitalized (Thalamus/Striatum)
        target_name = target.capitalize()
        return os.path.join(self.session_dir, "manual_roi", "fod_space",
                            f"roi_{target_name}_{hemi}_flip_fod.mif")

    def get_seed_roi(self, hemi, segment):
        return os.path.join(self.seg_dir, f"{self.prefix}_graft_{hemi}_{segment}_fod.mif")

    def get_output_filename(self, seed_hemi, segment, target, target_hemi):
        target_name = target.capitalize()
        return f"tracks_{self.prefix}_graft_{seed_hemi}_{segment}_to_{target_name}_{target_hemi}.tck"


class IsocortexStrategy(TractographyStrategy):
    """Strategy for CTRL subjects with atlas-based isocortex segments."""
    def __init__(self, session_dir, args):
        super().__init__(session_dir, args)
        self.seg_dir = os.path.join(self.dwi_dir, "isocortex_segments_v2")

    def get_segments(self):
        if self.args.segment:
            return [self.args.segment]
        return SEGMENTS_5_AP

    def get_target_roi(self, target, hemi):
        # Isocortex targets are lowercase (thalamus/striatum)
        target_name = target.lower()
        return os.path.join(self.session_dir, "atlas", "rois", "fod_space",
                            f"roi_{target_name}_{hemi}.nii.gz")

    def get_seed_roi(self, hemi, segment):
        return os.path.join(self.seg_dir, f"isocortex_{hemi}_{segment}.nii.gz")

    def get_output_filename(self, seed_hemi, segment, target, target_hemi):
        target_name = target.lower()
        return f"tracks_isocortex_{seed_hemi}_{segment}_to_{target_name}_{target_hemi}.tck"


def run_tckgen(cfg, fod_path, seed_roi, include_roi, output_tck, mask_path=None, dry_run=False):
    """Run tckgen from seed to include ROI."""
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
        "-stop",                  # Terminate when entering include region
        "-seed_unidirectional",   # Prevent backward growth
        "-force"
    ]
    if mask_path:
        cmd.extend(["-mask", mask_path])
    
    if dry_run:
        logger.info(f"      [DRY RUN] tckgen {os.path.basename(seed_roi)} -> {os.path.basename(output_tck)}")
        return True
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"      ERROR: {result.stderr}")
        return False
    return True


def run_tcksift2(fod_path, input_tck, output_weights, dry_run=False):
    """Run SIFT2 to compute streamline weights."""
    cmd = [
        "tcksift2", input_tck, fod_path, output_weights,
        "-force"
    ]
    if dry_run:
        logger.info(f"      [DRY RUN] tcksift2 {os.path.basename(input_tck)}")
        return True
    
    logger.info(f"      Running SIFT2...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"FAILED: {result.stderr[:200]}")
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
    """Run MATLAB AFQ cleaning on a tract file."""
    count = get_track_count(input_tck)
    if count < cfg['min_streamlines']:
        logger.info(f"      SKIP cleaning (only {count} fibers, need ≥{cfg['min_streamlines']})")
        return False
    
    if dry_run:
        logger.info(f"      [DRY RUN] Would run MATLAB AFQ cleaning")
        return True
    
    dwidir = os.path.dirname(input_tck)
    tckfilename = os.path.basename(input_tck)
    output_name = os.path.basename(output_tck)
    
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
    
    logger.info(f"      Cleaning with MATLAB AFQ...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=cfg['matlab_timeout'])
        if result.returncode != 0:
            logger.error(f"FAILED MATLAB: {result.stderr[:300]}")
            return False
        
        if os.path.exists(output_tck):
            cleaned_count = get_track_count(output_tck)
            logger.info(f"      → {cleaned_count} streamlines (was {count})")
            return True
        else:
            logger.error(f"FAILED (no output)")
            return False
    except subprocess.TimeoutExpired:
        logger.error(f"TIMEOUT after {cfg['matlab_timeout']}s")
        return False


def process_subject_tractography(cfg, sub, ses, group, args):
    """Unified tractography processing using Strategy pattern."""
    logger.info(f"{'='*60}")
    logger.info(f"Tractography {group}: {sub}/{ses}")
    logger.info(f"{'='*60}")

    session_dir = os.path.join(cfg['base_path'], sub, ses)
    
    # Initialize strategy based on group
    if group in ["AP", "APTR"]:
        strategy = GraftStrategy(session_dir, args)
    else:
        strategy = IsocortexStrategy(session_dir, args)

    fod_path = strategy.get_fod_path()
    if not os.path.exists(fod_path):
        logger.error(f"  ERROR: FOD not found: {fod_path}")
        return

    mask_path = strategy.get_mask_path()
    output_dir = strategy.get_output_dir()
    os.makedirs(output_dir, exist_ok=True)
    
    segments = strategy.get_segments()
    targets = [args.target] if args.target != "all" else ["thalamus", "striatum"]
    hemispheres = [args.hemisphere] if args.hemisphere != "both" else ["left", "right"]

    total_tracts = 0
    total_streamlines = 0

    for target in targets:
        for seed_hemi in hemispheres:
            target_hemi = seed_hemi  # Ipsilateral
            
            target_roi = strategy.get_target_roi(target, seed_hemi)
            if not os.path.exists(target_roi):
                logger.warning(f"  WARNING: Target ROI not found: {target_roi}")
                continue
            
            logger.info(f"  {seed_hemi.upper()} seed → {target_hemi.upper()} {target}:")
            
            for seg in segments:
                seed_roi = strategy.get_seed_roi(seed_hemi, seg)
                if not os.path.exists(seed_roi):
                    logger.info(f"    {seg}: SKIP (no seed segment)")
                    continue
                
                output_filename = strategy.get_output_filename(seed_hemi, seg, target, target_hemi)
                output_tck = os.path.join(output_dir, output_filename)
                
                logger.info(f"    {seg}...")
                
                if run_tckgen(cfg, fod_path, seed_roi, target_roi, output_tck, mask_path, args.dry_run):
                    if not args.dry_run:
                        count = get_track_count(output_tck)
                        logger.info(f"      -> {count} streamlines")
                        total_tracts += 1
                        total_streamlines += count
                        
                        # SIFT2
                        if args.sift2 and count > 0:
                            weights_txt = output_tck.replace('.tck', '_sift2_weights.txt')
                            run_tcksift2(fod_path, output_tck, weights_txt, dry_run=args.dry_run)
                            
                        # Cleaning
                        if args.clean and count > 0:
                            output_cleaned = output_tck.replace('.tck', '_cleaned.tck')
                            run_matlab_cleaning(cfg, output_tck, output_cleaned, dry_run=args.dry_run)
                    else:
                        if args.sift2: logger.info(f"      [DRY RUN] Would run SIFT2")
                        if args.clean: logger.info(f"      [DRY RUN] Would run MATLAB cleaning")
                else:
                    logger.error("FAILED")

    if not args.dry_run:
        logger.info(f"  SUMMARY: {total_tracts} tract files, {total_streamlines} total streamlines")


def main():
    parser = argparse.ArgumentParser(description="Tractography generation for flip1 pipeline")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--group", choices=["APTR", "AP", "CTRL", "all"], default=None,
                       help="Process all subjects in group")
    parser.add_argument("--subject", help="Process specific subject (sub-xxx)")
    parser.add_argument("--target", default="all",
                       help="Target ROI (thalamus, striatum, or all)")
    parser.add_argument("--hemisphere", choices=["left", "right", "both"], default="both",
                       help="Hemisphere to process (default: both)")
    parser.add_argument("--segment", help="Specific segment (e.g., anterior_dorsal)")
    parser.add_argument("--seg-mode", type=int, choices=[6, 9], default=6,
                       help="Segmentation mode: 6 (3x2) or 9 (3x3) segments (default: 6)")
    parser.add_argument("--clean", action="store_true",
                       help="Run MATLAB AFQ cleaning after tractography")
    parser.add_argument("--sift2", action="store_true",
                       help="Apply SIFT2 filtering to tractography output")
    parser.add_argument("--log-file", help="Path to log file")
    
    args = parser.parse_args()
    
    # Setup Logger
    global logger
    logger = setup_logger(__name__, args.log_file)
    
    if not args.subject and not args.group:
        parser.error("Either --subject or --group is required")
    
    cfg = get_config()
    subjects = get_subjects_from_config()
    
    logger.info("="*60)
    logger.info("TRACTOGRAPHY: flip1 Pipeline")
    logger.info("="*60)
    
    if args.dry_run:
        logger.info("*** DRY RUN MODE ***")
    
    # Build list of subjects
    subjects_to_process = []
    if args.subject:
        found = False
        for group_name in ["APTR", "AP", "CTRL"]:
            if group_name in subjects:
                for sub, ses in subjects.get(group_name, []):
                    if sub == args.subject:
                        subjects_to_process.append((sub, ses, group_name))
                        found = True
                        break
            if found: break
        if not found:
            logger.error(f"Subject {args.subject} not found in config")
            return
    else:
        groups = ["APTR", "AP", "CTRL"] if args.group == "all" else [args.group]
        for group_name in groups:
            for sub, ses in subjects.get(group_name, []):
                subjects_to_process.append((sub, ses, group_name))

    for sub, ses, group_name in subjects_to_process:
        process_subject_tractography(cfg, sub, ses, group_name, args)
    
    logger.info("="*60)
    logger.info("COMPLETED")

if __name__ == "__main__":
    main()
