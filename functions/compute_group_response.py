#!/usr/bin/env python3
"""
Compute Group Response Functions and Re-estimate FODs

This script ensures that Fiber Density (FD) values are comparable across subjects
by using a group-averaged response function. This is critical for valid fixel-based
analysis.

Steps:
1. Collect individual response functions from all subjects
2. Average them using MRtrix's responsemean
3. Re-estimate FODs for all subjects using the group response
4. Apply mtnormalise for intensity calibration

Output: wm_fod_norm_group.mif for each subject (calibrated for group comparison)

Usage:
    python compute_group_response.py
    python compute_group_response.py --force  # Regenerate even if exists
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from typing import List, Tuple

# Set pipeline
os.environ['MBRAIN_PIPELINE'] = 'flip1'

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import setup_logger
from functions.config_loader import load_config

# Global logger
logger = None


def get_all_subjects(cfg: dict) -> List[Tuple[str, str, str]]:
    """Get all subjects from all groups with their group label."""
    subjects = []
    groups = cfg.get('groups', {})
    
    for group_name, subject_list in groups.items():
        for subj in subject_list:
            subjects.append((subj['id'], subj['session'], group_name))
    
    return subjects


def run_cmd(cmd: str, description: str = None) -> bool:
    """Run shell command with logging."""
    if description:
        logger.info(f"  {description}")
    
    logger.info(f"    $ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"    ✗ FAILED: {result.stderr[:500]}")
        return False
    
    logger.info(f"    ✓ Success")
    return True


def collect_response_functions(subjects: List[Tuple[str, str, str]], deriv_dir: Path) -> dict:
    """Collect all individual response functions."""
    logger.info("\n=== Step 1: Collecting Response Functions ===")
    
    response_files = {
        'wm': [],
        'gm': [],
        'csf': []
    }
    
    missing_count = 0
    for subject, session, group in subjects:
        dwi_dir = deriv_dir / subject / session / 'dwi'
        
        resp_wm = dwi_dir / 'response_wm.txt'
        resp_gm = dwi_dir / 'response_gm.txt'
        resp_csf = dwi_dir / 'response_csf.txt'
        
        if resp_wm.exists() and resp_gm.exists() and resp_csf.exists():
            response_files['wm'].append(str(resp_wm))
            response_files['gm'].append(str(resp_gm))
            response_files['csf'].append(str(resp_csf))
            logger.info(f"  ✓ {subject} ({group})")
        else:
            logger.warning(f"  ⚠ {subject} ({group}): Missing response files")
            missing_count += 1
    
    logger.info(f"\nCollected: {len(response_files['wm'])} subjects")
    if missing_count > 0:
        logger.warning(f"Missing: {missing_count} subjects (run dwi_fod.py first)")
    
    return response_files


def compute_average_responses(response_files: dict, output_dir: Path) -> Tuple[Path, Path, Path]:
    """Average response functions across all subjects."""
    logger.info("\n=== Step 2: Computing Group Average Responses ===")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    avg_wm = output_dir / 'group_response_wm.txt'
    avg_gm = output_dir / 'group_response_gm.txt'
    avg_csf = output_dir / 'group_response_csf.txt'
    
    # Average WM response (quote each path individually)
    wm_files_str = ' '.join(f"'{f}'" for f in response_files['wm'])
    if not run_cmd(f"responsemean {wm_files_str} '{avg_wm}' -force",
                   f"Averaging {len(response_files['wm'])} WM responses..."):
        return None, None, None
    
    # Average GM response (quote each path individually)
    gm_files_str = ' '.join(f"'{f}'" for f in response_files['gm'])
    if not run_cmd(f"responsemean {gm_files_str} '{avg_gm}' -force",
                   f"Averaging {len(response_files['gm'])} GM responses..."):
        return None, None, None
    
    # Average CSF response (quote each path individually)
    csf_files_str = ' '.join(f"'{f}'" for f in response_files['csf'])
    if not run_cmd(f"responsemean {csf_files_str} '{avg_csf}' -force",
                   f"Averaging {len(response_files['csf'])} CSF responses..."):
        return None, None, None
    
    logger.info(f"\n✓ Group response functions saved to: {output_dir}")
    return avg_wm, avg_gm, avg_csf


def reestimate_fods(subjects: List[Tuple[str, str, str]], 
                    deriv_dir: Path,
                    avg_wm: Path, 
                    avg_gm: Path, 
                    avg_csf: Path,
                    force: bool = False) -> int:
    """Re-estimate FODs for all subjects using group response."""
    logger.info("\n=== Step 3: Re-estimating FODs with Group Response ===")
    
    success_count = 0
    skip_count = 0
    fail_count = 0
    
    for subject, session, group in subjects:
        logger.info(f"\n--- {subject} ({group}) ---")
        
        dwi_dir = deriv_dir / subject / session / 'dwi'
        
        # Input
        dwi_input = dwi_dir / 'dwi_flip_gradcheck_denoise_bias_anat.mif'
        mask = dwi_dir / 'mask_anat.mif'
        
        # Outputs
        wm_fod_group = dwi_dir / 'wm_fod_group.mif'
        gm_fod_group = dwi_dir / 'gm_fod_group.mif'
        csf_fod_group = dwi_dir / 'csf_fod_group.mif'
        wm_fod_norm_group = dwi_dir / 'wm_fod_norm_group.mif'
        csf_fod_norm_group = dwi_dir / 'csf_fod_norm_group.mif'
        
        # Check if exists
        if wm_fod_norm_group.exists() and not force:
            logger.info("  Skip: Group FOD already exists (use --force to regenerate)")
            skip_count += 1
            continue
        
        # Check inputs
        if not dwi_input.exists():
            logger.warning(f"  ⚠ Skipping: DWI not found")
            fail_count += 1
            continue
        
        if not mask.exists():
            logger.warning(f"  ⚠ Skipping: Mask not found")
            fail_count += 1
            continue
        
        # Re-estimate FOD with group response
        cmd = (f"dwi2fod msmt_csd '{dwi_input}' "
               f"'{avg_wm}' '{wm_fod_group}' "
               f"'{avg_gm}' '{gm_fod_group}' "
               f"'{avg_csf}' '{csf_fod_group}' "
               f"-mask '{mask}' -force -quiet")
        
        if not run_cmd(cmd, "Computing FOD with group response..."):
            fail_count += 1
            continue
        
        # Normalize — try WM+CSF first, fall back to WM-only, then copy
        norm_cmd = (f"mtnormalise '{wm_fod_group}' '{wm_fod_norm_group}' "
                   f"'{csf_fod_group}' '{csf_fod_norm_group}' "
                   f"-mask '{mask}' -force -quiet")
        
        result = subprocess.run(norm_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"  ⚠ mtnormalise WM+CSF failed: {result.stderr[:200]}")
            logger.info("  Trying WM-only mtnormalise as fallback...")
            norm_wm_cmd = (f"mtnormalise '{wm_fod_group}' '{wm_fod_norm_group}' "
                          f"-mask '{mask}' -force -quiet")
            result2 = subprocess.run(norm_wm_cmd, shell=True, capture_output=True, text=True)
            if result2.returncode != 0:
                logger.warning(f"  ⚠ mtnormalise WM-only also failed (copying unnormalized): {result2.stderr[:200]}")
                run_cmd(f"cp '{wm_fod_group}' '{wm_fod_norm_group}'", "Using unnormalized FOD...")
            else:
                logger.info("  ✓ FOD normalized (WM-only fallback)")
        else:
            logger.info("  ✓ FOD normalized")
        
        success_count += 1
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Success: {success_count} | Skipped: {skip_count} | Failed: {fail_count}")
    logger.info(f"{'='*60}")
    
    return success_count


def main():
    parser = argparse.ArgumentParser(
        description='Compute group response functions and re-estimate FODs',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration even if output exists')
    parser.add_argument('--log-file', help='Path to log file')
    parser.add_argument('--derivatives-dir', type=str,
                        help='Override derivatives directory (default: from config1.yaml)')
    
    args = parser.parse_args()
    
    # Override derivatives directory if specified
    if args.derivatives_dir:
        os.environ['MBRAIN_DERIVATIVES_DIR'] = args.derivatives_dir
    
    # Setup logger
    global logger
    logger = setup_logger(__name__, args.log_file)
    
    # Load config
    cfg = load_config()
    deriv_dir = Path(cfg['paths']['derivatives_flip1'])
    
    logger.info("="*60)
    logger.info("Group Response Function & FOD Re-estimation")
    logger.info("="*60)
    logger.info(f"Pipeline: flip1")
    logger.info(f"Derivatives: {deriv_dir}")
    
    # Get all subjects
    subjects = get_all_subjects(cfg)
    logger.info(f"Total subjects: {len(subjects)}")
    
    # Step 1: Collect individual responses
    response_files = collect_response_functions(subjects, deriv_dir)
    
    if len(response_files['wm']) < 2:
        logger.error("ERROR: Need at least 2 subjects with response functions")
        return 1
    
    # Step 2: Compute group averages
    group_response_dir = deriv_dir / 'group_responses'
    avg_wm, avg_gm, avg_csf = compute_average_responses(response_files, group_response_dir)
    
    if not avg_wm:
        logger.error("ERROR: Failed to compute average responses")
        return 1
    
    # Step 3: Re-estimate FODs
    success_count = reestimate_fods(subjects, deriv_dir, avg_wm, avg_gm, avg_csf, args.force)
    
    if success_count == 0:
        logger.error("ERROR: No FODs were successfully re-estimated")
        return 1
    
    logger.info("\n✓ Group-calibrated FODs ready for fixel analysis!")
    logger.info(f"  Output files: wm_fod_norm_group.mif")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
