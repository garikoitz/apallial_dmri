#!/usr/bin/env python3
"""
Batch Processing Script for FLIP1 Pipeline

This script processes all subjects through the complete flip1 pipeline:
1. Create folder structure in mrtrix_preproc_flip1
2. Flip T2w images
3. Flip manual ROIs (AP/APTR only)
4. Preprocess DWI (convert, flip, denoise, bias correct, register to T2w)
5. Generate FOD
6. Resample ROIs to FOD space (AP/APTR only)
7. Create graft segments (AP/APTR only)
8. Run tractography with MATLAB cleaning (AP/APTR only)

Usage:
    # Process a single subject:
    python functions/batch_process_flip1.py --subject sub-Apallial20R1F
    
    # Process all APTR subjects:
    python functions/batch_process_flip1.py --group APTR
    
    # Process all AP subjects:
    python functions/batch_process_flip1.py --group AP
    
    # Process all CTRL subjects:
    python functions/batch_process_flip1.py --group CTRL
    
    # Process all AP and APTR subjects:
    python functions/batch_process_flip1.py --group AP APTR
    
    # Dry run (show what would be done):
    python functions/batch_process_flip1.py --group APTR --dry-run
    
    # Run specific steps only:
    python functions/batch_process_flip1.py --group APTR --steps 1,2,3

    # Force reprocessing (ignore existing files):
    python functions/batch_process_flip1.py --subject sub-Apallial20R1F --force

The script reads all parameters from config1.yaml and logs progress.
"""

import os
import sys
import subprocess
import argparse
import time
from datetime import datetime

# Set pipeline to flip1
os.environ['MBRAIN_PIPELINE'] = 'flip1'

# Add functions directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config_loader import load_config

# Load config
config = load_config()


def log_message(msg, level="INFO"):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {msg}")


def run_cmd(cmd, dry_run=False, check=True, capture=False):
    """Run a shell command with logging."""
    if dry_run:
        print(f"  [DRY RUN] Would run: {cmd[:100]}..." if len(cmd) > 100 else f"  [DRY RUN] Would run: {cmd}")
        return True
    
    try:
        if capture:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=check)
            return result.stdout.strip() if result.returncode == 0 else None
        else:
            result = subprocess.run(cmd, shell=True, check=check)
            return result.returncode == 0
    except subprocess.CalledProcessError as e:
        log_message(f"Command failed: {e}", "ERROR")
        return False


def file_exists(path):
    """Check if file exists and has non-zero size."""
    return os.path.exists(path) and os.path.getsize(path) > 0


def step1_create_folders(subject, session, derivatives_dir, group, dry_run=False):
    """Step 1: Create folder structure."""
    log_message(f"Step 1: Creating folder structure")
    
    base_dir = os.path.join(derivatives_dir, subject, session)
    folders = [
        os.path.join(base_dir, "dwi"),
    ]
    
    if group != "CTRL":
        folders.extend([
            os.path.join(base_dir, "manual_roi"),
            os.path.join(base_dir, "manual_roi", "fod_space"),
        ])
    
    for folder in folders:
        if not os.path.exists(folder):
            if dry_run:
                print(f"  [DRY RUN] Would create: {folder}")
            else:
                os.makedirs(folder, exist_ok=True)
                print(f"  Created: {folder}")
        else:
            print(f"  Exists: {folder}")
    
    return True


def step2_flip_t2w(subject, session, bids_dir, derivatives_dir, dry_run=False, force=False):
    """Step 2: Flip T2w image using flip1 method (axes 2,1,0 + flip 1)."""
    log_message(f"Step 2: Flipping T2w image")
    
    # Find T2w in BIDS
    bids_anat = os.path.join(bids_dir, subject, session, "anat")
    t2w_input = os.path.join(bids_anat, f"{subject}_{session}_acq-highres_T2w.nii.gz")
    
    if not file_exists(t2w_input):
        log_message(f"T2w not found: {t2w_input}", "ERROR")
        return False
    
    dwi_dir = os.path.join(derivatives_dir, subject, session, "dwi")
    t2w_output = os.path.join(dwi_dir, "t2w_highres_flip.mif")
    
    if file_exists(t2w_output) and not force:
        print(f"  Already exists: {t2w_output}")
        return True
    
    # Flip1: axes 2,1,0 permutation + flip axis 1
    tmp_perm = os.path.join(dwi_dir, "tmp_t2w_perm.mif")
    
    cmds = [
        f"mrconvert '{t2w_input}' '{tmp_perm}' -axes 2,1,0 -strides 1,2,3 -force -quiet",
        f"mrtransform '{tmp_perm}' '{t2w_output}' -flip 1 -force -quiet",
        f"rm -f '{tmp_perm}'"
    ]
    
    for cmd in cmds:
        if not run_cmd(cmd, dry_run):
            return False
    
    if not dry_run:
        print(f"  Created: {t2w_output}")
    
    return True


def step3_flip_rois(subject, session, old_derivatives_dir, derivatives_dir, dry_run=False, force=False):
    """Step 3: Flip manual ROIs from old pipeline to flip1."""
    log_message(f"Step 3: Flipping manual ROIs")
    
    # Source: old pipeline manual_roi_flip (has L/R correction)
    old_roi_dir = os.path.join(old_derivatives_dir, subject, session, "manual_roi_flip")
    new_roi_dir = os.path.join(derivatives_dir, subject, session, "manual_roi")
    
    if not os.path.exists(old_roi_dir):
        log_message(f"Old ROI folder not found: {old_roi_dir}", "ERROR")
        return False
    
    # ROI names
    rois = ["Graft", "Striatum", "Thalamus", "Palecortex", "cerebellum"]
    hemispheres = ["left", "right"]
    
    success_count = 0
    for roi in rois:
        for hemi in hemispheres:
            input_roi = os.path.join(old_roi_dir, f"roi_{roi}_{hemi}.nii.gz")
            output_roi = os.path.join(new_roi_dir, f"roi_{roi}_{hemi}_flip.mif")
            
            if not file_exists(input_roi):
                print(f"  Not found (skipping): {input_roi}")
                continue
            
            if file_exists(output_roi) and not force:
                print(f"  Exists: roi_{roi}_{hemi}_flip.mif")
                success_count += 1
                continue
            
            # Flip1: axes 2,1,0 permutation + flip axis 1
            tmp_perm = os.path.join(new_roi_dir, f"tmp_{roi}_{hemi}_perm.mif")
            
            cmds = [
                f"mrconvert '{input_roi}' '{tmp_perm}' -axes 2,1,0 -strides 1,2,3 -force -quiet",
                f"mrtransform '{tmp_perm}' '{output_roi}' -flip 1 -force -quiet",
                f"rm -f '{tmp_perm}'"
            ]
            
            for cmd in cmds:
                if not run_cmd(cmd, dry_run):
                    continue
            
            if not dry_run:
                print(f"  Created: roi_{roi}_{hemi}_flip.mif")
            success_count += 1
    
    return success_count > 0


def step4_preprocess_dwi(subject, session, bids_dir, derivatives_dir, dry_run=False, force=False):
    """Step 4: Full DWI preprocessing (convert, flip, denoise, bias correct, register to T2w)."""
    log_message(f"Step 4: Preprocessing DWI")
    
    dwi_dir = os.path.join(derivatives_dir, subject, session, "dwi")
    bids_dwi = os.path.join(bids_dir, subject, session, "dwi")
    
    # Find DWI files
    dwi_nii = os.path.join(bids_dwi, f"{subject}_{session}_dwi.nii.gz")
    bvec = os.path.join(bids_dwi, f"{subject}_{session}_dwi.bvec")
    bval = os.path.join(bids_dwi, f"{subject}_{session}_dwi.bval")
    
    if not all(file_exists(f) for f in [dwi_nii, bvec, bval]):
        log_message(f"DWI files not found in BIDS", "ERROR")
        return False
    
    # Final output
    final_dwi = os.path.join(dwi_dir, "dwi_flip_gradcheck_denoise_bias_anat.mif")
    
    if file_exists(final_dwi) and not force:
        print(f"  Already exists: dwi_flip_gradcheck_denoise_bias_anat.mif")
        return True
    
    # Step 4a: Convert to MIF
    dwi_mif = os.path.join(dwi_dir, "dwi.mif")
    if not file_exists(dwi_mif) or force:
        cmd = f"mrconvert '{dwi_nii}' '{dwi_mif}' -fslgrad '{bvec}' '{bval}' -force -quiet"
        if not run_cmd(cmd, dry_run):
            return False
        print(f"  Created: dwi.mif")
    
    # Step 4b: Axis permutation (2,1,0,3 for 4D)
    dwi_perm = os.path.join(dwi_dir, "dwi_perm.mif")
    cmd = f"mrconvert '{dwi_mif}' '{dwi_perm}' -axes 2,1,0,3 -strides 1,2,3,4 -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    
    # Step 4c: Flip axis 1 (Y axis) - NO FLIP 2!
    dwi_flip = os.path.join(dwi_dir, "dwi_flip.mif")
    cmd = f"mrtransform '{dwi_perm}' '{dwi_flip}' -flip 1 -reorient_fod no -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    run_cmd(f"rm -f '{dwi_perm}'", dry_run)
    print(f"  Created: dwi_flip.mif")
    
    # Step 4d: Gradient check
    corrected_b = os.path.join(dwi_dir, "corrected.b")
    dwi_gradcheck = os.path.join(dwi_dir, "dwi_flip_gradcheck.mif")
    cmd = f"dwigradcheck '{dwi_flip}' -export_grad_mrtrix '{corrected_b}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    cmd = f"mrconvert '{dwi_flip}' '{dwi_gradcheck}' -grad '{corrected_b}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: dwi_flip_gradcheck.mif")
    
    # Step 4e: Denoise (MP-PCA)
    noise_mif = os.path.join(dwi_dir, "noise.mif")
    dwi_denoise = os.path.join(dwi_dir, "dwi_flip_gradcheck_denoise.mif")
    cmd = f"dwidenoise -extent 5,5,5 -noise '{noise_mif}' -estimator Exp2 '{dwi_gradcheck}' '{dwi_denoise}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: dwi_flip_gradcheck_denoise.mif")
    
    # Step 4f: Create mask
    mask_mif = os.path.join(dwi_dir, "mask.mif")
    cmd = f"dwi2mask '{dwi_denoise}' '{mask_mif}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: mask.mif")
    
    # Step 4g: Bias correction (ANTs)
    dwi_bias = os.path.join(dwi_dir, "dwi_flip_gradcheck_denoise_bias.mif")
    tmp_bias = os.path.join(dwi_dir, "tmp_bias")
    cmd = f"dwibiascorrect ants -ants.b '[200,3]' -ants.c '[10,1e-6]' -ants.s 4 -mask '{mask_mif}' '{dwi_denoise}' '{dwi_bias}' -scratch '{tmp_bias}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    run_cmd(f"rm -rf '{tmp_bias}'", dry_run)
    print(f"  Created: dwi_flip_gradcheck_denoise_bias.mif")
    
    # Step 4h: Extract mean low-b for registration (use b<50 as "b0" for mouse data)
    b0_bias = os.path.join(dwi_dir, "b0_bias.mif")
    cmd = f"dwiextract '{dwi_bias}' - -bvalue 50 | mrmath - mean '{b0_bias}' -axis 3 -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    
    # Step 4i: Convert to NIfTI for ANTs + apply mask
    t2w_flip = os.path.join(dwi_dir, "t2w_highres_flip.mif")
    t2w_nii = os.path.join(dwi_dir, "t2w_highres_flip.nii.gz")
    b0_brain_mif = os.path.join(dwi_dir, "b0_bias_brain.mif")
    b0_brain_nii = os.path.join(dwi_dir, "b0_bias_brain.nii.gz")
    
    cmds = [
        f"mrconvert '{t2w_flip}' '{t2w_nii}' -force -quiet",
        f"mrcalc '{b0_bias}' '{mask_mif}' -mult '{b0_brain_mif}' -force -quiet",
        f"mrconvert '{b0_brain_mif}' '{b0_brain_nii}' -force -quiet"
    ]
    for cmd in cmds:
        if not run_cmd(cmd, dry_run):
            return False
    
    # Step 4j: ANTs DWI→T2w rigid registration
    ants_prefix = os.path.join(dwi_dir, "dwi2anat_")
    cmd = f"antsRegistrationSyN.sh -d 3 -t r -f '{t2w_nii}' -m '{b0_brain_nii}' -o '{ants_prefix}'"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: ANTs registration files")
    
    # Step 4k: Convert ANTs transform to MRtrix format
    ants_mat = os.path.join(dwi_dir, "dwi2anat_0GenericAffine.mat")
    ants_txt = os.path.join(dwi_dir, "dwi2anat_affine.txt")
    mrtrix_txt = os.path.join(dwi_dir, "dwi2anat_mrtrix.txt")
    
    cmds = [
        f"ConvertTransformFile 3 '{ants_mat}' '{ants_txt}'",
        f"transformconvert '{ants_txt}' itk_import '{mrtrix_txt}' -force -quiet"
    ]
    for cmd in cmds:
        if not run_cmd(cmd, dry_run):
            return False
    
    # Step 4l: Apply transform to full DWI
    cmd = f"mrtransform '{dwi_bias}' '{final_dwi}' -linear '{mrtrix_txt}' -reorient_fod no -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: dwi_flip_gradcheck_denoise_bias_anat.mif")
    
    # Step 4m: Extract final low-b reference
    b0_final = os.path.join(dwi_dir, "b0_final.mif")
    cmd = f"dwiextract '{final_dwi}' - -bvalue 50 | mrmath - mean '{b0_final}' -axis 3 -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: b0_final.mif")
    
    return True


def step5_generate_fod(subject, session, derivatives_dir, dry_run=False, force=False):
    """Step 5: Generate FOD (response function estimation + FOD computation)."""
    log_message(f"Step 5: Generating FOD")
    
    dwi_dir = os.path.join(derivatives_dir, subject, session, "dwi")
    
    # Input
    dwi_anat = os.path.join(dwi_dir, "dwi_flip_gradcheck_denoise_bias_anat.mif")
    mask_mif = os.path.join(dwi_dir, "mask.mif")
    
    if not file_exists(dwi_anat):
        log_message(f"Preprocessed DWI not found: {dwi_anat}", "ERROR")
        return False
    
    # Final output
    wm_fod_norm = os.path.join(dwi_dir, "wm_fod_norm.mif")
    
    if file_exists(wm_fod_norm) and not force:
        print(f"  Already exists: wm_fod_norm.mif")
        return True
    
    # Step 5a: Response function estimation
    response_wm = os.path.join(dwi_dir, "response_wm.txt")
    response_gm = os.path.join(dwi_dir, "response_gm.txt")
    response_csf = os.path.join(dwi_dir, "response_csf.txt")
    response_voxels = os.path.join(dwi_dir, "response_voxels.mif")
    
    cmd = f"dwi2response dhollander '{dwi_anat}' '{response_wm}' '{response_gm}' '{response_csf}' -voxels '{response_voxels}' -mask '{mask_mif}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: response functions")
    
    # Step 5b: FOD estimation (msmt_csd)
    wm_fod = os.path.join(dwi_dir, "wm_fod.mif")
    gm_fod = os.path.join(dwi_dir, "gm_fod.mif")
    csf_fod = os.path.join(dwi_dir, "csf_fod.mif")
    
    cmd = f"dwi2fod msmt_csd '{dwi_anat}' '{response_wm}' '{wm_fod}' '{response_gm}' '{gm_fod}' '{response_csf}' '{csf_fod}' -mask '{mask_mif}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: FOD images")
    
    # Step 5c: FOD normalization (WM + CSF only, skip GM for mouse brain)
    csf_fod_norm = os.path.join(dwi_dir, "csf_fod_norm.mif")
    cmd = f"mtnormalise '{wm_fod}' '{wm_fod_norm}' '{csf_fod}' '{csf_fod_norm}' -mask '{mask_mif}' -force -quiet"
    if not run_cmd(cmd, dry_run):
        return False
    print(f"  Created: wm_fod_norm.mif (normalized)")
    
    return True


def step6_resample_rois(subject, session, derivatives_dir, dry_run=False, force=False):
    """Step 6: Resample manual ROIs to FOD space."""
    log_message(f"Step 6: Resampling ROIs to FOD space")
    
    base_dir = os.path.join(derivatives_dir, subject, session)
    dwi_dir = os.path.join(base_dir, "dwi")
    roi_dir = os.path.join(base_dir, "manual_roi")
    fod_space_dir = os.path.join(roi_dir, "fod_space")
    
    # Template for resampling
    fod_template = os.path.join(dwi_dir, "wm_fod_norm.mif")
    
    if not file_exists(fod_template):
        log_message(f"FOD template not found: {fod_template}", "ERROR")
        return False
    
    os.makedirs(fod_space_dir, exist_ok=True)
    
    # Find all ROI files
    import glob
    roi_files = glob.glob(os.path.join(roi_dir, "roi_*_flip.mif"))
    
    if not roi_files:
        log_message(f"No ROI files found in: {roi_dir}", "ERROR")
        return False
    
    success_count = 0
    for roi_path in roi_files:
        fname = os.path.basename(roi_path).replace(".mif", "_fod.mif")
        output_path = os.path.join(fod_space_dir, fname)
        
        if file_exists(output_path) and not force:
            print(f"  Exists: {fname}")
            success_count += 1
            continue
        
        cmd = f"mrgrid '{roi_path}' regrid '{output_path}' -template '{fod_template}' -interp nearest -force -quiet"
        if run_cmd(cmd, dry_run):
            if not dry_run:
                print(f"  Created: {fname}")
            success_count += 1
    
    return success_count > 0


def step7_create_segments(subject, session, derivatives_dir, dry_run=False, force=False):
    """Step 7: Create graft segments using segment_graft_slicewise.py."""
    log_message(f"Step 7: Creating graft segments (left hemisphere only, seg9)")
    
    dwi_dir = os.path.join(derivatives_dir, subject, session, "dwi")
    seg_dir = os.path.join(dwi_dir, "graft_segments_fod")
    
    # Check if left anterior_dorsal segment already exists
    target_seg = os.path.join(seg_dir, "seg9_graft_left_anterior_dorsal_fod.mif")
    if file_exists(target_seg) and not force:
        print(f"  Already exists: seg9_graft_left_anterior_dorsal_fod.mif")
        return True
    
    # Run segment_graft_slicewise.py for LEFT hemisphere only, 9 segments
    script_path = os.path.join(os.path.dirname(__file__), "segment_graft_slicewise.py")
    
    cmd = f"python '{script_path}' --subject {subject} --session {session} --segments 9 --hemisphere left"
    if not run_cmd(cmd, dry_run):
        return False
    
    return True


def step8_tractography(subject, session, derivatives_dir, dry_run=False, force=False):
    """Step 8: Run tractography with MATLAB cleaning (left anterior_dorsal only)."""
    log_message(f"Step 8: Running tractography (left anterior_dorsal → Thalamus + Striatum)")
    
    # Run generate_tracts.py for LEFT anterior_dorsal only
    script_path = os.path.join(os.path.dirname(__file__), "generate_tracts.py")
    
    cmd = f"python '{script_path}' --subject {subject} --seg-mode 9 --segment anterior_dorsal --hemisphere left --clean"
    if force:
        cmd += " --force"
    
    if not run_cmd(cmd, dry_run):
        return False
    
    return True


def process_subject(subject, session, config, group, steps="1,2,3,4,5,6,7,8", dry_run=False, force=False):
    """Process a single subject through all specified steps."""
    
    log_message(f"{'='*70}")
    log_message(f"Processing: {subject}/{session}")
    log_message(f"{'='*70}")
    
    bids_dir = config['paths']['bids_dir']
    derivatives_dir = config['paths']['derivatives_dir']
    old_derivatives_dir = config['paths']['derivatives_dir_old']
    
    # Parse steps
    step_list = [int(s) for s in steps.split(",")]
    
    results = {}
    
    if 1 in step_list:
        results[1] = step1_create_folders(subject, session, derivatives_dir, group, dry_run)
    
    if 2 in step_list:
        results[2] = step2_flip_t2w(subject, session, bids_dir, derivatives_dir, dry_run, force)
        if not results.get(2, True):
            log_message("Step 2 failed, stopping", "ERROR")
            return results
    
    if 3 in step_list and group != "CTRL":
        results[3] = step3_flip_rois(subject, session, old_derivatives_dir, derivatives_dir, dry_run, force)
        if not results.get(3, True):
            log_message("Step 3 failed, stopping", "ERROR")
            return results
    
    if 4 in step_list:
        results[4] = step4_preprocess_dwi(subject, session, bids_dir, derivatives_dir, dry_run, force)
        if not results.get(4, True):
            log_message("Step 4 failed, stopping", "ERROR")
            return results
    
    if 5 in step_list:
        results[5] = step5_generate_fod(subject, session, derivatives_dir, dry_run, force)
        if not results.get(5, True):
            log_message("Step 5 failed, stopping", "ERROR")
            return results
    
    if 6 in step_list and group != "CTRL":
        results[6] = step6_resample_rois(subject, session, derivatives_dir, dry_run, force)
        if not results.get(6, True):
            log_message("Step 6 failed, stopping", "ERROR")
            return results
    
    if 7 in step_list and group != "CTRL":
        results[7] = step7_create_segments(subject, session, derivatives_dir, dry_run, force)
        if not results.get(7, True):
            log_message("Step 7 failed, stopping", "ERROR")
            return results
    
    if 8 in step_list and group != "CTRL":
        results[8] = step8_tractography(subject, session, derivatives_dir, dry_run, force)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Batch processing for FLIP1 pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--subject", help="Single subject ID to process")
    parser.add_argument("--group", nargs="+", choices=["AP", "APTR", "CTRL"], 
                        help="Process all subjects in group(s)")
    parser.add_argument("--steps", default="1,2,3,4,5,6,7,8",
                        help="Comma-separated list of steps to run (default: all)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without executing")
    parser.add_argument("--force", action="store_true",
                        help="Force reprocessing even if files exist")
    args = parser.parse_args()
    
    if not args.subject and not args.group:
        parser.error("Must specify --subject or --group")
    
    # Build list of subjects to process
    subjects_to_process = []
    
    if args.subject:
        # Find session for this subject
        for group_name in ["AP", "APTR", "CTRL"]:
            for subj_info in config['groups'].get(group_name, []):
                if subj_info['id'] == args.subject:
                    subjects_to_process.append((args.subject, subj_info['session'], group_name))
                    break
    
    if args.group:
        for group_name in args.group:
            for subj_info in config['groups'].get(group_name, []):
                # Skip if already processed sub-Apallial17NmF (reference subject)
                if subj_info['id'] == 'sub-Apallial17NmF':
                    log_message(f"Skipping reference subject: {subj_info['id']}")
                    continue
                subjects_to_process.append((subj_info['id'], subj_info['session'], group_name))
    
    if not subjects_to_process:
        log_message("No subjects found to process", "ERROR")
        return 1
    
    log_message(f"{'='*70}")
    log_message(f"FLIP1 BATCH PROCESSING")
    log_message(f"{'='*70}")
    log_message(f"Subjects: {len(subjects_to_process)}")
    log_message(f"Steps: {args.steps}")
    log_message(f"Dry run: {args.dry_run}")
    log_message(f"Force: {args.force}")
    
    # Process each subject
    all_results = {}
    start_time = time.time()
    
    for i, (subject, session, group) in enumerate(subjects_to_process, 1):
        log_message(f"\n[{i}/{len(subjects_to_process)}] Subject: {subject} (Group: {group})")
        
        results = process_subject(
            subject, session, config, group,
            steps=args.steps,
            dry_run=args.dry_run,
            force=args.force
        )
        all_results[subject] = results
    
    # Summary
    elapsed = time.time() - start_time
    log_message(f"\n{'='*70}")
    log_message(f"PROCESSING COMPLETE")
    log_message(f"{'='*70}")
    log_message(f"Total time: {elapsed/60:.1f} minutes")
    
    for subject, results in all_results.items():
        status = "✓" if all(results.values()) else "✗"
        log_message(f"  {status} {subject}: steps {list(results.keys())}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
