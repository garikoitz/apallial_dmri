"""
MRtrix3 preprocessing pipeline for mouse brain DTI data.

This module follows the structure from rtp2-preproc:
- Each step is an independent function that takes the current filename (difm/t2wm)
- Each function returns the updated filename with new suffix
- A main function calls steps based on True/False configuration
- The naming accumulates based on which steps are run
- Steps can check if output exists and skip accordingly
- Uses -mask when available (following rtp2-preproc pattern)

Order of steps:
DWI Preprocessing:
    1. convert      - Convert NIfTI to MIF with gradient info
    2. flip         - Axis permutation and flip for orientation correction  
    3. gradcheck    - Check gradient orientation
    4. denoise      - PCA denoising (Marchenko-Pastur)
    5. degibbs      - Gibbs ringing artifact removal
    6. mask         - Brain mask estimation
    7. bias         - Bias field correction (ANTs or FSL)
    8. norm         - Intensity normalization (optional)

Anatomical T2w Processing:
    9a. convert_t2w (lowres)  - Convert lowres T2w to MIF
    9b. flip_t2w (lowres)     - Flip lowres T2w for orientation
    9c. convert_t2w (highres) - Convert highres T2w to MIF  
    9d. flip_t2w (highres)    - Flip highres T2w for orientation

Atlas Coregistration:
    10. atlas_coreg  - Coregister Allen CCF atlas to T2w
    11. check_rois   - Extract and verify ROI placement

DWI to Anatomical Coregistration (FINAL STEP):
    12. dwi2anat     - Coregister DWI to T2w (uses mask)

Export:
    13. export_nifti - Export final MIF to NIfTI (optional)
"""

import logging
import os
import subprocess as sp
from pathlib import Path
from typing import Optional, List, Dict, Any

from .config import Config

log = logging.getLogger(__name__)


def run_cmd(cmd: str, dry_run: bool = False) -> None:
    """Run a shell command with logging."""
    log.info("Running: %s", cmd)
    if not dry_run:
        sp.run(cmd, shell=True, check=True)


def do_convert(input_nii: str, bvec: str, bval: str, output_dir: str, 
               common: str = "-force", dry_run: bool = False) -> str:
    """
    Convert NIfTI to MIF format with gradient information.
    
    Arguments:
        input_nii: Path to input NIfTI file
        bvec: Path to bvec file  
        bval: Path to bval file
        output_dir: Output directory
        common: Common MRtrix options (default: "-force")
        dry_run: If True, only show commands
        
    Returns:
        str: Base filename for subsequent steps (difm)
    """
    # Determine base name from input
    basename = Path(input_nii).name.replace('.nii.gz', '').replace('.nii', '')
    difm = basename
    output_mif = os.path.join(output_dir, f"{difm}.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Convert: Output already exists, skipping: %s", output_mif)
        return difm
    
    cmd = (
        f"mrconvert {input_nii} {output_mif} "
        f"-fslgrad {bvec} {bval} {common}"
    )
    log.info("Converting NIfTI to MIF with:\n '%s'", cmd)
    run_cmd(cmd, dry_run)
    
    return difm


def do_flip(difm: str, output_dir: str, axes: str = "2,1,0,3", 
            flip_axis: str = "1", common: str = "-force", 
            dry_run: bool = False) -> str:
    """
    Apply axis permutation and flip for orientation correction.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        axes: Axis permutation string (default: "2,1,0,3" for 4D)
        flip_axis: Axis to flip (default: "1")
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _flip suffix
    """
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    temp_mif = os.path.join(output_dir, f"{difm}_perm.mif")
    output_mif = os.path.join(output_dir, f"{difm}_flip.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Flip: Output already exists, skipping: %s", output_mif)
        return f"{difm}_flip"
    
    # Step 1: Axis permutation
    cmd1 = (
        f"mrconvert {input_mif} {temp_mif} "
        f"-axes {axes} -strides 1,2,3,4 {common}"
    )
    log.info("Axis permutation with:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Step 2: Axis flip
    cmd2 = (
        f"mrtransform {temp_mif} {output_mif} "
        f"-flip {flip_axis} -reorient_fod no {common}"
    )
    log.info("Axis flip with:\n '%s'", cmd2)
    run_cmd(cmd2, dry_run)
    
    # Clean up temp file
    if os.path.isfile(temp_mif) and not dry_run:
        os.remove(temp_mif)
    
    difm = f"{difm}_flip"
    return difm


def do_gradcheck(difm: str, output_dir: str, common: str = "-force",
                 do_step: bool = True, dry_run: bool = False) -> str:
    """
    Check and correct gradient orientation using dwigradcheck.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        common: Common MRtrix options
        do_step: Whether to actually perform this step (True) or just pass through
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name (unchanged if gradcheck doesn't modify file)
    """
    if not do_step:
        return difm
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    grad_file = os.path.join(output_dir, "corrected.b")
    output_mif = os.path.join(output_dir, f"{difm}_gradcheck.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Gradcheck: Output already exists, skipping: %s", output_mif)
        return f"{difm}_gradcheck"
    
    # Run gradient check
    cmd1 = (
        f"dwigradcheck {input_mif} "
        f"-export_grad_mrtrix {grad_file} {common}"
    )
    log.info("Checking gradient orientation with:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Apply corrected gradients
    cmd2 = (
        f"mrconvert {input_mif} {output_mif} "
        f"-grad {grad_file} {common}"
    )
    log.info("Applying corrected gradients with:\n '%s'", cmd2)
    run_cmd(cmd2, dry_run)
    
    difm = f"{difm}_gradcheck"
    return difm


def do_denoise(difm: str, output_dir: str, common: str = "-force",
               do_step: bool = True, dry_run: bool = False) -> str:
    """
    Perform MRtrix's PCA denoising using Marchenko-Pastur distribution.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _denoise suffix
    """
    if not do_step:
        return difm
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_mif = os.path.join(output_dir, f"{difm}_denoise.mif")
    noise_mif = os.path.join(output_dir, f"{difm}_noise.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Denoise: Output already exists, skipping: %s", output_mif)
        return f"{difm}_denoise"
    
    cmd = (
        f"dwidenoise -extent 5,5,5 -noise {noise_mif} -estimator Exp2 "
        f"{input_mif} {output_mif} {common}"
    )
    log.info("Performing PCA denoising with:\n '%s'", cmd)
    run_cmd(cmd, dry_run)
    
    difm = f"{difm}_denoise"
    return difm


def do_degibbs(difm: str, output_dir: str, common: str = "-force",
               do_step: bool = True, dry_run: bool = False) -> str:
    """
    Perform MRtrix's Gibbs ringing correction.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _degibbs suffix
    """
    if not do_step:
        return difm
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_mif = os.path.join(output_dir, f"{difm}_degibbs.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Degibbs: Output already exists, skipping: %s", output_mif)
        return f"{difm}_degibbs"
    
    cmd = (
        f"mrdegibbs -nshifts 20 -minW 1 -maxW 3 "
        f"{input_mif} {output_mif} {common}"
    )
    log.info("Performing Gibbs ringing correction with:\n '%s'", cmd)
    run_cmd(cmd, dry_run)
    
    difm = f"{difm}_degibbs"
    return difm


def do_mask(difm: str, output_dir: str, mask_name: str = "mask",
            common: str = "-force", do_step: bool = True, 
            dry_run: bool = False) -> str:
    """
    Create brain mask from DWI data.
    
    For mouse ex-vivo data, dwi2mask may fail. This function includes
    a fallback to threshold-based masking.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        mask_name: Base name for mask file (default: "mask")
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        
    Returns:
        str: Mask file name (without extension)
    """
    if not do_step:
        return mask_name
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    mask_mif = os.path.join(output_dir, f"{mask_name}.mif")
    
    # Skip if output exists
    if os.path.isfile(mask_mif) and not dry_run:
        log.info("Mask: Output already exists, skipping: %s", mask_mif)
        return mask_name
    
    # Try dwi2mask first
    cmd = f"dwi2mask {input_mif} {mask_mif} {common}"
    log.info("Creating brain mask with:\n '%s'", cmd)
    
    try:
        run_cmd(cmd, dry_run)
    except sp.CalledProcessError:
        log.warning("dwi2mask failed (expected for ex-vivo/mouse data), using fallback")
        # Fallback: threshold-based mask
        cmd_fallback = f"mrthreshold {input_mif} {mask_mif} -abs 10 {common}"
        log.info("Creating threshold-based mask with:\n '%s'", cmd_fallback)
        run_cmd(cmd_fallback, dry_run)
    
    return mask_name


def do_bias(difm: str, mask: str, output_dir: str, 
            bias_method: str = "ants", common: str = "-force",
            do_step: bool = True, dry_run: bool = False,
            antsb: str = "[200,3]", antsc: str = "[1e-6,10]", 
            antss: str = "4") -> str:
    """
    Run ANTs or FSL bias correction.
    
    Arguments:
        difm: Current diffusion file base name
        mask: Mask file base name
        output_dir: Output directory
        bias_method: Bias correction method ("ants" or "fsl")
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        antsb: ANTs -b parameter
        antsc: ANTs -c parameter
        antss: ANTs -s parameter
        
    Returns:
        str: Updated file name with _bias suffix
    """
    if not do_step:
        return difm
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_mif = os.path.join(output_dir, f"{difm}_bias.mif")
    mask_mif = os.path.join(output_dir, f"{mask}.mif")
    tmp_dir = os.path.join(output_dir, "tmp_bias")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Bias: Output already exists, skipping: %s", output_mif)
        return f"{difm}_bias"
    
    # Create scratch directory
    os.makedirs(tmp_dir, exist_ok=True)
    
    if bias_method == "ants":
        if os.path.isfile(mask_mif):
            cmd = (
                f"dwibiascorrect ants -ants.b {antsb} -ants.c {antsc} "
                f"-ants.s {antss} -mask {mask_mif} "
                f"{input_mif} {output_mif} -scratch {tmp_dir} {common}"
            )
        else:
            cmd = (
                f"dwibiascorrect ants -ants.b {antsb} -ants.c {antsc} "
                f"-ants.s {antss} "
                f"{input_mif} {output_mif} -scratch {tmp_dir} {common}"
            )
        log.info("Performing bias correction with ANTs:\n '%s'", cmd)
    elif bias_method == "fsl":
        if os.path.isfile(mask_mif):
            cmd = (
                f"dwibiascorrect fsl -mask {mask_mif} "
                f"{input_mif} {output_mif} -scratch {tmp_dir} {common}"
            )
        else:
            cmd = (
                f"dwibiascorrect fsl "
                f"{input_mif} {output_mif} -scratch {tmp_dir} {common}"
            )
        log.info("Performing bias correction with FSL:\n '%s'", cmd)
    else:
        raise ValueError(f"Unknown bias method: {bias_method}")
    
    run_cmd(cmd, dry_run)
    
    difm = f"{difm}_bias"
    return difm


def do_norm(difm: str, mask: str, output_dir: str, nval: str = "1000",
            common: str = "-force", do_step: bool = True,
            dry_run: bool = False) -> str:
    """
    Perform intensity normalization of DWI data.
    
    Creates a WM mask based on FA and normalizes intensity to a target value.
    
    Arguments:
        difm: Current diffusion file base name
        mask: Mask file base name
        output_dir: Output directory
        nval: Target intensity value for normalization
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _norm suffix
    """
    if not do_step:
        return difm
        
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_mif = os.path.join(output_dir, f"{difm}_norm.mif")
    mask_mif = os.path.join(output_dir, f"{mask}.mif")
    wm_mif = os.path.join(output_dir, "wm.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("Norm: Output already exists, skipping: %s", output_mif)
        return f"{difm}_norm"
    
    # Create FA-based WM mask
    cmd1 = (
        f"dwi2tensor -mask {mask_mif} {input_mif} - {common} | "
        f"tensor2metric - -fa - {common} | "
        f"mrthreshold -abs 0.5 - {wm_mif} {common}"
    )
    log.info("Creating FA-based WM mask with:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Normalize intensity
    cmd2 = (
        f"dwinormalise individual -intensity {nval} "
        f"{input_mif} {wm_mif} {output_mif} {common}"
    )
    log.info("Performing intensity normalization with:\n '%s'", cmd2)
    run_cmd(cmd2, dry_run)
    
    difm = f"{difm}_norm"
    return difm


def do_export_nifti(difm: str, output_dir: str, common: str = "-force",
                    dry_run: bool = False) -> str:
    """
    Export final MIF to NIfTI format.
    
    Arguments:
        difm: Current diffusion file base name
        output_dir: Output directory
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        str: Path to output NIfTI file
    """
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_nii = os.path.join(output_dir, f"{difm}.nii.gz")
    
    # Skip if output exists
    if os.path.isfile(output_nii) and not dry_run:
        log.info("Export: Output already exists, skipping: %s", output_nii)
        return output_nii
    
    cmd = f"mrconvert {input_mif} {output_nii} {common}"
    log.info("Exporting to NIfTI with:\n '%s'", cmd)
    run_cmd(cmd, dry_run)
    
    return output_nii


def create_ref_images(difm: str, mask: str, output_dir: str, 
                      common: str = "-force", dry_run: bool = False) -> None:
    """
    Create b0 reference images.
    
    Arguments:
        difm: Current diffusion file base name
        mask: Mask file base name
        output_dir: Output directory
        common: Common MRtrix options
        dry_run: If True, only show commands
    """
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    mask_mif = os.path.join(output_dir, f"{mask}.mif")
    b0_mif = os.path.join(output_dir, "b0_dwi.mif")
    b0_nii = os.path.join(output_dir, "b0_dwi.nii.gz")
    mask_nii = os.path.join(output_dir, f"{mask}.nii.gz")
    
    log.info("Creating DWI space b0 reference images...")
    
    # Extract mean b0
    cmd1 = (
        f"dwiextract {input_mif} - -bzero {common} | "
        f"mrmath - mean {b0_mif} -axis 3 {common}"
    )
    log.info("Extracting mean b0 with:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Convert to NIfTI
    cmd2 = f"mrconvert {b0_mif} {b0_nii} {common}"
    cmd3 = f"mrconvert {mask_mif} {mask_nii} {common}"
    log.info("Converting to NIfTI:\n '%s'\n '%s'", cmd2, cmd3)
    run_cmd(cmd2, dry_run)
    run_cmd(cmd3, dry_run)


# =============================================================================
# ANATOMICAL T2w PROCESSING FUNCTIONS
# =============================================================================

def do_convert_t2w(input_nii: str, output_dir: str, acq_type: str = "lowres",
                   common: str = "-force", dry_run: bool = False) -> str:
    """
    Convert T2w anatomical NIfTI to MIF format.
    
    Arguments:
        input_nii: Path to input T2w NIfTI file
        output_dir: Output directory
        acq_type: Acquisition type - "lowres" or "highres"
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        str: Base filename for subsequent steps (t2wm)
    """
    t2wm = f"t2w_{acq_type}"
    output_mif = os.path.join(output_dir, f"{t2wm}.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("T2w Convert (%s): Output already exists, skipping: %s", acq_type, output_mif)
        return t2wm
    
    cmd = f"mrconvert {input_nii} {output_mif} {common}"
    log.info("Converting T2w %s to MIF with:\n '%s'", acq_type, cmd)
    run_cmd(cmd, dry_run)
    
    return t2wm


def do_flip_t2w(t2wm: str, output_dir: str, acq_type: str = "lowres",
                axes: str = "2,1,0", flip_axis: str = "1", 
                common: str = "-force", dry_run: bool = False) -> str:
    """
    Apply axis permutation and flip for T2w orientation correction.
    
    Arguments:
        t2wm: Current T2w file base name
        output_dir: Output directory
        acq_type: Acquisition type - "lowres" or "highres"
        axes: Axis permutation string (default: "2,1,0" for 3D)
        flip_axis: Axis to flip (default: "1")
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _flip suffix
    """
    input_mif = os.path.join(output_dir, f"{t2wm}.mif")
    temp_mif = os.path.join(output_dir, f"{t2wm}_perm.mif")
    output_mif = os.path.join(output_dir, f"{t2wm}_flip.mif")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("T2w Flip (%s): Output already exists, skipping: %s", acq_type, output_mif)
        return f"{t2wm}_flip"
    
    # Step 1: Axis permutation
    cmd1 = (
        f"mrconvert {input_mif} {temp_mif} "
        f"-axes {axes} -strides 1,2,3 {common}"
    )
    log.info("T2w %s axis permutation with:\n '%s'", acq_type, cmd1)
    run_cmd(cmd1, dry_run)
    
    # Step 2: Axis flip
    cmd2 = (
        f"mrtransform {temp_mif} {output_mif} "
        f"-flip {flip_axis} {common}"
    )
    log.info("T2w %s axis flip with:\n '%s'", acq_type, cmd2)
    run_cmd(cmd2, dry_run)
    
    # Clean up temp file
    if os.path.isfile(temp_mif) and not dry_run:
        os.remove(temp_mif)
    
    t2wm = f"{t2wm}_flip"
    return t2wm


# =============================================================================
# ATLAS TO T2w COREGISTRATION
# =============================================================================

def do_atlas_coreg(t2wm: str, output_dir: str, atlas_template: str,
                   atlas_annotation: str, ants_options: str = "-d 3 -t s",
                   common: str = "-force", dry_run: bool = False) -> Dict[str, str]:
    """
    Coregister Allen Mouse Brain CCF atlas to subject T2w image.
    
    Uses ANTs SyN registration to warp atlas template and annotation
    to subject space. Creates transforms for mapping ROIs.
    
    Arguments:
        t2wm: T2w file base name (target for registration)
        output_dir: Output directory
        atlas_template: Path to atlas template (Allen CCF average template)
        atlas_annotation: Path to atlas annotation (Allen CCF annotation)
        ants_options: ANTs registration options (default: "-d 3 -t s" for SyN)
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        dict: Dictionary with paths to warped atlas files and transforms
    """
    t2w_mif = os.path.join(output_dir, f"{t2wm}.mif")
    t2w_nii = os.path.join(output_dir, f"{t2wm}.nii.gz")
    atlas_output_dir = os.path.join(output_dir, "atlas")
    os.makedirs(atlas_output_dir, exist_ok=True)
    
    # Output files
    warped_template = os.path.join(atlas_output_dir, "atlas_template_warped.nii.gz")
    warped_annotation = os.path.join(atlas_output_dir, "atlas_annotation_warped.nii.gz")
    atlas2t2w_prefix = os.path.join(atlas_output_dir, "atlas2t2w_")
    
    results = {
        'warped_template': warped_template,
        'warped_annotation': warped_annotation,
        'transform_prefix': atlas2t2w_prefix,
    }
    
    # Skip if outputs exist
    if os.path.isfile(warped_template) and os.path.isfile(warped_annotation) and not dry_run:
        log.info("Atlas coreg: Outputs already exist, skipping")
        return results
    
    log.info("Performing atlas to T2w coregistration...")
    
    # Convert T2w to NIfTI for ANTs
    cmd1 = f"mrconvert {t2w_mif} {t2w_nii} {common}"
    log.info("Converting T2w to NIfTI:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Run ANTs registration: atlas template → T2w
    # Fixed: T2w (subject space), Moving: atlas template
    cmd2 = (
        f"antsRegistrationSyN.sh {ants_options} "
        f"-f {t2w_nii} -m {atlas_template} -o {atlas2t2w_prefix}"
    )
    log.info("ANTs registration (atlas → T2w):\n '%s'", cmd2)
    run_cmd(cmd2, dry_run)
    
    # Apply transform to atlas template
    warp_file = f"{atlas2t2w_prefix}1Warp.nii.gz"
    affine_file = f"{atlas2t2w_prefix}0GenericAffine.mat"
    
    cmd3 = (
        f"antsApplyTransforms -d 3 -i {atlas_template} -r {t2w_nii} "
        f"-t {warp_file} -t {affine_file} -o {warped_template}"
    )
    log.info("Applying transform to atlas template:\n '%s'", cmd3)
    run_cmd(cmd3, dry_run)
    
    # Apply transform to atlas annotation (use NearestNeighbor for labels)
    cmd4 = (
        f"antsApplyTransforms -d 3 -i {atlas_annotation} -r {t2w_nii} "
        f"-t {warp_file} -t {affine_file} -n NearestNeighbor -o {warped_annotation}"
    )
    log.info("Applying transform to atlas annotation:\n '%s'", cmd4)
    run_cmd(cmd4, dry_run)
    
    log.info("✓ Atlas coregistration complete")
    return results


def do_check_rois(atlas_annotation_warped: str, output_dir: str,
                  roi_ids: Optional[Dict[str, int]] = None,
                  common: str = "-force", dry_run: bool = False) -> Dict[str, str]:
    """
    Check that warped ROIs are in the correct location.
    
    Creates visualization files and extracts individual ROI masks.
    
    Arguments:
        atlas_annotation_warped: Path to warped atlas annotation
        output_dir: Output directory
        roi_ids: Dictionary of ROI names to structure IDs (default: cortex, striatum, thalamus)
        common: Common MRtrix options
        dry_run: If True, only show commands
        
    Returns:
        dict: Dictionary with paths to ROI mask files
    """
    if roi_ids is None:
        roi_ids = {
            'cortex': 315,
            'striatum': 477,
            'thalamus': 549,
            'hippocampus': 1089,
        }
    
    atlas_output_dir = os.path.join(output_dir, "atlas", "rois")
    os.makedirs(atlas_output_dir, exist_ok=True)
    
    roi_masks = {}
    
    log.info("Extracting ROI masks from warped atlas...")
    
    for roi_name, roi_id in roi_ids.items():
        roi_mask = os.path.join(atlas_output_dir, f"roi_{roi_name}.nii.gz")
        
        # Skip if exists
        if os.path.isfile(roi_mask) and not dry_run:
            log.info("ROI %s: Output already exists, skipping", roi_name)
            roi_masks[roi_name] = roi_mask
            continue
        
        # Extract ROI using mrcalc
        # Note: For proper ROI extraction with child structures, 
        # use atlas_manager.get_structure_mask instead
        cmd = (
            f"mrcalc {atlas_annotation_warped} {roi_id} -eq "
            f"{roi_mask} {common}"
        )
        log.info("Extracting %s ROI (ID=%d):\n '%s'", roi_name, roi_id, cmd)
        run_cmd(cmd, dry_run)
        
        roi_masks[roi_name] = roi_mask
    
    log.info("✓ ROI extraction complete")
    return roi_masks


# =============================================================================
# DWI TO ANATOMICAL T2w COREGISTRATION (FINAL STEP)
# =============================================================================

def do_dwi2anat(difm: str, mask: str, t2wm: str, output_dir: str,
                ants_options: str = "-d 3 -t r", common: str = "-force",
                do_step: bool = True, dry_run: bool = False) -> str:
    """
    Coregister DWI to anatomical T2w image.
    
    This should be the LAST step in the DWI preprocessing pipeline.
    Uses ANTs registration to align the DWI b0 to the T2w anatomical,
    then applies the transform to the full DWI data.
    
    Following rtp2-preproc pattern: uses mask when available.
    
    Arguments:
        difm: Current diffusion file base name
        mask: Mask file base name
        t2wm: T2w anatomical file base name (target)
        output_dir: Output directory
        ants_options: ANTs registration options (default: "-d 3 -t r" for rigid)
        common: Common MRtrix options
        do_step: Whether to actually perform this step
        dry_run: If True, only show commands
        
    Returns:
        str: Updated file name with _anat suffix
    """
    if not do_step:
        return difm
    
    input_mif = os.path.join(output_dir, f"{difm}.mif")
    output_mif = os.path.join(output_dir, f"{difm}_anat.mif")
    mask_mif = os.path.join(output_dir, f"{mask}.mif")
    t2w_mif = os.path.join(output_dir, f"{t2wm}.mif")
    
    # NIfTI versions for ANTs
    b0_nii = os.path.join(output_dir, "b0_dwi_for_coreg.nii.gz")
    b0_brain_nii = os.path.join(output_dir, "b0_dwi_brain_for_coreg.nii.gz")
    mask_nii = os.path.join(output_dir, f"{mask}_for_coreg.nii.gz")
    t2w_nii = os.path.join(output_dir, f"{t2wm}_for_coreg.nii.gz")
    
    # Transform files
    transform_prefix = os.path.join(output_dir, "dwi2anat_")
    transform_mat = f"{transform_prefix}0GenericAffine.mat"
    transform_txt = os.path.join(output_dir, "dwi2anat_affine.txt")
    transform_mrtrix = os.path.join(output_dir, "dwi2anat_mrtrix.txt")
    
    # Skip if output exists
    if os.path.isfile(output_mif) and not dry_run:
        log.info("DWI2Anat: Output already exists, skipping: %s", output_mif)
        return f"{difm}_anat"
    
    # Check that T2w exists
    if not os.path.isfile(t2w_mif) and not dry_run:
        log.warning("DWI2Anat: T2w file not found: %s, skipping coregistration", t2w_mif)
        return difm
    
    log.info("Performing DWI to T2w coregistration...")
    
    # Extract mean b0 from DWI
    b0_mif = os.path.join(output_dir, "b0_dwi_for_coreg.mif")
    cmd1 = (
        f"dwiextract {input_mif} - -bzero {common} | "
        f"mrmath - mean {b0_mif} -axis 3 {common}"
    )
    log.info("Extracting mean b0:\n '%s'", cmd1)
    run_cmd(cmd1, dry_run)
    
    # Convert b0 to NIfTI
    cmd2 = f"mrconvert {b0_mif} {b0_nii} {common}"
    log.info("Converting b0 to NIfTI:\n '%s'", cmd2)
    run_cmd(cmd2, dry_run)
    
    # Convert T2w to NIfTI
    cmd3 = f"mrconvert {t2w_mif} {t2w_nii} {common}"
    log.info("Converting T2w to NIfTI:\n '%s'", cmd3)
    run_cmd(cmd3, dry_run)
    
    # Apply mask to b0 if available (following rtp2-preproc pattern)
    if os.path.isfile(mask_mif) or dry_run:
        cmd4a = f"mrconvert {mask_mif} {mask_nii} {common}"
        log.info("Converting mask to NIfTI:\n '%s'", cmd4a)
        run_cmd(cmd4a, dry_run)
        
        cmd4b = f"fslmaths {b0_nii} -mas {mask_nii} {b0_brain_nii}"
        log.info("Applying mask to b0:\n '%s'", cmd4b)
        run_cmd(cmd4b, dry_run)
        
        moving_image = b0_brain_nii
    else:
        log.warning("Mask not found, using unmasked b0 for registration")
        moving_image = b0_nii
    
    # ANTs registration: DWI b0 → T2w
    cmd5 = (
        f"antsRegistrationSyN.sh {ants_options} "
        f"-f {t2w_nii} -m {moving_image} -o {transform_prefix}"
    )
    log.info("ANTs registration (DWI → T2w):\n '%s'", cmd5)
    run_cmd(cmd5, dry_run)
    
    # Convert ANTs transform to MRtrix format
    cmd6 = f"ConvertTransformFile 3 {transform_mat} {transform_txt}"
    log.info("Converting transform:\n '%s'", cmd6)
    run_cmd(cmd6, dry_run)
    
    cmd7 = f"transformconvert {transform_txt} itk_import {transform_mrtrix} -force"
    log.info("Importing to MRtrix:\n '%s'", cmd7)
    run_cmd(cmd7, dry_run)
    
    # Apply transform to full DWI
    # Note: -reorient_fod no for raw DWI (not FOD)
    cmd8 = (
        f"mrtransform {input_mif} {output_mif} "
        f"-linear {transform_mrtrix} -reorient_fod no {common}"
    )
    log.info("Applying transform to DWI:\n '%s'", cmd8)
    run_cmd(cmd8, dry_run)
    
    log.info("✓ DWI to T2w coregistration complete")
    
    difm = f"{difm}_anat"
    return difm


class DWIPreprocessor:
    """
    DWI preprocessing pipeline for mouse brain data.
    
    This class coordinates the preprocessing steps, calling independent
    functions based on configuration flags.
    """
    
    def __init__(self, bids_dir: Optional[Path] = None, 
                 derivatives_dir: Optional[Path] = None):
        """
        Initialize the DWI preprocessor.
        
        Parameters
        ----------
        bids_dir : Path, optional
            Path to BIDS dataset
        derivatives_dir : Path, optional
            Path to derivatives output directory
        """
        self.bids_dir = Path(bids_dir) if bids_dir else Config.BIDS_DIR
        self.derivatives_dir = Path(derivatives_dir) if derivatives_dir else Config.DERIVATIVES_DIR
        
        # Create preprocessing directory
        self.preproc_dir = self.derivatives_dir / "mrtrix_preproc"
        self.preproc_dir.mkdir(parents=True, exist_ok=True)
        
        # Track processing statistics
        self.stats = {
            'subjects_processed': 0,
            'steps_completed': [],
            'errors': [],
            'warnings': [],
            'subjects_with_issues': {},
            'commands_run': [],
        }
        
        # Default step configuration
        self.default_steps = {
            # DWI preprocessing steps
            'convert': True,
            'flip': True,
            'gradcheck': True,
            'denoise': True,
            'degibbs': False,
            'mask': True,
            'bias': True,
            'norm': False,
            # T2w anatomical processing
            'convert_t2w': True,
            'flip_t2w': True,
            # Atlas coregistration
            'atlas_coreg': True,
            'check_rois': True,
            # DWI to T2w coregistration (final step)
            'dwi2anat': True,
            # Export
            'export_nifti': False,
        }
    
    def preprocess_dwi(self, 
                       input_nii: str,
                       bvec: str,
                       bval: str,
                       output_dir: str,
                       subject_id: str = "",
                       session_id: str = "",
                       t2w_lowres: Optional[str] = None,
                       t2w_highres: Optional[str] = None,
                       atlas_template: Optional[str] = None,
                       atlas_annotation: Optional[str] = None,
                       steps: Optional[Dict[str, bool]] = None,
                       dry_run: bool = False,
                       force_redo: bool = False) -> Dict[str, Any]:
        """
        Preprocess a single DWI acquisition.
        
        Parameters
        ----------
        input_nii : str
            Path to input NIfTI file
        bvec : str
            Path to bvec file
        bval : str
            Path to bval file
        output_dir : str
            Output directory for this subject/session
        subject_id : str
            Subject ID for logging
        session_id : str
            Session ID for logging
        t2w_lowres : str, optional
            Path to lowres T2w anatomical
        t2w_highres : str, optional
            Path to highres T2w anatomical
        atlas_template : str, optional
            Path to atlas template (e.g., Allen CCF)
        atlas_annotation : str, optional
            Path to atlas annotation volume
        steps : dict, optional
            Dictionary of step names to True/False values
            If None, uses default_steps
        dry_run : bool
            If True, show commands without executing
        force_redo : bool
            If True, redo steps even if output exists
            
        Returns
        -------
        dict
            Processing results including final filename and steps completed
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Merge with defaults
        step_config = self.default_steps.copy()
        if steps:
            step_config.update(steps)
        
        common = "-force" if force_redo else "-force"
        mask_name = "mask"
        results = {
            'subject_id': subject_id,
            'session_id': session_id,
            'input': input_nii,
            'output_dir': output_dir,
            'steps_run': [],
            'final_file': None,
            'success': True,
        }
        
        log.info("=" * 60)
        log.info("Processing: sub-%s ses-%s", subject_id, session_id)
        log.info("Input: %s", input_nii)
        log.info("Output dir: %s", output_dir)
        log.info("Steps: %s", {k: v for k, v in step_config.items() if v})
        log.info("=" * 60)
        
        try:
            # Step 1: Convert NIfTI to MIF
            if step_config.get('convert', True):
                difm = do_convert(input_nii, bvec, bval, output_dir, 
                                  common=common, dry_run=dry_run)
                results['steps_run'].append('convert')
            else:
                # Assume input is already MIF
                difm = Path(input_nii).stem
            
            # Step 2: Axis flip for orientation correction
            if step_config.get('flip', True):
                difm = do_flip(difm, output_dir, common=common, dry_run=dry_run)
                results['steps_run'].append('flip')
            
            # Step 3: Gradient check
            if step_config.get('gradcheck', True):
                difm = do_gradcheck(difm, output_dir, common=common, 
                                    do_step=True, dry_run=dry_run)
                results['steps_run'].append('gradcheck')
            
            # Step 4: Denoising
            if step_config.get('denoise', True):
                difm = do_denoise(difm, output_dir, common=common,
                                  do_step=True, dry_run=dry_run)
                results['steps_run'].append('denoise')
            
            # Step 5: Gibbs ringing removal
            if step_config.get('degibbs', False):
                difm = do_degibbs(difm, output_dir, common=common,
                                  do_step=True, dry_run=dry_run)
                results['steps_run'].append('degibbs')
            
            # Step 6: Create mask
            if step_config.get('mask', True):
                mask_name = do_mask(difm, output_dir, mask_name=mask_name,
                                    common=common, do_step=True, dry_run=dry_run)
                results['steps_run'].append('mask')
            
            # Step 7: Bias correction
            if step_config.get('bias', True):
                difm = do_bias(difm, mask_name, output_dir, 
                               bias_method="ants", common=common,
                               do_step=True, dry_run=dry_run)
                results['steps_run'].append('bias')
            
            # Step 8: Intensity normalization (optional)
            if step_config.get('norm', False):
                difm = do_norm(difm, mask_name, output_dir,
                               common=common, do_step=True, dry_run=dry_run)
                results['steps_run'].append('norm')
            
            # Create reference images
            create_ref_images(difm, mask_name, output_dir, 
                              common=common, dry_run=dry_run)
            
            # =================================================================
            # T2w ANATOMICAL PROCESSING
            # =================================================================
            t2w_lowres_name = None
            t2w_highres_name = None
            
            # Step 9a: Convert and flip lowres T2w
            if step_config.get('convert_t2w', True) and t2w_lowres:
                log.info("Processing lowres T2w anatomical...")
                t2w_lowres_name = do_convert_t2w(t2w_lowres, output_dir, 
                                                  acq_type="lowres",
                                                  common=common, dry_run=dry_run)
                results['steps_run'].append('convert_t2w_lowres')
                
                if step_config.get('flip_t2w', True):
                    t2w_lowres_name = do_flip_t2w(t2w_lowres_name, output_dir,
                                                   acq_type="lowres",
                                                   common=common, dry_run=dry_run)
                    results['steps_run'].append('flip_t2w_lowres')
            
            # Step 9b: Convert and flip highres T2w
            if step_config.get('convert_t2w', True) and t2w_highres:
                log.info("Processing highres T2w anatomical...")
                t2w_highres_name = do_convert_t2w(t2w_highres, output_dir,
                                                   acq_type="highres",
                                                   common=common, dry_run=dry_run)
                results['steps_run'].append('convert_t2w_highres')
                
                if step_config.get('flip_t2w', True):
                    t2w_highres_name = do_flip_t2w(t2w_highres_name, output_dir,
                                                    acq_type="highres",
                                                    common=common, dry_run=dry_run)
                    results['steps_run'].append('flip_t2w_highres')
            
            # Choose which T2w to use for atlas and DWI coregistration
            # Prefer highres if available, fallback to lowres
            t2w_for_coreg = t2w_highres_name or t2w_lowres_name
            
            # =================================================================
            # ATLAS COREGISTRATION
            # =================================================================
            # Step 10: Atlas to T2w coregistration
            if step_config.get('atlas_coreg', True) and t2w_for_coreg and atlas_template:
                log.info("Performing atlas to T2w coregistration...")
                atlas_results = do_atlas_coreg(
                    t2wm=t2w_for_coreg,
                    output_dir=output_dir,
                    atlas_template=atlas_template,
                    atlas_annotation=atlas_annotation or "",
                    common=common,
                    dry_run=dry_run
                )
                results['steps_run'].append('atlas_coreg')
                results['atlas_results'] = atlas_results
                
                # Step 11: Check ROIs
                if step_config.get('check_rois', True) and atlas_annotation:
                    log.info("Checking ROI placement...")
                    roi_results = do_check_rois(
                        atlas_annotation_warped=atlas_results['warped_annotation'],
                        output_dir=output_dir,
                        common=common,
                        dry_run=dry_run
                    )
                    results['steps_run'].append('check_rois')
                    results['roi_masks'] = roi_results
            elif step_config.get('atlas_coreg', True) and not t2w_for_coreg:
                log.warning("Atlas coregistration requested but no T2w available")
            elif step_config.get('atlas_coreg', True) and not atlas_template:
                log.warning("Atlas coregistration requested but no atlas template provided")
            
            # =================================================================
            # DWI TO T2w COREGISTRATION (FINAL STEP)
            # =================================================================
            # Step 12: DWI to T2w coregistration - THIS IS THE LAST STEP
            if step_config.get('dwi2anat', True) and t2w_for_coreg:
                log.info("Performing DWI to T2w coregistration (final step)...")
                difm = do_dwi2anat(
                    difm=difm,
                    mask=mask_name,
                    t2wm=t2w_for_coreg,
                    output_dir=output_dir,
                    common=common,
                    do_step=True,
                    dry_run=dry_run
                )
                results['steps_run'].append('dwi2anat')
            elif step_config.get('dwi2anat', True) and not t2w_for_coreg:
                log.warning("DWI to T2w coregistration requested but no T2w available")
            
            # Export to NIfTI (optional)
            if step_config.get('export_nifti', False):
                final_nii = do_export_nifti(difm, output_dir, 
                                            common=common, dry_run=dry_run)
                results['final_file'] = final_nii
                results['steps_run'].append('export_nifti')
            else:
                results['final_file'] = os.path.join(output_dir, f"{difm}.mif")
            
            log.info("✓ Preprocessing complete for sub-%s ses-%s", 
                     subject_id, session_id)
            log.info("Final file: %s", results['final_file'])
            
        except Exception as e:
            log.error("❌ Error processing sub-%s ses-%s: %s", 
                      subject_id, session_id, e)
            results['success'] = False
            results['error'] = str(e)
            self.stats['errors'].append(f"{subject_id}/{session_id}: {e}")
        
        return results
    
    def preprocess_all_subjects(self,
                                derivatives_dir: Optional[Path] = None,
                                steps: Optional[Dict[str, bool]] = None,
                                subjects: Optional[List[str]] = None,
                                first_subject_only: bool = False,
                                dry_run: bool = False,
                                force_redo: bool = False,
                                input_dwi_name: str = "_dwi.nii.gz") -> Dict[str, Any]:
        """
        Preprocess all subjects in the BIDS dataset.
        
        Parameters
        ----------
        derivatives_dir : Path, optional
            Override derivatives directory
        steps : dict, optional
            Step configuration
        subjects : list, optional
            List of specific subject IDs to process
        first_subject_only : bool
            If True, only process first subject (for testing)
        dry_run : bool
            If True, show commands without executing
        force_redo : bool
            If True, redo all steps
        input_dwi_name : str
            Pattern to match DWI files
            
        Returns
        -------
        dict
            Processing statistics
        """
        if derivatives_dir:
            self.preproc_dir = Path(derivatives_dir)
        
        # Find all subjects
        bids_dir = self.bids_dir
        subject_dirs = sorted([d for d in bids_dir.iterdir() 
                               if d.is_dir() and d.name.startswith('sub-')])
        
        if subjects:
            subject_dirs = [d for d in subject_dirs 
                           if d.name.replace('sub-', '') in subjects]
        
        if first_subject_only and subject_dirs:
            subject_dirs = subject_dirs[:1]
        
        log.info("Found %d subjects to process", len(subject_dirs))
        
        for subject_dir in subject_dirs:
            subject_id = subject_dir.name.replace('sub-', '')
            
            # Find sessions
            session_dirs = sorted([d for d in subject_dir.iterdir()
                                   if d.is_dir() and d.name.startswith('ses-')])
            
            for session_dir in session_dirs:
                session_id = session_dir.name.replace('ses-', '')
                
                # Find DWI files
                dwi_dir = session_dir / "dwi"
                if not dwi_dir.exists():
                    log.warning("No DWI directory for sub-%s ses-%s", 
                                subject_id, session_id)
                    continue
                
                dwi_files = list(dwi_dir.glob(f"*{input_dwi_name}"))
                if not dwi_files:
                    log.warning("No DWI files matching %s in %s", 
                                input_dwi_name, dwi_dir)
                    continue
                
                # Find T2w anatomical files
                anat_dir = session_dir / "anat"
                t2w_lowres = None
                t2w_highres = None
                
                if anat_dir.exists():
                    lowres_files = list(anat_dir.glob("*_acq-lowres_T2w.nii.gz"))
                    highres_files = list(anat_dir.glob("*_acq-highres_T2w.nii.gz"))
                    
                    if lowres_files:
                        t2w_lowres = str(lowres_files[0])
                        log.info("Found lowres T2w: %s", t2w_lowres)
                    if highres_files:
                        t2w_highres = str(highres_files[0])
                        log.info("Found highres T2w: %s", t2w_highres)
                else:
                    log.warning("No anat directory for sub-%s ses-%s", 
                                subject_id, session_id)
                
                # Get atlas files if available
                atlas_template = None
                atlas_annotation = None
                try:
                    from .atlas_manager import AllenAtlasManager
                    atlas_manager = AllenAtlasManager()
                    if atlas_manager.is_downloaded():
                        atlas_files = atlas_manager.get_atlas_files()
                        atlas_template = str(atlas_files['template'])
                        atlas_annotation = str(atlas_files['annotation'])
                        log.info("Found atlas files")
                except ImportError:
                    log.warning("Atlas manager not available")
                
                for dwi_file in dwi_files:
                    # Find bvec/bval
                    base = str(dwi_file).replace('.nii.gz', '').replace('.nii', '')
                    bvec = f"{base}.bvec"
                    bval = f"{base}.bval"
                    
                    if not os.path.isfile(bvec) or not os.path.isfile(bval):
                        log.warning("Missing bvec/bval for %s", dwi_file)
                        continue
                    
                    # Output directory
                    output_dir = str(self.preproc_dir / f"sub-{subject_id}" / 
                                     f"ses-{session_id}" / "dwi")
                    
                    # Process
                    result = self.preprocess_dwi(
                        input_nii=str(dwi_file),
                        bvec=bvec,
                        bval=bval,
                        output_dir=output_dir,
                        subject_id=subject_id,
                        session_id=session_id,
                        t2w_lowres=t2w_lowres,
                        t2w_highres=t2w_highres,
                        atlas_template=atlas_template,
                        atlas_annotation=atlas_annotation,
                        steps=steps,
                        dry_run=dry_run,
                        force_redo=force_redo
                    )
                    
                    if result['success']:
                        self.stats['subjects_processed'] += 1
                        self.stats['steps_completed'].extend(result['steps_run'])
        
        return self.stats
