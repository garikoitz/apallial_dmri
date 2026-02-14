#!/usr/bin/env python3
"""
Mouse Brain DWI Processing Pipeline

Complete pipeline for preprocessing and tractography of mouse brain DWI data.
Processes multiple subjects with intermediate quality checks.

Usage:
    python run_pipeline.py --subjects sub-596-C-OO-NmM sub-769-C-OO-NmM
    python run_pipeline.py --all
    python run_pipeline.py --subjects sub-596-C-OO-NmM --steps preproc
    python run_pipeline.py --subjects sub-596-C-OO-NmM --steps tractography --skip-qc
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger(__name__)


class PipelineConfig:
    """Configuration for the pipeline."""
    
    def __init__(self, base_dir: str = None):
        self.base_dir = Path(base_dir or "/Users/glerma/soft/mbrain")
        self.bids_dir = self.base_dir / "data" / "bids"
        self.deriv_dir = self.base_dir / "data" / "derivatives" / "mrtrix_preproc"
        self.atlas_dir = self.base_dir / "local" / "atlas" / "allen_ccf"
        
        # Atlas files
        self.atlas_template = self.atlas_dir / "average_template_25.nii.gz"
        self.atlas_annotation = self.atlas_dir / "annotation_25.nii.gz"
        self.atlas_rois_dir = self.atlas_dir / "rois"
        
        # ROI names to process
        self.roi_names = ["isocortex", "striatum", "thalamus", "hippocampus", "cerebellum"]
        
        # MRtrix config
        self.bzero_threshold = 20
        
        # Processing parameters
        self.denoise_extent = "5,5,5"
        self.degibbs_nshifts = 20
        self.bias_bspline = "[200,3]"
        self.bias_convergence = "[10,1e-6]"
        self.bias_shrink = 4
        
        # Tractography parameters
        self.tck_select = 2000
        self.tck_cutoff = 0.1
        self.tck_minlength = 5
        self.tck_maxlength = 25
        self.tck_step = 0.05
        self.tck_angle = 45
        self.sift_term_number = 50
        
        # FA threshold for WM mask
        self.fa_threshold = 0.2


class PipelineRunner:
    """Main pipeline runner."""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.common = "-force -quiet"
        
    def run_cmd(self, cmd: str, description: str = None, check: bool = True) -> subprocess.CompletedProcess:
        """Run a shell command."""
        if description:
            log.info(f"  {description}")
        log.debug(f"  CMD: {cmd}")
        
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True
        )
        
        if check and result.returncode != 0:
            log.error(f"Command failed: {cmd}")
            log.error(f"STDERR: {result.stderr}")
            raise RuntimeError(f"Command failed with exit code {result.returncode}")
        
        return result
    
    def file_exists(self, path: Path) -> bool:
        """Check if file exists."""
        return path.exists()
    
    def get_subjects(self) -> List[str]:
        """Get list of all subjects in BIDS directory."""
        subjects = []
        for item in self.config.bids_dir.iterdir():
            if item.is_dir() and item.name.startswith("sub-"):
                subjects.append(item.name)
        return sorted(subjects)
    
    def get_sessions(self, subject: str) -> List[str]:
        """Get list of sessions for a subject."""
        subj_dir = self.config.bids_dir / subject
        sessions = []
        for item in subj_dir.iterdir():
            if item.is_dir() and item.name.startswith("ses-"):
                sessions.append(item.name)
        return sorted(sessions)
    
    def setup_output_dir(self, subject: str, session: str) -> Path:
        """Create output directory structure."""
        out_dir = self.config.deriv_dir / subject / session / "dwi"
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "atlas" / "rois").mkdir(parents=True, exist_ok=True)
        (out_dir / "qc").mkdir(parents=True, exist_ok=True)
        return out_dir

    # =========================================================================
    # PART 1: T2w ANATOMICAL PROCESSING
    # =========================================================================
    
    def process_t2w(self, subject: str, session: str, out_dir: Path) -> Dict[str, Path]:
        """Process T2w anatomical images."""
        log.info("=== Processing T2w anatomical images ===")
        
        bids_anat = self.config.bids_dir / subject / session / "anat"
        outputs = {}
        
        # Find T2w files (highres and lowres)
        for suffix in ["highres", "lowres"]:
            t2w_nii = bids_anat / f"{subject}_{session}_acq-{suffix}_T2w.nii.gz"
            
            if not t2w_nii.exists():
                # Try without acq- label
                t2w_nii = bids_anat / f"{subject}_{session}_T2w.nii.gz"
                if not t2w_nii.exists():
                    log.warning(f"  T2w {suffix} not found, skipping")
                    continue
            
            # Convert to MIF
            t2w_mif = out_dir / f"t2w_{suffix}.mif"
            self.run_cmd(
                f"mrconvert {t2w_nii} {t2w_mif} {self.common}",
                f"Converting T2w {suffix} to MIF"
            )
            
            # Step 1: Axis permutation (for 3D: 2,1,0)
            t2w_perm = out_dir / f"t2w_{suffix}_perm.mif"
            self.run_cmd(
                f"mrconvert {t2w_mif} {t2w_perm} -axes 2,1,0 -strides 1,2,3 {self.common}",
                f"Permuting T2w {suffix} axes"
            )
            
            # Step 2: Axis flip
            t2w_flip = out_dir / f"t2w_{suffix}_flip.mif"
            self.run_cmd(
                f"mrtransform {t2w_perm} {t2w_flip} -flip 1 {self.common}",
                f"Flipping T2w {suffix}"
            )
            
            # Clean up intermediate file
            t2w_perm.unlink(missing_ok=True)
            
            # Also save as NIfTI for ANTs
            t2w_flip_nii = out_dir / f"t2w_{suffix}_flip.nii.gz"
            self.run_cmd(
                f"mrconvert {t2w_flip} {t2w_flip_nii} {self.common}",
                f"Converting flipped T2w {suffix} to NIfTI"
            )
            
            outputs[suffix] = t2w_flip_nii
        
        return outputs

    # =========================================================================
    # PART 2: ATLAS TO T2w COREGISTRATION
    # =========================================================================
    
    def register_atlas_to_t2w(self, out_dir: Path, t2w_file: Path) -> Dict[str, Path]:
        """Register Allen CCF atlas to subject T2w space."""
        log.info("=== Registering atlas to T2w ===")
        
        atlas_dir = out_dir / "atlas"
        outputs = {}
        
        # Run ANTs registration: atlas template -> T2w
        prefix = atlas_dir / "atlas2t2w_"
        self.run_cmd(
            f"antsRegistrationSyN.sh -d 3 -t s -f {t2w_file} "
            f"-m {self.config.atlas_template} -o {prefix}",
            "Running ANTs SyN registration (atlas -> T2w)"
        )
        
        outputs["warp"] = atlas_dir / "atlas2t2w_1Warp.nii.gz"
        outputs["affine"] = atlas_dir / "atlas2t2w_0GenericAffine.mat"
        
        # Warp annotation to T2w space
        annotation_warped = atlas_dir / "annotation_to_t2w.nii.gz"
        self.run_cmd(
            f"antsApplyTransforms -d 3 -i {self.config.atlas_annotation} "
            f"-r {t2w_file} -o {annotation_warped} "
            f"-t {outputs['warp']} -t {outputs['affine']} "
            f"-n NearestNeighbor",
            "Warping annotation to T2w space"
        )
        outputs["annotation"] = annotation_warped
        
        # Warp template to T2w space (for QC)
        template_warped = atlas_dir / "template_to_t2w.nii.gz"
        self.run_cmd(
            f"antsApplyTransforms -d 3 -i {self.config.atlas_template} "
            f"-r {t2w_file} -o {template_warped} "
            f"-t {outputs['warp']} -t {outputs['affine']}",
            "Warping template to T2w space"
        )
        outputs["template"] = template_warped
        
        return outputs
    
    def warp_rois_to_subject(self, out_dir: Path, t2w_file: Path, 
                             warp: Path, affine: Path) -> Dict[str, Path]:
        """Warp atlas ROIs to subject T2w space."""
        log.info("=== Warping ROIs to subject space ===")
        
        rois_dir = out_dir / "atlas" / "rois"
        outputs = {}
        
        for roi_name in self.config.roi_names:
            roi_src = self.config.atlas_rois_dir / f"roi_{roi_name}.nii.gz"
            
            if not roi_src.exists():
                log.warning(f"  ROI {roi_name} not found in atlas, skipping")
                continue
            
            roi_dst = rois_dir / f"roi_{roi_name}.nii.gz"
            self.run_cmd(
                f"antsApplyTransforms -d 3 -i {roi_src} -r {t2w_file} "
                f"-o {roi_dst} -t {warp} -t {affine} -n NearestNeighbor",
                f"Warping ROI: {roi_name}"
            )
            outputs[roi_name] = roi_dst
        
        return outputs

    # =========================================================================
    # PART 3: DWI PREPROCESSING
    # =========================================================================
    
    def preprocess_dwi(self, subject: str, session: str, out_dir: Path) -> Dict[str, Path]:
        """Complete DWI preprocessing pipeline."""
        log.info("=== DWI Preprocessing ===")
        
        bids_dwi = self.config.bids_dir / subject / session / "dwi"
        outputs = {}
        
        # Find DWI files
        dwi_nii = bids_dwi / f"{subject}_{session}_dwi.nii.gz"
        dwi_bvec = bids_dwi / f"{subject}_{session}_dwi.bvec"
        dwi_bval = bids_dwi / f"{subject}_{session}_dwi.bval"
        
        if not all(f.exists() for f in [dwi_nii, dwi_bvec, dwi_bval]):
            raise FileNotFoundError(f"DWI files not found for {subject}/{session}")
        
        # Track current filename
        difm = out_dir / "dwi"
        
        # 1. Convert to MIF with gradients
        self.run_cmd(
            f"mrconvert {dwi_nii} {difm}.mif -fslgrad {dwi_bvec} {dwi_bval} {self.common}",
            "Converting DWI to MIF"
        )
        
        # 2. Flip orientation (two-step process for 4D data)
        # Step 2a: Axis permutation (for 4D: 2,1,0,3)
        difm_perm = out_dir / "dwi_perm"
        self.run_cmd(
            f"mrconvert {difm}.mif {difm_perm}.mif -axes 2,1,0,3 -strides 1,2,3,4 {self.common}",
            "Permuting DWI axes"
        )
        
        # Step 2b: Axis flip
        difm_flip = out_dir / "dwi_flip"
        self.run_cmd(
            f"mrtransform {difm_perm}.mif {difm_flip}.mif -flip 1 -reorient_fod no {self.common}",
            "Flipping DWI orientation"
        )
        
        # Clean up intermediate file
        (difm_perm.with_suffix(".mif")).unlink(missing_ok=True)
        difm = difm_flip
        
        # 3. Gradient check
        difm_gradcheck = out_dir / "dwi_flip_gradcheck"
        corrected_grad = out_dir / "corrected.b"
        self.run_cmd(
            f"dwigradcheck {difm}.mif -export_grad_mrtrix {corrected_grad} {self.common}",
            "Checking gradients"
        )
        self.run_cmd(
            f"mrconvert {difm}.mif {difm_gradcheck}.mif -grad {corrected_grad} {self.common}",
            "Applying corrected gradients"
        )
        difm = difm_gradcheck
        
        # 4. Denoise
        difm_denoise = out_dir / "dwi_flip_gradcheck_denoise"
        self.run_cmd(
            f"dwidenoise -extent {self.config.denoise_extent} {difm}.mif {difm_denoise}.mif {self.common}",
            "Denoising"
        )
        difm = difm_denoise
        
        # 5. Create mask
        mask = out_dir / "mask"
        self.run_cmd(
            f"dwi2mask {difm}.mif {mask}.mif {self.common}",
            "Creating brain mask"
        )
        self.run_cmd(
            f"mrconvert {mask}.mif {mask}.nii.gz {self.common}",
            "Converting mask to NIfTI"
        )
        outputs["mask"] = mask
        
        # 6. Extract b0 (need for bias correction)
        b0 = out_dir / "b0_dwi"
        self.run_cmd(
            f"dwiextract {difm}.mif - -bzero -config BZeroThreshold {self.config.bzero_threshold} {self.common} | "
            f"mrmath - mean {b0}.mif -axis 3 {self.common}",
            "Extracting mean b0"
        )
        self.run_cmd(
            f"mrconvert {b0}.mif {b0}.nii.gz {self.common}",
            "Converting b0 to NIfTI"
        )
        
        # 7. Bias correction
        difm_bias = out_dir / "dwi_flip_gradcheck_denoise_bias"
        scratch_dir = out_dir / "tmp_bias"
        scratch_dir.mkdir(exist_ok=True)
        self.run_cmd(
            f"dwibiascorrect ants -mask {mask}.mif "
            f"-ants.b \"{self.config.bias_bspline}\" "
            f"-ants.c \"{self.config.bias_convergence}\" "
            f"-ants.s {self.config.bias_shrink} "
            f"{difm}.mif {difm_bias}.mif "
            f"-scratch {scratch_dir} {self.common}",
            "Bias field correction (ANTs N4)"
        )
        difm = difm_bias
        
        # Clean up scratch
        shutil.rmtree(scratch_dir, ignore_errors=True)
        
        outputs["dwi"] = difm
        outputs["b0"] = b0
        
        return outputs
    
    def register_dwi_to_t2w(self, out_dir: Path, dwi_file: Path, 
                            mask_file: Path, t2w_file: Path) -> Dict[str, Path]:
        """Register DWI to T2w anatomical space."""
        log.info("=== Registering DWI to T2w ===")
        
        outputs = {}
        
        # Extract b0 from bias-corrected DWI
        b0_bias = out_dir / "b0_dwi_bias"
        self.run_cmd(
            f"dwiextract {dwi_file}.mif - -bzero -config BZeroThreshold {self.config.bzero_threshold} {self.common} | "
            f"mrmath - mean {b0_bias}.mif -axis 3 {self.common}",
            "Extracting b0 from bias-corrected DWI"
        )
        self.run_cmd(
            f"mrconvert {b0_bias}.mif {b0_bias}.nii.gz {self.common}",
            "Converting b0 to NIfTI"
        )
        
        # Apply mask to b0
        b0_brain = out_dir / "b0_dwi_bias_brain"
        self.run_cmd(
            f"mrcalc {b0_bias}.nii.gz {mask_file}.nii.gz -mult {b0_brain}.nii.gz {self.common}",
            "Masking b0"
        )
        
        # ANTs registration: b0 -> T2w
        prefix = out_dir / "dwi2anat_"
        self.run_cmd(
            f"antsRegistrationSyN.sh -d 3 -t r -f {t2w_file} -m {b0_brain}.nii.gz -o {prefix}",
            "Running ANTs rigid registration (DWI -> T2w)"
        )
        
        outputs["affine_mat"] = out_dir / "dwi2anat_0GenericAffine.mat"
        
        # Convert transform for MRtrix
        affine_txt = out_dir / "dwi2anat_affine.txt"
        mrtrix_txt = out_dir / "dwi2anat_mrtrix.txt"
        
        self.run_cmd(
            f"ConvertTransformFile 3 {outputs['affine_mat']} {affine_txt}",
            "Converting ANTs transform"
        )
        self.run_cmd(
            f"transformconvert {affine_txt} itk_import {mrtrix_txt} {self.common}",
            "Converting to MRtrix format"
        )
        outputs["mrtrix_transform"] = mrtrix_txt
        
        # Apply transform to DWI
        dwi_anat = Path(str(dwi_file) + "_anat")
        self.run_cmd(
            f"mrtransform {dwi_file}.mif {dwi_anat}.mif "
            f"-linear {mrtrix_txt} -reorient_fod no {self.common}",
            "Applying transform to DWI"
        )
        outputs["dwi_anat"] = dwi_anat
        
        return outputs

    # =========================================================================
    # PART 4: TRACTOGRAPHY
    # =========================================================================
    
    def run_tractography(self, out_dir: Path, dwi_file: Path, mask_file: Path) -> Dict[str, Path]:
        """Run tractography pipeline."""
        log.info("=== Running Tractography ===")
        
        outputs = {}
        
        # 1. Fit tensor
        tensor = out_dir / "tensor"
        self.run_cmd(
            f"dwi2tensor -mask {mask_file}.mif {dwi_file}.mif {tensor}.mif {self.common}",
            "Fitting diffusion tensor"
        )
        outputs["tensor"] = tensor
        
        # 2. Extract FA
        fa = out_dir / "fa"
        self.run_cmd(
            f"tensor2metric {tensor}.mif -fa {fa}.mif -mask {mask_file}.mif {self.common}",
            "Extracting FA"
        )
        outputs["fa"] = fa
        
        # 3. Extract MD
        md = out_dir / "md"
        self.run_cmd(
            f"tensor2metric {tensor}.mif -adc {md}.mif -mask {mask_file}.mif {self.common}",
            "Extracting MD"
        )
        outputs["md"] = md
        
        # 4. Extract AD (Axial Diffusivity)
        ad = out_dir / "ad"
        self.run_cmd(
            f"tensor2metric {tensor}.mif -ad {ad}.mif -mask {mask_file}.mif {self.common}",
            "Extracting AD"
        )
        outputs["ad"] = ad
        
        # 5. Extract RD (Radial Diffusivity)
        rd = out_dir / "rd"
        self.run_cmd(
            f"tensor2metric {tensor}.mif -rd {rd}.mif -mask {mask_file}.mif {self.common}",
            "Extracting RD"
        )
        outputs["rd"] = rd
        
        # 6. Create WM mask from FA
        wm_mask = out_dir / "wm_mask"
        self.run_cmd(
            f"mrthreshold {fa}.mif -abs {self.config.fa_threshold} {wm_mask}.mif {self.common}",
            f"Creating WM mask (FA > {self.config.fa_threshold})"
        )
        outputs["wm_mask"] = wm_mask
        
        # 5. Estimate response function
        response = out_dir / "response_wm.txt"
        response_voxels = out_dir / "response_voxels"
        self.run_cmd(
            f"dwi2response tournier {dwi_file}.mif {response} "
            f"-mask {mask_file}.mif -voxels {response_voxels}.mif {self.common}",
            "Estimating response function"
        )
        outputs["response"] = response
        
        # 6. Compute FODs
        fod = out_dir / "fod_wm"
        self.run_cmd(
            f"dwi2fod csd {dwi_file}.mif {response} {fod}.mif "
            f"-mask {mask_file}.mif {self.common}",
            "Computing FODs (CSD)"
        )
        
        # 7. Normalize FODs
        fod_norm = out_dir / "fod_wm_norm"
        self.run_cmd(
            f"mtnormalise {fod}.mif {fod_norm}.mif -mask {mask_file}.mif {self.common}",
            "Normalizing FODs"
        )
        outputs["fod"] = fod_norm
        
        # 8. Whole-brain tractography
        tracks_wb = out_dir / "tracks_wholebrain"
        self.run_cmd(
            f"tckgen {fod_norm}.mif {tracks_wb}.tck "
            f"-algorithm iFOD2 -seed_image {wm_mask}.mif -mask {mask_file}.mif "
            f"-select 100000 -cutoff {self.config.tck_cutoff} "
            f"-minlength {self.config.tck_minlength} -maxlength {self.config.tck_maxlength} "
            f"-step {self.config.tck_step} -angle {self.config.tck_angle} {self.common}",
            "Generating whole-brain tractogram"
        )
        outputs["tracks_wholebrain"] = tracks_wb
        
        return outputs
    
    def run_roi_tractography(self, out_dir: Path, fod_file: Path, mask_file: Path,
                              seed_roi: Path, include_roi: Path, 
                              tract_name: str) -> Dict[str, Path]:
        """Run ROI-to-ROI tractography with SIFT."""
        log.info(f"=== ROI Tractography: {tract_name} ===")
        
        outputs = {}
        
        # Generate raw tracks
        tracks_raw = out_dir / f"tracks_{tract_name}_raw"
        self.run_cmd(
            f"tckgen {fod_file}.mif {tracks_raw}.tck "
            f"-algorithm iFOD2 -seed_image {seed_roi} -include {include_roi} "
            f"-mask {mask_file}.mif -select {self.config.tck_select} "
            f"-cutoff {self.config.tck_cutoff} "
            f"-minlength {self.config.tck_minlength} -maxlength {self.config.tck_maxlength} "
            f"-step {self.config.tck_step} -angle {self.config.tck_angle} {self.common}",
            f"Generating {tract_name} tracks"
        )
        outputs["raw"] = tracks_raw
        
        # Apply SIFT
        tracks_sift = out_dir / f"tracks_{tract_name}_sift{self.config.sift_term_number}"
        self.run_cmd(
            f"tcksift {tracks_raw}.tck {fod_file}.mif {tracks_sift}.tck "
            f"-term_number {self.config.sift_term_number} {self.common}",
            f"Applying SIFT (keeping {self.config.sift_term_number} tracks)"
        )
        outputs["sift"] = tracks_sift
        
        return outputs
    
    def create_tract_density(self, tracks_file: Path, template: Path, 
                              output: Path) -> Path:
        """Create tract density map in template space."""
        self.run_cmd(
            f"tckmap {tracks_file}.tck {output}.mif -template {template} {self.common}",
            "Creating tract density map"
        )
        self.run_cmd(
            f"mrconvert {output}.mif {output}.nii.gz {self.common}",
            "Converting to NIfTI"
        )
        return output

    # =========================================================================
    # QUALITY CONTROL
    # =========================================================================
    
    def generate_qc_report(self, subject: str, session: str, out_dir: Path,
                           outputs: Dict[str, Any]) -> Path:
        """Generate QC report for a subject."""
        log.info("=== Generating QC Report ===")
        
        qc_dir = out_dir / "qc"
        report_file = qc_dir / f"qc_report_{subject}_{session}.md"
        
        report = []
        report.append(f"# QC Report: {subject} / {session}")
        report.append(f"\nGenerated: {datetime.now().isoformat()}")
        report.append("\n## File Inventory\n")
        
        # Check key files
        key_files = [
            ("T2w highres", out_dir / "t2w_highres_flip.nii.gz"),
            ("DWI preprocessed", out_dir / "dwi_flip_gradcheck_denoise_bias_anat.mif"),
            ("Brain mask", out_dir / "mask.mif"),
            ("FA map", out_dir / "fa.mif"),
            ("FOD", out_dir / "fod_wm_norm.mif"),
            ("Atlas annotation", out_dir / "atlas" / "annotation_to_t2w.nii.gz"),
            ("Thalamus ROI", out_dir / "atlas" / "rois" / "roi_thalamus.nii.gz"),
            ("Thalamus-cortex tracks", out_dir / "tracks_thalamus_cortex_sift50.tck"),
        ]
        
        for name, path in key_files:
            status = "✓" if path.exists() else "✗"
            size = f"({path.stat().st_size / 1024 / 1024:.1f} MB)" if path.exists() else "(missing)"
            report.append(f"- [{status}] {name}: `{path.name}` {size}")
        
        # Get FA stats if available
        fa_file = out_dir / "fa.mif"
        if fa_file.exists():
            result = self.run_cmd(
                f"mrstats {fa_file} -mask {out_dir / 'mask.mif'} -output mean,std,min,max",
                check=False
            )
            if result.returncode == 0:
                report.append("\n## FA Statistics\n")
                report.append(f"```\n{result.stdout}\n```")
        
        # Track info
        tracks_file = out_dir / "tracks_thalamus_cortex_sift50.tck"
        if tracks_file.exists():
            result = self.run_cmd(f"tckinfo {tracks_file}", check=False)
            if result.returncode == 0:
                report.append("\n## Thalamus-Cortex Track Info\n")
                report.append(f"```\n{result.stdout}\n```")
        
        report.append("\n## Visual QC Commands\n")
        report.append("```bash")
        report.append(f"# View FA with ROIs")
        report.append(f"mrview {out_dir}/fa.mif \\")
        report.append(f"    -overlay.load {out_dir}/atlas/rois/roi_thalamus.nii.gz \\")
        report.append(f"    -overlay.load {out_dir}/atlas/rois/roi_isocortex.nii.gz &")
        report.append("")
        report.append(f"# View tracks")
        report.append(f"mrview {out_dir}/fa.mif \\")
        report.append(f"    -tractography.load {out_dir}/tracks_thalamus_cortex_sift50.tck &")
        report.append("```")
        
        # Write report
        with open(report_file, 'w') as f:
            f.write('\n'.join(report))
        
        log.info(f"  QC report saved to: {report_file}")
        return report_file

    # =========================================================================
    # MAIN PIPELINE
    # =========================================================================
    
    def process_subject(self, subject: str, session: str = None,
                        steps: List[str] = None, skip_qc: bool = False) -> Dict[str, Any]:
        """Process a single subject through the complete pipeline."""
        
        log.info(f"\n{'='*60}")
        log.info(f"Processing: {subject}")
        log.info(f"{'='*60}")
        
        # Get session if not specified
        if session is None:
            sessions = self.get_sessions(subject)
            if not sessions:
                raise ValueError(f"No sessions found for {subject}")
            session = sessions[0]
            log.info(f"Using session: {session}")
        
        # Setup output directory
        out_dir = self.setup_output_dir(subject, session)
        log.info(f"Output directory: {out_dir}")
        
        # Default steps
        if steps is None:
            steps = ["preproc", "tractography"]
        
        outputs = {"subject": subject, "session": session, "out_dir": out_dir}
        
        try:
            # === PREPROCESSING ===
            if "preproc" in steps:
                log.info("\n>>> PREPROCESSING PHASE <<<")
                
                # 1. Process T2w
                t2w_outputs = self.process_t2w(subject, session, out_dir)
                outputs["t2w"] = t2w_outputs
                
                t2w_file = t2w_outputs.get("highres") or t2w_outputs.get("lowres")
                if not t2w_file:
                    raise FileNotFoundError("No T2w image found")
                
                # 2. Register atlas to T2w
                atlas_outputs = self.register_atlas_to_t2w(out_dir, t2w_file)
                outputs["atlas"] = atlas_outputs
                
                # 3. Warp ROIs to subject
                roi_outputs = self.warp_rois_to_subject(
                    out_dir, t2w_file,
                    atlas_outputs["warp"], atlas_outputs["affine"]
                )
                outputs["rois"] = roi_outputs
                
                # 4. Preprocess DWI
                dwi_outputs = self.preprocess_dwi(subject, session, out_dir)
                outputs["dwi_preproc"] = dwi_outputs
                
                # 5. Register DWI to T2w
                dwi_reg_outputs = self.register_dwi_to_t2w(
                    out_dir, dwi_outputs["dwi"],
                    dwi_outputs["mask"], t2w_file
                )
                outputs["dwi_reg"] = dwi_reg_outputs
            
            # === TRACTOGRAPHY ===
            if "tractography" in steps:
                log.info("\n>>> TRACTOGRAPHY PHASE <<<")
                
                # Get DWI file (either from preproc or existing)
                dwi_anat = out_dir / "dwi_flip_gradcheck_denoise_bias_anat"
                mask = out_dir / "mask"
                
                if not (dwi_anat.with_suffix(".mif")).exists():
                    raise FileNotFoundError(
                        "Preprocessed DWI not found. Run preproc step first."
                    )
                
                # 1. Run tractography
                tract_outputs = self.run_tractography(out_dir, dwi_anat, mask)
                outputs["tractography"] = tract_outputs
                
                # 2. ROI tractography (thalamus -> cortex)
                thalamus_roi = out_dir / "atlas" / "rois" / "roi_thalamus.nii.gz"
                cortex_roi = out_dir / "atlas" / "rois" / "roi_isocortex.nii.gz"
                
                if thalamus_roi.exists() and cortex_roi.exists():
                    roi_tract_outputs = self.run_roi_tractography(
                        out_dir, tract_outputs["fod"], mask,
                        thalamus_roi, cortex_roi, "thalamus_cortex"
                    )
                    outputs["roi_tractography"] = roi_tract_outputs
                    
                    # 3. Create tract density map
                    density = out_dir / "thalamus_cortex_density"
                    self.create_tract_density(
                        roi_tract_outputs["sift"],
                        self.config.atlas_template,
                        density
                    )
                    outputs["tract_density"] = density
            
            # === QC ===
            if not skip_qc:
                qc_report = self.generate_qc_report(subject, session, out_dir, outputs)
                outputs["qc_report"] = qc_report
            
            outputs["status"] = "success"
            log.info(f"\n✓ Successfully processed {subject}")
            
        except Exception as e:
            outputs["status"] = "failed"
            outputs["error"] = str(e)
            log.error(f"\n✗ Failed to process {subject}: {e}")
            raise
        
        return outputs
    
    def process_all_subjects(self, subjects: List[str] = None,
                              steps: List[str] = None, 
                              skip_qc: bool = False) -> List[Dict[str, Any]]:
        """Process multiple subjects."""
        
        if subjects is None:
            subjects = self.get_subjects()
        
        log.info(f"Processing {len(subjects)} subjects: {subjects}")
        
        results = []
        for subject in subjects:
            try:
                result = self.process_subject(subject, steps=steps, skip_qc=skip_qc)
                results.append(result)
            except Exception as e:
                log.error(f"Failed to process {subject}: {e}")
                results.append({
                    "subject": subject,
                    "status": "failed",
                    "error": str(e)
                })
        
        # Summary
        success = sum(1 for r in results if r.get("status") == "success")
        failed = len(results) - success
        
        log.info(f"\n{'='*60}")
        log.info(f"PIPELINE COMPLETE")
        log.info(f"  Successful: {success}/{len(results)}")
        log.info(f"  Failed: {failed}/{len(results)}")
        log.info(f"{'='*60}")
        
        return results


def main():
    parser = argparse.ArgumentParser(
        description="Mouse Brain DWI Processing Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process specific subjects
  python run_pipeline.py --subjects sub-596-C-OO-NmM sub-769-C-OO-NmM

  # Process all subjects
  python run_pipeline.py --all

  # Run only preprocessing
  python run_pipeline.py --subjects sub-596-C-OO-NmM --steps preproc

  # Run only tractography (requires preprocessing done)
  python run_pipeline.py --subjects sub-596-C-OO-NmM --steps tractography

  # Skip QC report generation
  python run_pipeline.py --subjects sub-596-C-OO-NmM --skip-qc

  # Verbose output
  python run_pipeline.py --subjects sub-596-C-OO-NmM -v
        """
    )
    
    parser.add_argument(
        "--subjects", "-s",
        nargs="+",
        help="Subject IDs to process"
    )
    parser.add_argument(
        "--all", "-a",
        action="store_true",
        help="Process all subjects in BIDS directory"
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=["preproc", "tractography"],
        default=["preproc", "tractography"],
        help="Pipeline steps to run (default: all)"
    )
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        help="Skip QC report generation"
    )
    parser.add_argument(
        "--base-dir", "-b",
        default="/Users/glerma/soft/mbrain",
        help="Base directory for the project"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize pipeline
    config = PipelineConfig(base_dir=args.base_dir)
    runner = PipelineRunner(config)
    
    # Determine subjects
    if args.all:
        subjects = runner.get_subjects()
    elif args.subjects:
        subjects = args.subjects
    else:
        parser.error("Must specify --subjects or --all")
    
    # Run pipeline
    results = runner.process_all_subjects(
        subjects=subjects,
        steps=args.steps,
        skip_qc=args.skip_qc
    )
    
    # Exit with appropriate code
    failed = sum(1 for r in results if r.get("status") != "success")
    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()
