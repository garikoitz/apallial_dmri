"""
MRtrix3 preprocessing pipeline for mouse brain DTI data.
"""

import subprocess
from pathlib import Path
from typing import Optional, List
import logging

from .config import Config

logger = logging.getLogger(__name__)


class DWIPreprocessor:
    """Handles DWI preprocessing using MRtrix3."""
    def _step_gradcheck(self, subject_id, session_id, current_file, output_dir, basename, dry_run):
        if not current_file.exists():
            logger.error(f"    ❌ Input file not found: {current_file}")
            logger.error("    Run previous steps first or check file paths")
            return current_file
        logger.info("    Step 1: Checking gradient orientation...")
        self._run_command([
            'dwigradcheck',
            str(current_file),
            '-export_grad_mrtrix', 'raw.b',
            '-force && mrconvert', str(current_file), '-grad raw.b {difm}.mif'
        ], subject_id=subject_id, step='gradcheck')
        self.stats['steps_completed'].append(f"gradcheck_{subject_id}_{session_id}")
        return current_file

    def _step_denoise(self, subject_id, session_id, current_file, output_dir, basename, dry_run):
        denoised_file = output_dir / f"{basename}_denoised.mif"
        if denoised_file.exists():
            logger.info("    Step 2: Denoised file already exists, skipping")
            return denoised_file
        if not current_file.exists():
            logger.error(f"    ❌ Input file not found: {current_file}")
            logger.error("    Run previous steps first or check file paths")
            return current_file
        noise_map = output_dir / f"{basename}_noise.mif"
        logger.info("    Step 2: Denoising...")
        self._run_command([
            'dwidenoise',
            str(current_file),
            str(denoised_file),
            '-noise', str(noise_map),
            '-force'
        ], subject_id=subject_id, step='denoise')
        self.stats['steps_completed'].append(f"denoise_{subject_id}_{session_id}")
        return denoised_file

    def _step_unring(self, subject_id, session_id, current_file, output_dir, basename, dry_run):
        unringed_file = output_dir / f"{basename}_unringed.mif"
        if unringed_file.exists():
            logger.info("    Step 3: Unringed file already exists, skipping")
            return unringed_file
        if not current_file.exists():
            logger.error(f"    ❌ Input file not found: {current_file}")
            logger.error("    Run previous steps first or check file paths")
            return current_file
        logger.info("    Step 3: Unringing (Gibbs artifact removal)...")
        self._run_command([
            'mrdegibbs',
            str(current_file),
            str(unringed_file),
            '-force'
        ], subject_id=subject_id, step='unring')
        self.stats['steps_completed'].append(f"unring_{subject_id}_{session_id}")
        return unringed_file

    def _step_mask(self, subject_id, session_id, current_file, output_dir, basename, dry_run):
        mask_file = output_dir / f"{basename}_mask.mif"
        if mask_file.exists():
            logger.info("    Step 4: Mask file already exists, skipping")
            return
        if not current_file.exists():
            logger.error(f"    ❌ Input file not found: {current_file}")
            logger.error("    Run previous steps first or check file paths")
            return
        logger.info("    Step 4: Brain mask estimation...")
        try:
            self._run_command([
                'dwi2mask',
                str(current_file),
                str(mask_file),
                '-force'
            ], subject_id=subject_id, step='mask')
            self.stats['steps_completed'].append(f"mask_{subject_id}_{session_id}")
        except subprocess.CalledProcessError as e:
            logger.warning(f"    dwi2mask failed (expected for ex-vivo/mouse data): {e}")
            warning_msg = f"dwi2mask failed for {subject_id}/{session_id} (using fallback mask)"
            self.stats['warnings'].append(warning_msg)
            if subject_id not in self.stats['subjects_with_issues']:
                self.stats['subjects_with_issues'][subject_id] = []
            self.stats['subjects_with_issues'][subject_id].append(warning_msg)
            logger.info("    Creating simple threshold-based mask instead...")
            self._create_simple_mask(current_file, mask_file, subject_id)

    def _step_bias_correct(self, subject_id, session_id, current_file, output_dir, basename, dry_run):
        corrected_file = output_dir / f"{basename}_preprocessed.mif"
        mask_file = output_dir / f"{basename}_mask.mif"
        if corrected_file.exists():
            logger.info("    Step 5: Bias corrected file already exists, skipping")
            return corrected_file
        if not current_file.exists():
            logger.error(f"    ❌ Input file not found: {current_file}")
            logger.error("    Run previous steps first or check file paths")
            return current_file
        logger.info("    Step 5: Bias field correction...")
        if mask_file.exists():
            self._run_command([
                'dwibiascorrect', 'ants',
                str(current_file),
                str(corrected_file),
                '-mask', str(mask_file),
                '-force'
            ], subject_id=subject_id, step='bias_correct')
        else:
            logger.warning("    No mask available, running bias correction without mask")
            self._run_command([
                'dwibiascorrect', 'ants',
                str(current_file),
                str(corrected_file),
                '-force'
            ], subject_id=subject_id, step='bias_correct')
        self.stats['steps_completed'].append(f"bias_correct_{subject_id}_{session_id}")
        return corrected_file

    def __init__(self, bids_dir: Optional[Path] = None, derivatives_dir: Optional[Path] = None):
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
        }

    def _preprocess_dwi(self, subject_id: str, 
                        session_id: str, 
                        dwi_file: Path, 
                        steps: List[str], 
                        dry_run: bool = False,
                        input_dwi_name="_dwi_flip.mif"):
        """
        Preprocess a single DWI acquisition.
        """
        if steps is None:
            steps = ['gradcheck', 'denoise', 'unring', 'mask', 'bias_correct']
        output_dir = self.preproc_dir / f"sub-{subject_id}" / f"ses-{session_id}" / "dwi"
        output_dir.mkdir(parents=True, exist_ok=True)
        basename = dwi_file.name.replace(input_dwi_name, '')
        logger.info(f"  Processing: {dwi_file.name}")
        current_file = dwi_file

        if 'gradcheck' in steps:
            current_file = self._step_gradcheck(subject_id, session_id, current_file, output_dir, basename, dry_run)
        if 'denoise' in steps:
            current_file = self._step_denoise(subject_id, session_id, current_file, output_dir, basename, dry_run)
        if 'unring' in steps:
            current_file = self._step_unring(subject_id, session_id, current_file, output_dir, basename, dry_run)
        if 'mask' in steps:
            self._step_mask(subject_id, session_id, current_file, output_dir, basename, dry_run)
        if 'bias_correct' in steps:
            current_file = self._step_bias_correct(subject_id, session_id, current_file, output_dir, basename, dry_run)

        # ...existing code...
        final_nifti = output_dir / f"{basename}_preprocessed.nii.gz"
        logger.info(f"    Converting to NIfTI: {final_nifti.name}")
        self._run_command([
            'mrconvert',
            str(current_file),
            str(final_nifti),
            '-force'
        ], subject_id=subject_id, step='export_nifti')
        logger.info(f"    ✓ Preprocessing complete for {dwi_file.name}")

    
    def _preprocess_subject(self, 
                            derivatives_dir: str,
                            subject_id: str, 
                            steps: List[str], 
                            dry_run: bool = False,
                            input_dwi_name="_dwi_flip.mif"):
        """
        Preprocess a single subject.
        
        Parameters
        ----------
        subject_id : str
            Subject ID
        steps : List[str]
            Preprocessing steps to run
        dry_run : bool, optional
            If True, show commands without executing them
        """
        logger.info(f"Processing subject: sub-{subject_id}")
        
        # Find all sessions for this subject
        subject_dir = Path(derivatives_dir) / f"sub-{subject_id}"
        session_dirs = sorted([d for d in subject_dir.iterdir() 
                              if d.is_dir() and d.name.startswith('ses-')])
        
        for session_dir in session_dirs:
            session_id = session_dir.name.replace('ses-', '')
            
            # Check if DWI data exists
            dwi_dir = session_dir / "dwi"
            if not dwi_dir.exists():
                logger.warning(f"No DWI data for sub-{subject_id} ses-{session_id}")
                continue
            
            cmd = f'For each subject/session, find the DWI file ending with {input_dwi_name} in '
            print(f"{cmd}{dwi_dir}...")
            dwi_files = list(dwi_dir.glob(f"*{input_dwi_name}"))
            print(f"Files found: {dwi_files}...")
            if not dwi_files:
                logger.warning(f"No DWI files ending with {input_dwi_name} found in {dwi_dir}")
                continue
            for dwi_file in dwi_files:
                # Now dwi_file is the correct .mif file to use for further processing
                self._preprocess_dwi(subject_id, session_id, dwi_file, steps, dry_run=dry_run)

            if not dwi_files:
                logger.warning(f"No DWI files found in {dwi_dir}")
                continue


    def _preprocess_dwi(self, subject_id: str, 
                        session_id: str, 
                        dwi_file: Path, 
                        steps: List[str], 
                        dry_run: bool = False,
                        input_dwi_name="_dwi_flip.mif"):
        """
        Preprocess a single DWI acquisition.
        
        Parameters
        ----------
        subject_id : str
            Subject ID
        session_id : str
            Session ID
        dwi_file : Path
            Path to DWI NIfTI file
        steps : List[str]
            Preprocessing steps to run
        dry_run : bool, optional
            If True, show commands without executing them
        """
        # Ensure steps is a list
        if steps is None:
            steps = ['gradcheck', 'denoise', 
                     'unring', 'mask', 'bias_correct']
        # , 'dwi2response', 'dwi2fod', 'tckgen'
        # Create output directory for this subject/session
        output_dir = self.preproc_dir / f"sub-{subject_id}" / f"ses-{session_id}" / "dwi"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get basename for output files
        basename = dwi_file.name.replace(input_dwi_name, '')
        
        logger.info(f"  Processing: {dwi_file.name}")
        
        # Track the current file through the pipeline
        current_file = dwi_file
        
        
        # Step 1: Gradient orientation check
        if 'gradcheck' in steps:
            if not current_file.exists():
                logger.error(f"    ❌ Input file not found: {current_file}")
                logger.error("    Run previous steps first or check file paths")
                return
            
            logger.info("    Step 1: Checking gradient orientation...")
            
            self._run_command([
                'dwigradcheck',
                str(current_file),
                '-export_grad_mrtrix', 'raw.b',
                '-force && mrconvert', str(current_file), '-grad raw.b {difm}.mif'
            ], subject_id=subject_id, step='gradcheck')
            
            
            
            self.stats['steps_completed'].append(f"gradcheck_{subject_id}_{session_id}")
        
        # Step 2: Denoising
        if 'denoise' in steps:
            denoised_file = output_dir / f"{basename}_denoised.mif"
            
            # Check if already exists
            if denoised_file.exists():
                logger.info("    Step 2: Denoised file already exists, skipping")
                current_file = denoised_file
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                noise_map = output_dir / f"{basename}_noise.mif"
                logger.info("    Step 2: Denoising...")
                self._run_command([
                    'dwidenoise',
                    str(current_file),
                    str(denoised_file),
                    '-noise', str(noise_map),
                    '-force'
                ], subject_id=subject_id, step='denoise')
                current_file = denoised_file
                self.stats['steps_completed'].append(f"denoise_{subject_id}_{session_id}")
        
        # Step 3: Remove Gibbs ringing artifacts
        if 'unring' in steps:
            unringed_file = output_dir / f"{basename}_unringed.mif"
            
            # Check if already exists
            if unringed_file.exists():
                logger.info("    Step 3: Unringed file already exists, skipping")
                current_file = unringed_file
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                logger.info("    Step 3: Unringing (Gibbs artifact removal)...")
                self._run_command([
                    'mrdegibbs',
                    str(current_file),
                    str(unringed_file),
                    '-force'
                ], subject_id=subject_id, step='unring')
                current_file = unringed_file
                self.stats['steps_completed'].append(f"unring_{subject_id}_{session_id}")
        
        # Step 4: Create brain mask
        if 'mask' in steps:
            mask_file = output_dir / f"{basename}_mask.mif"
            
            # Check if already exists
            if mask_file.exists():
                logger.info("    Step 4: Mask file already exists, skipping")
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                logger.info("    Step 4: Brain mask estimation...")
                # Note: dwi2mask may not work well for ex-vivo or mouse data
                # We'll try it, but may need manual correction
                try:
                    self._run_command([
                        'dwi2mask',
                        str(current_file),
                        str(mask_file),
                        '-force'
                    ], subject_id=subject_id, step='mask')
                    self.stats['steps_completed'].append(f"mask_{subject_id}_{session_id}")
                except subprocess.CalledProcessError as e:
                    logger.warning(f"    dwi2mask failed (expected for ex-vivo/mouse data): {e}")
                    warning_msg = f"dwi2mask failed for {subject_id}/{session_id} (using fallback mask)"
                    self.stats['warnings'].append(warning_msg)
                    if subject_id not in self.stats['subjects_with_issues']:
                        self.stats['subjects_with_issues'][subject_id] = []
                self.stats['subjects_with_issues'][subject_id].append(warning_msg)
                logger.info("    Creating simple threshold-based mask instead...")
                # Create a simple threshold-based mask as fallback
                self._create_simple_mask(current_file, mask_file, subject_id)
        
        # Step 5: Bias field correction
        if 'bias_correct' in steps:
            corrected_file = output_dir / f"{basename}_preprocessed.mif"
            mask_file = output_dir / f"{basename}_mask.mif"
            
            # Check if already exists
            if corrected_file.exists():
                logger.info("    Step 5: Bias corrected file already exists, skipping")
                current_file = corrected_file
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                logger.info("    Step 5: Bias field correction...")
                if mask_file.exists():
                    self._run_command([
                        'dwibiascorrect', 'ants',
                        str(current_file),
                        str(corrected_file),
                        '-mask', str(mask_file),
                        '-force'
                    ], subject_id=subject_id, step='bias_correct')
                else:
                    logger.warning("    No mask available, running bias correction without mask")
                    self._run_command([
                        'dwibiascorrect', 'ants',
                        str(current_file),
                        str(corrected_file),
                        '-force'
                    ], subject_id=subject_id, step='bias_correct')
                current_file = corrected_file
                self.stats['steps_completed'].append(f"bias_correct_{subject_id}_{session_id}")
        
        # Step 6: Two-stage coregistration with anatomical T2w
        # Stage 1: DWI → lowres T2w (matching resolution)
        # Stage 2: lowres T2w → highres T2w (for high-res ROI drawing)
        # Creates transform chain: highres → lowres → DWI for ROI mapping
        if 'coreg' in steps:
            # Find anatomical T2w files
            subject_dir = self.bids_dir / f"sub-{subject_id}"
            session_dir = subject_dir / f"ses-{session_id}"
            anat_dir = session_dir / "anat"
            
            t2w_lowres = None
            t2w_highres = None
            
            if anat_dir.exists():
                lowres_files = list(anat_dir.glob("*_acq-lowres_T2w.nii.gz"))
                highres_files = list(anat_dir.glob("*_acq-highres_T2w.nii.gz"))
                
                if lowres_files:
                    t2w_lowres = lowres_files[0]
                if highres_files:
                    t2w_highres = highres_files[0]
            
            if not t2w_lowres:
                logger.warning("    No lowres T2w anatomical found for coregistration, skipping")
            else:
                logger.info("    Step 6: Two-stage coregistration")
                logger.info(f"      Stage 1: DWI → lowres T2w ({t2w_lowres.name})")
                if t2w_highres:
                    logger.info(f"      Stage 2: lowres → highres T2w ({t2w_highres.name})")
                    logger.info("      → ROIs can be drawn in highres and mapped to DWI space")
                
                # Files for coregistration
                coreg_file = output_dir / f"{basename}_coreg.mif"
                lowres_mif = output_dir / f"{basename}_t2w_lowres.mif"
                highres_mif = output_dir / f"{basename}_t2w_highres.mif" if t2w_highres else None
                
                # Transforms for ROI mapping
                dwi_to_lowres_mrtrix = output_dir / f"{basename}_dwi2lowres.txt"
                lowres_to_highres_mrtrix = output_dir / f"{basename}_lowres2highres.txt" if t2w_highres else None
                highres_to_dwi_mrtrix = output_dir / f"{basename}_highres2dwi.txt" if t2w_highres else None
                
                # Check if already exists
                if coreg_file.exists() and not dry_run:
                    logger.info("    Step 6: Coregistered file already exists, skipping")
                    current_file = coreg_file
                else:
                    if not current_file.exists():
                        logger.error(f"    ❌ Input file not found: {current_file}")
                        logger.error("    Run previous steps first or check file paths")
                        return
                    
                    if coreg_file.exists():
                        logger.info("    Step 6: Coregistered file already exists (dry-run: showing commands)")
                    
                    # ============================================================
                    # STAGE 1: DWI → lowres T2w
                    # ============================================================
                    logger.info(f"    Stage 1: DWI → lowres T2w...")
                    
                    # Convert lowres T2w to MIF
                    self._run_command([
                        'mrconvert',
                        str(t2w_lowres),
                        str(lowres_mif),
                        '-force'
                    ], subject_id=subject_id, step='coreg_convert_lowres', dry_run=dry_run)
                    
                    # Extract b0 from DWI for registration
                    b0_file = output_dir / f"{basename}_b0.mif"
                    logger.info(f"      Extracting mean b0...")
                    self._run_command([
                        'dwiextract',
                        str(current_file),
                        '-bzero', '-',
                        '|', 'mrmath', '-', 'mean',
                        str(b0_file),
                        '-axis', '3',
                        '-force'
                    ], subject_id=subject_id, step='coreg_extract_b0', dry_run=dry_run)
                    
                    # Convert b0 to NIfTI for ANTs
                    b0_nii = output_dir / f"{basename}_b0.nii.gz"
                    self._run_command([
                        'mrconvert',
                        str(b0_file),
                        str(b0_nii),
                        '-force'
                    ], subject_id=subject_id, step='coreg_b0_to_nii', dry_run=dry_run)
                    
                    # ANTs rigid registration: DWI → lowres
                    logger.info(f"      Computing DWI→lowres registration...")
                    transform_prefix_dwi = output_dir / f"{basename}_dwi2lowres_ants"
                    self._run_command([
                        'antsRegistrationSyN.sh',
                        '-d', '3',
                        '-t', 'r',  # Rigid transformation
                        '-f', str(t2w_lowres),  # Fixed: lowres T2w
                        '-m', str(b0_nii),       # Moving: DWI b0
                        '-o', str(transform_prefix_dwi)
                    ], subject_id=subject_id, step='coreg_dwi2lowres_register', dry_run=dry_run)
                    
                    # Convert transform to MRtrix format
                    transform_mat_dwi = str(transform_prefix_dwi) + '0GenericAffine.mat'
                    transform_txt_dwi = output_dir / f"{basename}_dwi2lowres_ants0GenericAffine.txt"
                    
                    self._run_command([
                        'ConvertTransformFile',
                        '3',
                        str(transform_mat_dwi),
                        str(transform_txt_dwi)
                    ], subject_id=subject_id, step='coreg_dwi2lowres_convert', dry_run=dry_run)
                    
                    self._run_command([
                        'transformconvert',
                        str(transform_txt_dwi),
                        'itk_import',
                        str(dwi_to_lowres_mrtrix),
                        '-force'
                    ], subject_id=subject_id, step='coreg_dwi2lowres_import', dry_run=dry_run)
                    
                    # Apply DWI → lowres transformation
                    # Note: -reorient_fod no because this is raw DWI, not FOD
                    self._run_command([
                        'mrtransform',
                        str(current_file),
                        str(coreg_file),
                        '-linear', str(dwi_to_lowres_mrtrix),
                        '-reorient_fod', 'no',  # Not an FOD image, just raw DWI
                        '-force'
                    ], subject_id=subject_id, step='coreg_dwi2lowres_apply', dry_run=dry_run)
                    
                    logger.info(f"      ✓ DWI→lowres complete")
                    
                    # ============================================================
                    # STAGE 2: lowres → highres T2w (if highres exists)
                    # ============================================================
                    if t2w_highres:
                        logger.info(f"    Stage 2: lowres → highres T2w...")
                        
                        # Convert highres T2w to MIF
                        self._run_command([
                            'mrconvert',
                            str(t2w_highres),
                            str(highres_mif),
                            '-force'
                        ], subject_id=subject_id, step='coreg_convert_highres', dry_run=dry_run)
                        
                        # ANTs rigid registration: lowres → highres
                        logger.info(f"      Computing lowres→highres registration...")
                        transform_prefix_anat = output_dir / f"{basename}_lowres2highres_ants"
                        self._run_command([
                            'antsRegistrationSyN.sh',
                            '-d', '3',
                            '-t', 'r',  # Rigid transformation
                            '-f', str(t2w_highres),  # Fixed: highres T2w
                            '-m', str(t2w_lowres),   # Moving: lowres T2w
                            '-o', str(transform_prefix_anat)
                        ], subject_id=subject_id, step='coreg_lowres2highres_register', dry_run=dry_run)
                        
                        # Convert transform to MRtrix format
                        transform_mat_anat = str(transform_prefix_anat) + '0GenericAffine.mat'
                        transform_txt_anat = output_dir / f"{basename}_lowres2highres_ants0GenericAffine.txt"
                        
                        self._run_command([
                            'ConvertTransformFile',
                            '3',
                            str(transform_mat_anat),
                            str(transform_txt_anat)
                        ], subject_id=subject_id, step='coreg_lowres2highres_convert', dry_run=dry_run)
                        
                        self._run_command([
                            'transformconvert',
                            str(transform_txt_anat),
                            'itk_import',
                            str(lowres_to_highres_mrtrix),
                            '-force'
                        ], subject_id=subject_id, step='coreg_lowres2highres_import', dry_run=dry_run)
                        
                        # Create inverse transform: highres → lowres
                        highres_to_lowres_mrtrix = output_dir / f"{basename}_highres2lowres.txt"
                        self._run_command([
                            'transformcalc',
                            str(lowres_to_highres_mrtrix),
                            'invert',
                            str(highres_to_lowres_mrtrix),
                            '-force'
                        ], subject_id=subject_id, step='coreg_invert_highres2lowres', dry_run=dry_run)
                        
                        # Compose full chain: highres → lowres → DWI
                        logger.info(f"      Creating transform chain: highres → lowres → DWI...")
                        self._run_command([
                            'transformcompose',
                            str(highres_to_lowres_mrtrix),
                            str(dwi_to_lowres_mrtrix),
                            str(highres_to_dwi_mrtrix),
                            '-force'
                        ], subject_id=subject_id, step='coreg_compose_highres2dwi', dry_run=dry_run)
                        
                        logger.info(f"      ✓ lowres→highres complete")
                        logger.info(f"      ✓ Transform chain created for ROI mapping")
                        logger.info(f"        → Use {highres_to_dwi_mrtrix.name} to map ROIs from highres to DWI")
                    
                    current_file = coreg_file
                    if not dry_run:
                        self.stats['steps_completed'].append(f"coreg_{subject_id}_{session_id}")
                    
                    logger.info(f"    ✓ Coregistration complete")
        
        # Step 7: Atlas-based coregistration (Allen Mouse Brain CCF)
        if 'atlas_coreg' in steps:
            logger.info("    Step 7: Atlas-based coregistration")
            
            # Check if atlas is downloaded
            from .atlas_manager import AllenAtlasManager
            atlas_manager = AllenAtlasManager()
            
            if not atlas_manager.is_downloaded():
                logger.info("      Downloading Allen Mouse Brain CCF atlas...")
                atlas_manager.download()
            
            atlas_files = atlas_manager.get_atlas_files()
            atlas_template = atlas_files['template']
            
            # Create atlas-specific output directory
            atlas_output_dir = output_dir / 'atlas'
            atlas_output_dir.mkdir(exist_ok=True)
            
            # Find T2w reference (prefer highres, fallback to lowres)
            t2w_ref = output_dir / f"{basename}_t2w_highres.mif"
            if not t2w_ref.exists():
                t2w_ref = output_dir / f"{basename}_t2w_lowres.mif"
            
            if not t2w_ref.exists():
                logger.warning("      No T2w reference found. Run coregistration step first.")
                logger.warning("      Skipping atlas coregistration")
            else:
                logger.info(f"      Using T2w reference: {t2w_ref.name}")
                logger.info("      Preparing T2w for atlas registration...")
                
                # Convert template to MIF for MRtrix3 compatibility
                atlas_template_mif = atlas_output_dir / 'allen_template_25um.mif'
                if not atlas_template_mif.exists() and not dry_run:
                    self._run_command([
                        'mrconvert',
                        str(atlas_template),
                        str(atlas_template_mif),
                        '-force'
                    ], subject_id=subject_id, step='atlas_convert_template')
                
                # Convert T2w to NIfTI for ANTs processing
                t2w_nii = atlas_output_dir / f"{basename}_t2w_ref.nii.gz"
                if not t2w_nii.exists() or dry_run:
                    self._run_command([
                        'mrconvert',
                        str(t2w_ref),
                        str(t2w_nii),
                        '-force'
                    ], subject_id=subject_id, step='atlas_convert_t2w')
                
                # Step 7.1: N4 bias field correction on T2w
                t2w_n4 = atlas_output_dir / f"{basename}_t2w_n4.nii.gz"
                if not t2w_n4.exists() or dry_run:
                    logger.info("      Step 7.1: N4 bias field correction on T2w...")
                    self._run_command([
                        'N4BiasFieldCorrection',
                        '-d', '3',
                        '-i', str(t2w_nii),
                        '-o', str(t2w_n4)
                    ], subject_id=subject_id, step='atlas_n4_bias_correct')
                
                # Step 7.2: Brain extraction using simplified template-based approach
                # Note: antsBrainExtraction.sh consistently segfaults during Atropos segmentation
                # This simplified approach avoids Atropos while still using Allen atlas as prior
                t2w_brain = atlas_output_dir / f"{basename}_BrainExtractionBrain.nii.gz"
                t2w_brainmask = atlas_output_dir / f"{basename}_BrainExtractionMask.nii.gz"
                
                if not t2w_brain.exists() or dry_run:
                    logger.info("      Step 7.2: Brain extraction using Allen atlas prior (simplified approach)...")
                    
                    # Ensure atlas NIfTI exists with correct orientation
                    atlas_nii = atlas_output_dir / 'allen_template_25um.nii.gz'
                    if not atlas_nii.exists() or dry_run:
                        self._run_command([
                            'mrconvert',
                            str(atlas_template),
                            str(atlas_nii),
                            '-strides', '1,2,3',
                            '-force'
                        ], subject_id=subject_id, step='atlas_convert_to_nii')
                    
                    # Create binary mask from Allen template
                    atlas_mask = atlas_output_dir / 'allen_template_mask.nii.gz'
                    if not atlas_mask.exists() or dry_run:
                        logger.info("        Creating brain mask from Allen template...")
                        self._run_command([
                            'ThresholdImage',
                            '3',
                            str(atlas_nii),
                            str(atlas_mask),
                            '1', 'Inf',
                            '1', '0'
                        ], subject_id=subject_id, step='atlas_create_mask')
                        
                        # Dilate mask to ensure we don't cut off brain tissue
                        self._run_command([
                            'ImageMath',
                            '3',
                            str(atlas_mask),
                            'MD',  # Morphological dilation
                            str(atlas_mask),
                            '3'  # Dilate by 3 voxels
                        ], subject_id=subject_id, step='atlas_dilate_mask')
                    
                    # Quick affine registration: template → subject
                    logger.info("        Registering Allen atlas to subject for mask warping...")
                    template_to_subject_prefix = atlas_output_dir / f"{basename}_template2subject_"
                    if not Path(str(template_to_subject_prefix) + '0GenericAffine.mat').exists() or dry_run:
                        self._run_command([
                            'antsRegistrationSyN.sh',
                            '-d', '3',
                            '-t', 'a',  # Affine only
                            '-f', str(t2w_n4),
                            '-m', str(atlas_nii),
                            '-o', str(template_to_subject_prefix)
                        ], subject_id=subject_id, step='atlas_register_for_mask')
                    
                    # Warp dilated mask to subject space
                    logger.info("        Applying mask to subject...")
                    if not t2w_brainmask.exists() or dry_run:
                        self._run_command([
                            'antsApplyTransforms',
                            '-d', '3',
                            '-i', str(atlas_mask),
                            '-r', str(t2w_n4),
                            '-t', str(template_to_subject_prefix) + '0GenericAffine.mat',
                            '-o', str(t2w_brainmask),
                            '-n', 'NearestNeighbor'
                        ], subject_id=subject_id, step='atlas_warp_mask')
                    
                    # Apply mask to create brain-extracted image
                    if not t2w_brain.exists() or dry_run:
                        self._run_command([
                            'ImageMath',
                            '3',
                            str(t2w_brain),
                            'm',
                            str(t2w_n4),
                            str(t2w_brainmask)
                        ], subject_id=subject_id, step='atlas_apply_mask')
                
                # Use brain-extracted T2w for registration
                t2w_for_reg = t2w_brain if t2w_brain.exists() or dry_run else t2w_n4
                
                # Step 7.3: Register subject to atlas using ANTs
                # Fixed: atlas template (target), Moving: subject T2w (source)
                # This gives us forward transform (subject→atlas) and inverse (atlas→subject)
                subject_to_atlas_prefix = atlas_output_dir / f"{basename}_subject2atlas_"
                subject_to_atlas_mat = Path(str(subject_to_atlas_prefix) + '0GenericAffine.mat')
                atlas_to_subject_mat = Path(str(subject_to_atlas_prefix) + '0InverseGenericAffine.mat')
                
                if not subject_to_atlas_mat.exists() or dry_run:
                    logger.info("      Step 7.3: Computing subject → atlas registration (SyNBold - less aggressive non-linear)...")
                    logger.info("      This will generate inverse warp for mapping ROIs to subject...")
                    
                    # Ensure atlas NIfTI exists
                    atlas_nii = atlas_output_dir / 'allen_template_25um.nii.gz'
                    if not atlas_nii.exists() or dry_run:
                        self._run_command([
                            'mrconvert',
                            str(atlas_template),
                            str(atlas_nii),
                            '-force'
                        ], subject_id=subject_id, step='atlas_convert_to_nii')
                    
                    self._run_command([
                        'antsRegistrationSyN.sh',
                        '-d', '3',
                        '-t', 'b',  # SyNBold - less aggressive non-linear than 's'
                        '-f', str(atlas_nii),  # Fixed: atlas template
                        '-m', str(t2w_for_reg),  # Moving: brain-extracted T2w
                        '-o', str(subject_to_atlas_prefix)
                    ], subject_id=subject_id, step='atlas_register')
                
                # Transform ROIs to subject space using ANTs
                # Use antsApplyTransforms with inverse flag to bring ROIs from atlas→subject
                roi_dir = atlas_output_dir / 'rois'
                roi_dir.mkdir(exist_ok=True)
                
                atlas_roi_dir = atlas_manager.atlas_dir / 'rois'
                if not atlas_roi_dir.exists():
                    logger.info("      Extracting atlas ROIs...")
                    atlas_manager.extract_all_rois(atlas_roi_dir)
                
                # Transform each ROI (cortex, striatum, thalamus) with left/right hemispheres
                for roi_base in ['cortex', 'striatum', 'thalamus']:
                    for hemi in ['left', 'right']:
                        roi_name = f'{roi_base}_{hemi}'
                        atlas_roi = atlas_roi_dir / f'atlas_{roi_name}.nii.gz'
                        subject_roi_nii = roi_dir / f"{basename}_{roi_name}.nii.gz"
                        subject_roi_mif = roi_dir / f"{basename}_{roi_name}.mif"
                        
                        if atlas_roi.exists() and (not subject_roi_mif.exists() or dry_run):
                            # Use antsApplyTransforms with inverse warp to map atlas→subject
                            logger.info(f"      Transforming {roi_name} ROI to subject space (ANTs inverse warp)...")
                            if not subject_roi_nii.exists() or dry_run:
                                # For SyN, we need to apply inverse warp and affine in correct order
                                warp_file = Path(str(subject_to_atlas_prefix) + '1Warp.nii.gz')
                                inverse_warp = str(warp_file).replace('1Warp', '1InverseWarp')
                                self._run_command([
                                    'antsApplyTransforms',
                                    '-d', '3',
                                    '-i', str(atlas_roi),  # Input: atlas ROI
                                    '-r', str(t2w_nii),  # Reference: subject T2w
                                    '-t', f'[{str(subject_to_atlas_mat)},1]',  # Inverse affine (first)
                                    '-t', inverse_warp,  # Inverse warp (second)
                                    '-o', str(subject_roi_nii),
                                    '-n', 'NearestNeighbor'  # Preserve labels
                                ], subject_id=subject_id, step=f'atlas_transform_{roi_name}')
                            
                            # Convert to MIF format
                            if not subject_roi_mif.exists() or dry_run:
                                self._run_command([
                                    'mrconvert',
                                    str(subject_roi_nii),
                                    str(subject_roi_mif),
                                    '-force'
                                ], subject_id=subject_id, step=f'atlas_roi_to_mif_{roi_name}')
                            
                            # Also transform ROI to DWI space (if coregistration was done)
                            highres2dwi = output_dir / f"{basename}_highres2dwi.txt"
                            if highres2dwi.exists():
                                subject_roi_dwi = roi_dir / f"{basename}_{roi_name}_dwi.mif"
                                if not subject_roi_dwi.exists() or dry_run:
                                    logger.info(f"      Transforming {roi_name} ROI to DWI space...")
                                    self._run_command([
                                        'mrtransform',
                                        str(subject_roi_mif),
                                        str(subject_roi_dwi),
                                        '-linear', str(highres2dwi),
                                        '-inverse',  # Apply inverse transform (highres→DWI)
                                        '-template', str(current_file),  # Use current DWI as template
                                        '-interp', 'nearest',  # Preserve labels
                                        '-force'
                                    ], subject_id=subject_id, step=f'atlas_roi_to_dwi_{roi_name}')
                
                logger.info(f"    ✓ Atlas coregistration complete")
                logger.info(f"      Brain-extracted T2w: {t2w_brain.name if t2w_brain.exists() else 'N/A'}")
                logger.info(f"      ROIs available in subject space: {roi_dir}")
                if not dry_run:
                    self.stats['steps_completed'].append(f"atlas_coreg_{subject_id}_{session_id}")
        
        # Step 8: Response function estimation
        if 'dwi2response' in steps:
            response_wm = output_dir / f"{basename}_response_wm.txt"
            response_gm = output_dir / f"{basename}_response_gm.txt"
            response_csf = output_dir / f"{basename}_response_csf.txt"
            
            # Check if already exists
            if response_wm.exists():
                logger.info("    Step 8: Response functions already exist, skipping")
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                logger.info("    Step 8: Estimating response functions...")
                mask_file = output_dir / f"{basename}_mask.mif"
                
                if mask_file.exists():
                    self._run_command([
                        'dwi2response', 'dhollander',
                        str(current_file),
                        str(response_wm),
                        str(response_gm),
                        str(response_csf),
                        '-mask', str(mask_file),
                        '-force'
                    ], subject_id=subject_id, step='dwi2response')
                else:
                    logger.warning("    No mask available, running without mask")
                    self._run_command([
                        'dwi2response', 'dhollander',
                        str(current_file),
                        str(response_wm),
                        str(response_gm),
                        str(response_csf),
                        '-force'
                    ], subject_id=subject_id, step='dwi2response')
                
                self.stats['steps_completed'].append(f"dwi2response_{subject_id}_{session_id}")
        
        # Step 8: Fiber Orientation Distribution estimation
        if 'dwi2fod' in steps:
            fod_wm = output_dir / f"{basename}_fod_wm.mif"
            fod_gm = output_dir / f"{basename}_fod_gm.mif"
            fod_csf = output_dir / f"{basename}_fod_csf.mif"
            
            # Check if already exists
            if fod_wm.exists():
                logger.info("    Step 9: FOD files already exist, skipping")
            else:
                if not current_file.exists():
                    logger.error(f"    ❌ Input file not found: {current_file}")
                    logger.error("    Run previous steps first or check file paths")
                    return
                
                # Check for response functions
                response_wm = output_dir / f"{basename}_response_wm.txt"
                response_gm = output_dir / f"{basename}_response_gm.txt"
                response_csf = output_dir / f"{basename}_response_csf.txt"
                
                if not all([response_wm.exists(), response_gm.exists(), response_csf.exists()]):
                    logger.error("    ❌ Response functions not found")
                    logger.error("    Run dwi2response step first")
                    return
                
                logger.info("    Step 9: Computing fiber orientation distributions...")
                mask_file = output_dir / f"{basename}_mask.mif"
                
                if mask_file.exists():
                    self._run_command([
                        'dwi2fod', 'msmt_csd',
                        str(current_file),
                        str(response_wm), str(fod_wm),
                        str(response_gm), str(fod_gm),
                        str(response_csf), str(fod_csf),
                        '-mask', str(mask_file),
                        '-force'
                    ], subject_id=subject_id, step='dwi2fod')
                else:
                    logger.warning("    No mask available, running without mask")
                    self._run_command([
                        'dwi2fod', 'msmt_csd',
                        str(current_file),
                        str(response_wm), str(fod_wm),
                        str(response_gm), str(fod_gm),
                        str(response_csf), str(fod_csf),
                        '-force'
                    ], subject_id=subject_id, step='dwi2fod')
                
                self.stats['steps_completed'].append(f"dwi2fod_{subject_id}_{session_id}")
        
        # Step 10: Tractography
        if 'tckgen' in steps:
            tracks_file = output_dir / f"{basename}_tracks.tck"
            
            # Check if already exists (skip only in non-dry-run mode)
            if tracks_file.exists() and not dry_run:
                logger.info("    Step 10: Tractography file already exists, skipping")
            else:
                # Check for FOD file
                fod_wm = output_dir / f"{basename}_fod_wm.mif"
                
                if not fod_wm.exists():
                    logger.error(f"    ❌ WM FOD file not found: {fod_wm}")
                    logger.error("    Run dwi2fod step first")
                    return
                
                if tracks_file.exists():
                    logger.info("    Step 10: Tractography file already exists (dry-run: showing command)")
                else:
                    logger.info("    Step 9: Performing tractography...")
                
                mask_file = output_dir / f"{basename}_mask.mif"
                
                # Mouse brain parameters (smaller than human):
                # - Smaller step size (0.02mm) for fine mouse brain structures
                # - Shorter min/max length (1-50mm) appropriate for mouse brain size
                # - More lenient angle (60°) to capture complex mouse white matter
                # - Lower cutoff (0.05) for mouse ex-vivo data
                # - More streamlines (10M) for good coverage of small structures
                
                if mask_file.exists():
                    self._run_command([
                        'tckgen',
                        str(fod_wm),
                        str(tracks_file),
                        '-algorithm', 'iFOD2',
                        '-seed_image', str(mask_file),
                        '-select', '10000000',  # 10M streamlines
                        '-step', '0.02',        # 0.02mm step size for mouse brain
                        '-angle', '60',         # 60° max curvature angle
                        '-minlength', '1',      # 1mm minimum length
                        '-maxlength', '50',     # 50mm maximum length
                        '-cutoff', '0.05',      # Lower cutoff for ex-vivo
                        '-force'
                    ], subject_id=subject_id, step='tckgen', dry_run=dry_run)
                else:
                    logger.warning("    No mask available, using whole brain seeding")
                    self._run_command([
                        'tckgen',
                        str(fod_wm),
                        str(tracks_file),
                        '-algorithm', 'iFOD2',
                        '-select', '10000000',
                        '-step', '0.02',
                        '-angle', '60',
                        '-minlength', '1',
                        '-maxlength', '50',
                        '-cutoff', '0.05',
                        '-force'
                    ], subject_id=subject_id, step='tckgen', dry_run=dry_run)
                
                if not dry_run:
                    self.stats['steps_completed'].append(f"tckgen_{subject_id}_{session_id}")
        
        # Convert final output to NIfTI for compatibility
        final_nifti = output_dir / f"{basename}_preprocessed.nii.gz"
        logger.info(f"    Converting to NIfTI: {final_nifti.name}")
        self._run_command([
            'mrconvert',
            str(current_file),
            str(final_nifti),
            '-force'
        ], subject_id=subject_id, step='export_nifti')
        
        logger.info(f"    ✓ Preprocessing complete for {dwi_file.name}")
    
    def _correct_bvals(self, bval_file: Path, output_file: Path = None):
        """
        Correct b-values: set values < 20 to 0, round values close to 1000 to 1000.
        
        Parameters
        ----------
        bval_file : Path
            Path to input bval file
        output_file : Path, optional
            Path to output bval file. If None, overwrites input file.
        """
        import numpy as np
        
        # Read bvals
        bvals = np.loadtxt(bval_file)
        
        # Correct values
        bvals[bvals < 20] = 0
        bvals[(bvals >= 900) & (bvals <= 1100)] = 1000
        
        # Write to output
        if output_file is None:
            output_file = bval_file
        np.savetxt(output_file, bvals.reshape(1, -1), fmt='%g')
        logger.info(f"    Corrected b-values: {np.unique(bvals)}")
    
    def _create_simple_mask(self, dwi_file: Path, mask_file: Path, subject_id: str = None):
        """
        Create a simple threshold-based mask as fallback.
        
        Parameters
        ----------
        dwi_file : Path
            Input DWI file
        mask_file : Path
            Output mask file
        subject_id : str, optional
            Subject ID for logging
        """
        # Use median of non-zero values as threshold
        self._run_command([
            'mrthreshold',
            str(dwi_file),
            str(mask_file),
            '-abs', '10',  # Minimum threshold
            '-force'
        ], subject_id=subject_id, step='mask_fallback')
    
    def _split_long_string(self, s: str, max_width: int = 75) -> List[str]:
        """
        Split a long string (like a path) into chunks for display.
        Uses quoted strings that bash will concatenate.
        
        Parameters
        ----------
        s : str
            String to split
        max_width : int
            Maximum width per chunk
            
        Returns
        -------
        List[str]
            List of quoted string chunks
        """
        if len(s) <= max_width:
            return [s]
        
        # For paths, try to split at directory boundaries
        if '/' in s:
            parts = []
            current = ''
            segments = s.split('/')
            
            for i, seg in enumerate(segments):
                # Handle absolute paths (first element is empty string)
                if i == 0 and seg == '':
                    current = '/'
                    continue
                
                # Build test string with separator
                if current == '/':
                    test = '/' + seg
                elif current:
                    test = current + '/' + seg
                else:
                    test = seg
                
                # Check if adding this segment would exceed max width
                if len(test) > max_width and current and current != '/':
                    # Save current chunk (add trailing / for clarity)
                    parts.append(f"'{current}/'")
                    current = seg
                else:
                    current = test
            
            # Add final chunk
            if current:
                parts.append(f"'{current}'")
            
            return parts if len(parts) > 1 else [s]
        else:
            # Non-path string, just split at max_width
            chunks = []
            for i in range(0, len(s), max_width):
                chunks.append(f"'{s[i:i+max_width]}'")
            return chunks

    def _format_command_for_display(self, cmd: List[str], max_width: int = 75) -> str:
        """
        Format command with backslashes for multi-line display.
        Split long paths across multiple lines to prevent terminal wrapping.
        Handles piped commands specially to preserve pipeline structure.
        
        Parameters
        ----------
        cmd : List[str]
            Command and arguments as list
        max_width : int
            Maximum line width (default: 75)
        
        Returns
        -------
        str
            Formatted command string with line breaks
        """
        if len(cmd) <= 3:
            # Short commands on one line
            return ' '.join(cmd)
        
        # Check if this is a piped command
        if '|' in cmd:
            # Split into separate commands at pipes
            pipe_commands = []
            current_cmd = []
            for item in cmd:
                if item == '|':
                    if current_cmd:
                        pipe_commands.append(current_cmd)
                        current_cmd = []
                else:
                    current_cmd.append(item)
            if current_cmd:
                pipe_commands.append(current_cmd)
            
            # Format each command separately and join with pipes
            formatted_cmds = []
            for pcmd in pipe_commands:
                formatted_cmds.append(self._format_single_command(pcmd, max_width, indent='  '))
            
            # Join with pipe separators
            result = []
            for i, fcmd in enumerate(formatted_cmds):
                if i == 0:
                    result.append(fcmd)
                else:
                    result.append('  | \\')
                    result.append(fcmd)
            return '\n'.join(result)
        else:
            return self._format_single_command(cmd, max_width)
    
    def _format_single_command(self, cmd: List[str], max_width: int = 75, indent: str = '') -> str:
        """
        Format a single command (no pipes).
        
        Parameters
        ----------
        cmd : List[str]
            Command and arguments as list
        max_width : int
            Maximum line width
        indent : str
            Indentation for continuation lines
        
        Returns
        -------
        str
            Formatted command string
        """
        lines = []
        
        # First line: command name
        lines.append(indent + cmd[0])
        i = 1
        
        # Add positional arguments (before any options)
        while i < len(cmd) and not cmd[i].startswith('-'):
            arg = cmd[i]
            # Keep full paths on single line for easy copy-paste
            lines.append(indent + '  ' + arg)
            i += 1
        
        # Process remaining arguments (options and their values)
        while i < len(cmd):
            if cmd[i].startswith('-'):
                # Check if next argument exists and is a value (not another option)
                if i + 1 < len(cmd) and not cmd[i + 1].startswith('-'):
                    # Option with value
                    value = cmd[i + 1]
                    option_line = indent + '  ' + cmd[i]
                    
                    # Keep full paths on single line for easy copy-paste
                    if len(option_line + ' ' + value) <= max_width + len(indent):
                        lines.append(option_line + ' ' + value)
                    else:
                        # Value on separate line (full path, no splitting)
                        lines.append(option_line)
                        lines.append(indent + '    ' + value)
                    i += 2
                else:
                    # Option without value (flag)
                    lines.append(indent + '  ' + cmd[i])
                    i += 1
            else:
                # Standalone argument
                lines.append(indent + '  ' + cmd[i])
                i += 1
        
        # Join with backslash continuation
        return ' \\\n'.join(lines)
    
    def _run_command(self, cmd: List[str], subject_id: str = None, step: str = None, dry_run: bool = False):
        """
        Run a shell command and log output.
        
        Parameters
        ----------
        cmd : List[str]
            Command and arguments
        subject_id : str, optional
            Subject ID for logging
        step : str, optional
            Step name for logging
        dry_run : bool, optional
            If True, only show command without executing
        """
        cmd_str = ' '.join(cmd)
        
        # Format command with line breaks for readability
        formatted_cmd = self._format_command_for_display(cmd)
        
        # Print command clearly to console
        print(f"\n{'='*80}")
        if subject_id and step:
            print(f"📍 Subject: {subject_id} | Step: {step}")
        print("📋 Copy-paste command:")
        print(f"{formatted_cmd}")
        print(f"{'='*80}")
        
        logger.debug(f"Running: {cmd_str}")
        
        # Track command in stats
        self.stats['commands_run'].append({
            'subject': subject_id,
            'step': step,
            'command': cmd_str
        })
        
        # Skip execution in dry-run mode
        if dry_run or (hasattr(self, 'dry_run') and self.dry_run):
            print("🔍 DRY RUN - Command not executed\n")
            return
        
        # Check if command contains pipes - if so, use shell execution
        if '|' in cmd:
            # For piped commands, use shell=True with the string version
            try:
                result = subprocess.run(
                    cmd_str,
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Print only essential output
                if result.stdout:
                    lines = [line for line in result.stdout.strip().split('\n') 
                            if line and not any(char in line for char in ['%', '\r']) 
                            and not line.strip().startswith('[')]
                    if lines:
                        print(f"✓ Output: {chr(10).join(lines)}")
                
                if result.stderr:
                    lines = [line for line in result.stderr.strip().split('\n')
                            if line and not any(char in line for char in ['%', '\r'])]
                    if lines:
                        print(f"ℹ️  Info: {chr(10).join(lines)}")
                    logger.debug(result.stderr)
                
                print("✅ Completed successfully\n")
                
            except subprocess.CalledProcessError as e:
                print(f"❌ ERROR: Command failed with exit code {e.returncode}\n")
                print("📋 Copy-paste this command to run in terminal:")
                print(f"{cmd_str}\n")
                if e.stdout:
                    print(f"Standard Output:\n{e.stdout}")
                if e.stderr:
                    print(f"Standard Error:\n{e.stderr}")
                raise
            return
            return
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Print only essential output (filter out progress lines)
            if result.stdout:
                # Filter out progress percentage lines
                lines = [line for line in result.stdout.strip().split('\n') 
                        if line and not any(char in line for char in ['%', '\r']) 
                        and not line.strip().startswith('[')]
                if lines:
                    print(f"✓ Output: {chr(10).join(lines)}")
            
            if result.stderr:
                # Filter MRtrix progress messages
                lines = [line for line in result.stderr.strip().split('\n')
                        if line and not any(char in line for char in ['%', '\r'])
                        and 'dwigradcheck:' not in line
                        and not line.strip().startswith('[')]
                if lines:
                    print(f"ℹ️  Info: {chr(10).join(lines)}")
                logger.debug(result.stderr)
            
            print("✅ Completed successfully\n")
            
        except subprocess.CalledProcessError as e:
            print(f"❌ ERROR: Command failed with exit code {e.returncode}\n")
            print("📋 Copy-paste this command to run in terminal:")
            print(f"{cmd_str}\n")
            if e.stdout:
                print(f"Standard Output:\n{e.stdout}")
            if e.stderr:
                print(f"Standard Error:\n{e.stderr}")
            raise    
    def transform_roi_to_dwi(self, subject_id: str, session_id: str, roi_highres: Path, 
                             output_roi: Optional[Path] = None, dry_run: bool = False) -> Path:
        """
        Transform an ROI from highres T2w space to native DWI space.
        
        This uses the pre-computed transformation chain from coregistration:
        highres → lowres → DWI
        
        Parameters
        ----------
        subject_id : str
            Subject ID (without 'sub-' prefix)
        session_id : str
            Session ID (without 'ses-' prefix)
        roi_highres : Path
            Path to ROI mask in highres T2w space (NIfTI or MIF)
        output_roi : Path, optional
            Output path for transformed ROI. If None, saves to derivatives/mrtrix_preproc
        dry_run : bool, optional
            If True, show commands without executing
        
        Returns
        -------
        Path
            Path to transformed ROI in DWI space
        
        Examples
        --------
        >>> preprocessor = DWIPreprocessor()
        >>> # Draw ROI in highres T2w space (e.g., using mrview)
        >>> roi_highres = Path("my_roi_highres.nii.gz")
        >>> # Transform to DWI space for tractography
        >>> roi_dwi = preprocessor.transform_roi_to_dwi("596-C-OO-NmM", "240911103349", 
        ...                                             roi_highres)
        >>> # Use roi_dwi for tracking: tckgen ... -seed_image roi_dwi
        """
        # Find the transform file
        output_dir = self.preproc_dir / f"sub-{subject_id}" / f"ses-{session_id}" / "dwi"
        basename = f"sub-{subject_id}_ses-{session_id}"
        
        transform_file = output_dir / f"{basename}_highres2dwi.txt"
        
        if not transform_file.exists() and not dry_run:
            raise FileNotFoundError(
                f"Transform file not found: {transform_file}\n"
                f"Run coregistration first: python main.py --preprocess --steps coreg"
            )
        
        # Set output path
        if output_roi is None:
            roi_name = Path(roi_highres).stem
            output_roi = output_dir / f"{basename}_{roi_name}_dwi.mif"
        
        logger.info(f"Transforming ROI: {roi_highres.name} → DWI space")
        logger.info(f"  Using transform: {transform_file.name}")
        logger.info(f"  Output: {output_roi.name}")
        
        # Convert ROI to MIF if needed (for easier handling)
        roi_mif = output_dir / "temp_roi_highres.mif"
        if not str(roi_highres).endswith('.mif'):
            self._run_command([
                'mrconvert',
                str(roi_highres),
                str(roi_mif),
                '-force'
            ], subject_id=subject_id, step='roi_convert_to_mif', dry_run=dry_run)
        else:
            roi_mif = roi_highres
        
        # Apply transformation using mrtransform
        # Use nearest-neighbor interpolation for ROI masks to preserve integer labels
        self._run_command([
            'mrtransform',
            str(roi_mif),
            str(output_roi),
            '-linear', str(transform_file),
            '-template', str(output_dir / f"{basename}_raw.mif"),  # Use DWI as template
            '-interp', 'nearest',  # Preserve ROI labels
            '-force'
        ], subject_id=subject_id, step='roi_transform', dry_run=dry_run)
        
        # Clean up temp file
        if not dry_run and roi_mif != roi_highres:
            roi_mif.unlink(missing_ok=True)
        
        logger.info(f"✓ ROI transformed to DWI space: {output_roi}")
        return output_roi