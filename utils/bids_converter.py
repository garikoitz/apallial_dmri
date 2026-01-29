"""
BIDS converter for Flywheel data structure.
"""

import json
import shutil
from pathlib import Path
from typing import Optional, Dict, List
import logging

from .config import Config
from .file_utils import (
    find_nifti_files,
    find_diffusion_files,
    clean_subject_id,
    extract_session_info,
    determine_modality
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FlywheelToBIDSConverter:
    """Converts Flywheel data structure to BIDS format."""
    
    def __init__(self, flywheel_root: Optional[Path] = None, bids_root: Optional[Path] = None):
        """
        Initialize the converter.
        
        Parameters
        ----------
        flywheel_root : Path, optional
            Root directory of Flywheel data. Defaults to Config.FLYWHEEL_DIR
        bids_root : Path, optional
            Output BIDS directory. Defaults to Config.BIDS_DIR
        """
        self.flywheel_root = Path(flywheel_root) if flywheel_root else Config.FLYWHEEL_DIR
        self.bids_root = Path(bids_root) if bids_root else Config.BIDS_DIR
        
        # Track conversion statistics
        self.stats = {
            'subjects': 0,
            'sessions': 0,
            'acquisitions': 0,
            'files_copied': 0,
            'errors': []
        }
    
    def convert(self, dry_run: bool = False) -> Dict:
        """
        Convert all Flywheel data to BIDS format.
        
        Parameters
        ----------
        dry_run : bool
            If True, only report what would be done without copying files
            
        Returns
        -------
        Dict
            Conversion statistics
        """
        logger.info(f"Starting conversion from {self.flywheel_root} to {self.bids_root}")
        logger.info(f"Dry run: {dry_run}")
        
        # Create BIDS root and dataset_description.json
        if not dry_run:
            Config.ensure_directories()
            self._create_dataset_description()
        
        # Find all subjects
        subjects_dir = self.flywheel_root / "SUBJECTS"
        if not subjects_dir.exists():
            logger.error(f"Subjects directory not found: {subjects_dir}")
            return self.stats
        
        subject_dirs = [d for d in subjects_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        logger.info(f"Found {len(subject_dirs)} subjects")
        
        # Process each subject
        for subject_dir in subject_dirs:
            try:
                self._convert_subject(subject_dir, dry_run=dry_run)
                self.stats['subjects'] += 1
            except Exception as e:
                error_msg = f"Error processing subject {subject_dir.name}: {str(e)}"
                logger.error(error_msg)
                self.stats['errors'].append(error_msg)
        
        # Create participants.tsv
        if not dry_run:
            self._create_participants_file()
        
        logger.info("Conversion complete!")
        logger.info(f"Statistics: {self.stats}")
        
        return self.stats
    
    def _create_dataset_description(self):
        """Create BIDS dataset_description.json file."""
        desc_file = self.bids_root / "dataset_description.json"
        with open(desc_file, 'w') as f:
            json.dump(Config.DATASET_DESCRIPTION, f, indent=2)
        logger.info(f"Created {desc_file}")
    
    def _create_participants_file(self):
        """Create BIDS participants.tsv file."""
        participants_file = self.bids_root / "participants.tsv"
        
        # Collect all subject IDs
        subject_dirs = [d for d in self.bids_root.iterdir() if d.is_dir() and d.name.startswith('sub-')]
        
        with open(participants_file, 'w') as f:
            f.write("participant_id\n")
            for sub_dir in sorted(subject_dirs):
                f.write(f"{sub_dir.name}\n")
        
        logger.info(f"Created {participants_file} with {len(subject_dirs)} subjects")
    
    def _convert_subject(self, subject_dir: Path, dry_run: bool = False):
        """
        Convert a single subject's data.
        
        Parameters
        ----------
        subject_dir : Path
            Path to subject directory in Flywheel structure
        dry_run : bool
            If True, only report what would be done
        """
        subject_id = clean_subject_id(subject_dir.name)
        logger.info(f"Processing subject: {subject_id}")
        
        # Find sessions
        sessions_dir = subject_dir / "SESSIONS"
        if not sessions_dir.exists():
            logger.warning(f"No SESSIONS directory found for subject {subject_id}")
            return
        
        session_dirs = [d for d in sessions_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        
        for session_dir in session_dirs:
            try:
                self._convert_session(subject_id, session_dir, dry_run=dry_run)
                self.stats['sessions'] += 1
            except Exception as e:
                error_msg = f"Error processing session {session_dir.name} for subject {subject_id}: {str(e)}"
                logger.error(error_msg)
                self.stats['errors'].append(error_msg)
    
    def _convert_session(self, subject_id: str, session_dir: Path, dry_run: bool = False):
        """
        Convert a single session's data.
        
        Parameters
        ----------
        subject_id : str
            BIDS-formatted subject ID
        session_dir : Path
            Path to session directory
        dry_run : bool
            If True, only report what would be done
        """
        session_date, session_id = extract_session_info(session_dir.name)
        logger.info(f"  Processing session: {session_id} (date: {session_date})")
        
        # Find acquisitions
        acquisitions_dir = session_dir / "ACQUISITIONS"
        if not acquisitions_dir.exists():
            logger.warning(f"No ACQUISITIONS directory found for session {session_dir.name}")
            return
        
        acquisition_dirs = [d for d in acquisitions_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        
        for acq_dir in acquisition_dirs:
            try:
                self._convert_acquisition(subject_id, session_id, acq_dir, dry_run=dry_run)
                self.stats['acquisitions'] += 1
            except Exception as e:
                error_msg = f"Error processing acquisition {acq_dir.name}: {str(e)}"
                logger.error(error_msg)
                self.stats['errors'].append(error_msg)
    
    def _convert_acquisition(self, subject_id: str, session_id: str, acq_dir: Path, dry_run: bool = False):
        """
        Convert a single acquisition's data.
        
        Parameters
        ----------
        subject_id : str
            BIDS-formatted subject ID
        session_id : str
            BIDS-formatted session ID
        acq_dir : Path
            Path to acquisition directory
        dry_run : bool
            If True, only report what would be done
        """
        acquisition_name = acq_dir.name
        logger.info(f"    Processing acquisition: {acquisition_name}")
        
        # Determine modality and suffix
        modality, suffix, acq_label = determine_modality(acquisition_name)
        
        # Find FILES directory
        files_dir = acq_dir / "FILES"
        if not files_dir.exists():
            logger.warning(f"No FILES directory in {acq_dir.name}")
            return
        
        # Create BIDS output directory
        bids_subject_dir = self.bids_root / f"sub-{subject_id}"
        bids_session_dir = bids_subject_dir / f"ses-{session_id}"
        bids_modality_dir = bids_session_dir / modality
        
        if not dry_run:
            bids_modality_dir.mkdir(parents=True, exist_ok=True)
        
        # Build BIDS filename prefix
        prefix = f"sub-{subject_id}_ses-{session_id}"
        if acq_label:
            prefix += f"_{acq_label}"
        
        # Process based on modality
        if modality == 'dwi':
            self._convert_dwi_acquisition(files_dir, bids_modality_dir, prefix, suffix, dry_run)
        else:
            self._convert_anat_acquisition(files_dir, bids_modality_dir, prefix, suffix, dry_run)
    
    def _convert_dwi_acquisition(self, files_dir: Path, output_dir: Path, prefix: str, suffix: str, dry_run: bool):
        """Convert diffusion-weighted imaging acquisition."""
        dwi_files = find_diffusion_files(files_dir)
        
        # Copy NIfTI file
        if dwi_files['nifti']:
            output_nifti = output_dir / f"{prefix}_{suffix}.nii.gz"
            if dry_run:
                logger.info(f"      Would copy: {dwi_files['nifti'].name} -> {output_nifti}")
            else:
                shutil.copy2(dwi_files['nifti'], output_nifti)
                logger.info(f"      Copied: {output_nifti.name}")
                self.stats['files_copied'] += 1
            
            # Copy JSON sidecar if exists
            json_sidecar = self._find_json_sidecar(dwi_files['nifti'])
            if json_sidecar:
                output_json = output_dir / f"{prefix}_{suffix}.json"
                if dry_run:
                    logger.info(f"      Would copy: {json_sidecar.name} -> {output_json}")
                else:
                    shutil.copy2(json_sidecar, output_json)
                    logger.info(f"      Copied: {output_json.name}")
                    self.stats['files_copied'] += 1
        
        # Copy bval file
        if dwi_files['bval']:
            output_bval = output_dir / f"{prefix}_{suffix}.bval"
            if dry_run:
                logger.info(f"      Would copy: {dwi_files['bval'].name} -> {output_bval}")
            else:
                shutil.copy2(dwi_files['bval'], output_bval)
                logger.info(f"      Copied: {output_bval.name}")
                self.stats['files_copied'] += 1
        
        # Copy and rename bvec/mvec file
        if dwi_files['bvec']:
            output_bvec = output_dir / f"{prefix}_{suffix}.bvec"
            if dry_run:
                logger.info(f"      Would copy: {dwi_files['bvec'].name} -> {output_bvec}")
            else:
                shutil.copy2(dwi_files['bvec'], output_bvec)
                logger.info(f"      Copied: {output_bvec.name}")
                self.stats['files_copied'] += 1
    
    def _find_json_sidecar(self, nifti_file: Path) -> Optional[Path]:
        """
        Find JSON sidecar for a NIfTI file.
        
        Parameters
        ----------
        nifti_file : Path
            Path to NIfTI file
            
        Returns
        -------
        Optional[Path]
            Path to JSON sidecar if found, None otherwise
        """
        # Remove .nii.gz or .nii extension
        stem = nifti_file.name.replace('.nii.gz', '').replace('.nii', '')
        json_file = nifti_file.parent / f"{stem}.json"
        
        if json_file.exists():
            # Check if it's not a flywheel metadata file
            if not json_file.name.endswith('.flywheel.json'):
                return json_file
        
        return None
    
    def _convert_anat_acquisition(self, files_dir: Path, output_dir: Path, prefix: str, suffix: str, dry_run: bool):
        """Convert anatomical acquisition."""
        nifti_files = find_nifti_files(files_dir)
        
        for i, nifti_file in enumerate(nifti_files):
            # If multiple files, add run number
            if len(nifti_files) > 1:
                run_suffix = f"_run-{i+1:02d}"
            else:
                run_suffix = ""
            
            output_nifti = output_dir / f"{prefix}{run_suffix}_{suffix}.nii.gz"
            
            if dry_run:
                logger.info(f"      Would copy: {nifti_file.name} -> {output_nifti}")
            else:
                shutil.copy2(nifti_file, output_nifti)
                logger.info(f"      Copied: {output_nifti.name}")
                self.stats['files_copied'] += 1
            
            # Copy JSON sidecar if exists
            json_sidecar = self._find_json_sidecar(nifti_file)
            if json_sidecar:
                output_json = output_dir / f"{prefix}{run_suffix}_{suffix}.json"
                if dry_run:
                    logger.info(f"      Would copy: {json_sidecar.name} -> {output_json}")
                else:
                    shutil.copy2(json_sidecar, output_json)
                    logger.info(f"      Copied: {output_json.name}")
                    self.stats['files_copied'] += 1
