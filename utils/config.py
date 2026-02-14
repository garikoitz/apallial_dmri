"""
Configuration settings for the mbrain project.

This module provides backward-compatible Config class that now loads
settings from config1.yaml via config_loader.
"""

import os
import sys
from pathlib import Path

# Add parent directory to path to import config_loader
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from config_loader import load_config, get_paths
    _cfg = load_config()
    _paths = get_paths()
    _USE_CONFIG = True
except ImportError:
    _USE_CONFIG = False


class Config:
    """Configuration class for mbrain project."""
    
    # Project paths - load from config1.yaml if available
    PROJECT_ROOT = Path(__file__).parent.parent
    
    if _USE_CONFIG:
        DATA_DIR = Path(_paths.get('base_dir', PROJECT_ROOT)) / "local"
        FLYWHEEL_DIR = DATA_DIR / "flywheel" / "fMRI"
        BIDS_DIR = Path(_paths.get('bids_dir', PROJECT_ROOT / "data" / "bids"))
        DERIVATIVES_DIR = Path(_paths.get('derivatives_dir', PROJECT_ROOT / "data" / "derivatives")).parent
    else:
        DATA_DIR = PROJECT_ROOT / "local"
        FLYWHEEL_DIR = DATA_DIR / "flywheel" / "fMRI"
        BIDS_DIR = PROJECT_ROOT / "data" / "bids"
        DERIVATIVES_DIR = PROJECT_ROOT / "data" / "derivatives"
    
    # BIDS dataset information
    DATASET_NAME = "mouse_brain_fmri"
    DATASET_DESCRIPTION = {
        "Name": "Mouse Brain fMRI and DTI Dataset",
        "BIDSVersion": "1.9.0",
        "DatasetType": "raw",
        "Authors": [""],
        "License": "",
        "Acknowledgements": "",
        "HowToAcknowledge": "",
        "Funding": [""],
        "ReferencesAndLinks": [""],
        "DatasetDOI": ""
    }
    
    # Modality mappings from Flywheel acquisition names to BIDS modalities
    MODALITY_MAPPING = {
        "T2w": "anat",
        "T2star": "anat",
        "DTI": "dwi",
        "Localizer": "anat"
    }
    
    # Suffix mappings
    SUFFIX_MAPPING = {
        "T2w_lowres": "T2w",
        "T2w_highres": "T2w",
        "T2star_FID_EPI": "T2starw",
        "DTI_EPI": "dwi",
        "Localizer": "T2w"
    }
    
    @classmethod
    def ensure_directories(cls):
        """Create necessary directories if they don't exist."""
        cls.BIDS_DIR.mkdir(parents=True, exist_ok=True)
        cls.DERIVATIVES_DIR.mkdir(parents=True, exist_ok=True)
        
    @classmethod
    def get_subject_dir(cls, subject_id):
        """Get BIDS subject directory path."""
        return cls.BIDS_DIR / f"sub-{subject_id}"
    
    @classmethod
    def get_session_dir(cls, subject_id, session_id):
        """Get BIDS session directory path."""
        return cls.get_subject_dir(subject_id) / f"ses-{session_id}"
