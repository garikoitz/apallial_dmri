"""
File utility functions for finding and processing data files.
"""

import os
from pathlib import Path
import re
from typing import List, Dict, Optional, Tuple


def find_nifti_files(directory: Path) -> List[Path]:
    """
    Find all NIfTI files in a directory.
    
    Parameters
    ----------
    directory : Path
        Directory to search for NIfTI files
        
    Returns
    -------
    List[Path]
        List of paths to NIfTI files
    """
    nifti_files = []
    for ext in ['*.nii.gz', '*.nii']:
        nifti_files.extend(directory.glob(ext))
    return sorted(nifti_files)


def find_diffusion_files(directory: Path) -> Dict[str, Optional[Path]]:
    """
    Find diffusion-related files (nifti, bval, bvec) in a directory.
    
    Parameters
    ----------
    directory : Path
        Directory to search for diffusion files
        
    Returns
    -------
    Dict[str, Optional[Path]]
        Dictionary with keys 'nifti', 'bval', 'bvec' pointing to file paths
    """
    files = {
        'nifti': None,
        'bval': None,
        'bvec': None
    }
    
    # Find NIfTI file
    nifti_files = find_nifti_files(directory)
    if nifti_files:
        # Prefer files with gradient info (mvec/bval)
        for nf in nifti_files:
            stem = nf.name.replace('.nii.gz', '').replace('.nii', '')
            if (directory / f"{stem}.bval").exists() or (directory / f"{stem}.mvec").exists():
                files['nifti'] = nf
                break
        if files['nifti'] is None:
            files['nifti'] = nifti_files[0]
    
    # Find bval file
    bval_files = list(directory.glob('*.bval'))
    if bval_files:
        files['bval'] = bval_files[0]
    
    # Find bvec/mvec file (mvec is Bruker's version of bvec)
    bvec_files = list(directory.glob('*.bvec'))
    mvec_files = list(directory.glob('*.mvec'))
    
    if bvec_files:
        files['bvec'] = bvec_files[0]
    elif mvec_files:
        files['bvec'] = mvec_files[0]
    
    return files


def clean_subject_id(subject_id: str) -> str:
    """
    Clean subject ID to make it BIDS-compatible.
    
    Parameters
    ----------
    subject_id : str
        Raw subject ID from Flywheel
        
    Returns
    -------
    str
        Cleaned subject ID (alphanumeric and hyphens only)
    """
    # Remove spaces and special characters except hyphens
    cleaned = re.sub(r'[^a-zA-Z0-9-]', '', subject_id)
    return cleaned


def extract_session_info(session_name: str) -> Tuple[str, str]:
    """
    Extract session date and ID from Flywheel session name.
    
    Parameters
    ----------
    session_name : str
        Session name from Flywheel (e.g., "Pasca_596-C-OO-NmM_9-11-24 - 240911103349")
        
    Returns
    -------
    Tuple[str, str]
        (session_date, session_id) where session_date is like "20240911"
    """
    # Try to extract the numeric session ID (e.g., "240911103349")
    match = re.search(r'(\d{12,})', session_name)
    if match:
        session_id = match.group(1)
        # Extract date portion (first 6 digits are YYMMDD)
        if len(session_id) >= 6:
            date_part = session_id[:6]
            session_date = f"20{date_part[:2]}{date_part[2:4]}{date_part[4:6]}"
        else:
            session_date = "unknown"
    else:
        session_id = "01"
        session_date = "unknown"
    
    return session_date, session_id


def determine_modality(acquisition_name: str) -> Tuple[str, str]:
    """
    Determine BIDS modality and suffix from acquisition name.
    
    Parameters
    ----------
    acquisition_name : str
        Acquisition name from Flywheel
        
    Returns
    -------
    Tuple[str, str]
        (modality, suffix) e.g., ("anat", "T2w") or ("dwi", "dwi")
    """
    name_lower = acquisition_name.lower()
    
    if 't2w' in name_lower:
        if 'highres' in name_lower:
            return 'anat', 'T2w', 'acq-highres'
        elif 'lowres' in name_lower:
            return 'anat', 'T2w', 'acq-lowres'
        else:
            return 'anat', 'T2w', ''
    elif 't2star' in name_lower:
        return 'anat', 'T2starw', ''
    elif 'dti' in name_lower or 'dwi' in name_lower:
        return 'dwi', 'dwi', ''
    elif 'localizer' in name_lower:
        return 'anat', 'T2w', 'acq-localizer'
    else:
        return 'anat', 'T2w', ''
