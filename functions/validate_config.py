#!/usr/bin/env python3
"""
Configuration Validation Script for MBRAIN Pipeline.

Checks:
1.  Structure of config1.yaml (required keys).
2.  Existence of defined paths (base_dir, derivative_dir, atlas_dir).
3.  Accessibility of external tools (tckgen, tcksift2, matlab).
4.  Existence of input data for all defined subjects.
5.  Existence of Atlas ROI files (for CTRL).

Usage:
    python functions/validate_config.py
    
    # Or import in other scripts:
    from functions.validate_config import validate_pipeline_config
    if not validate_pipeline_config():
        sys.exit(1)
"""

import os
import sys
import shutil
from pathlib import Path

# Add parent dir to path to import config_loader
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config, get_paths, get_subject_groups

def check_executable(tool_name):
    """Check if a command line tool is available."""
    if shutil.which(tool_name):
        return True
    return False

def validate_pipeline_config(config_path=None):
    """
    Run comprehensive validation of the pipeline configuration and environment.
    
    Returns:
        bool: True if all checks pass, False otherwise.
    """
    print("="*60)
    print("MBRAIN Pipeline Configuration Validation")
    print("="*60)
    
    issues = []
    
    # 1. Load Config
    try:
        cfg = load_config(config_path)
        print("✅ Config loaded successfully.")
    except Exception as e:
        print(f"❌ Failed to load config: {e}")
        return False

    # 2. Check Paths
    paths = get_paths()
    required_paths = ['base_dir', 'derivatives_dir', 'atlas_dir']
    
    print("\nChecking Directories:")
    for p_name in required_paths:
        path_val = paths.get(p_name)
        if not path_val:
            issues.append(f"Missing config path: paths.{p_name}")
            print(f"  ❌ {p_name}: Not defined")
            continue
            
        p = Path(path_val)
        if p.exists() and p.is_dir():
            print(f"  ✅ {p_name}: {p}")
        else:
            issues.append(f"Directory not found: {p_name} ({p})")
            print(f"  ❌ {p_name}: NOT FOUND ({p})")

    # 3. Check Tools
    print("\nChecking External Tools:")
    required_tools = ['tckgen', 'tcksift2', 'mrstats', 'mrview']
    for tool in required_tools:
        if check_executable(tool):
            print(f"  ✅ {tool}: Found")
        else:
            issues.append(f"Missing executable: {tool}")
            print(f"  ❌ {tool}: NOT FOUND in PATH")
            
    # Check MATLAB
    matlab_bin = paths.get('matlab_bin')
    if matlab_bin and os.path.exists(matlab_bin):
         print(f"  ✅ matlab: Found at {matlab_bin}")
    elif shutil.which('matlab'):
         print(f"  ✅ matlab: Found in PATH")
    else:
         issues.append("MATLAB executable not found")
         print(f"  ❌ matlab: NOT FOUND")

    # 4. Check Subjects
    print("\nChecking Subjects:")
    groups = get_subject_groups()
    base_path = Path(paths['derivatives_dir'])
    
    total_subjects = 0
    missing_subjects = 0
    
    for group, subject_list in groups.items():
        print(f"  Group {group}:")
        for subj_entry in subject_list:
            if isinstance(subj_entry, dict):
                sid = subj_entry['id']
                sess = subj_entry.get('session')
            else:
                sid = subj_entry
                sess = None
            
            total_subjects += 1
            
            # Check subject dir
            s_path = base_path / sid
            if not s_path.exists():
                issues.append(f"Missing subject data: {sid}")
                print(f"    ❌ {sid}: Directory not found")
                missing_subjects += 1
                continue
                
            # Check session if specified
            if sess:
                sess_path = s_path / sess
                if not sess_path.exists():
                    issues.append(f"Missing session data: {sid}/{sess}")
                    print(f"    ❌ {sid}/{sess}: Session not found")
                    missing_subjects += 1
            
    if missing_subjects == 0:
        print(f"  ✅ All {total_subjects} subjects found.")
    else:
        print(f"  ❌ {missing_subjects}/{total_subjects} subjects missing.")

    # 5. Check Atlas ROIs (for CTRL)
    print("\nChecking Atlas ROIs:")
    atlas_roi_dir = Path(paths.get('atlas_rois_dir', ''))
    if atlas_roi_dir.exists():
        # Check a few key ones (lowercase as per file system)
        key_rois = ['thalamus', 'striatum', 'isocortex']
        found = 0
        for roi in key_rois:
            matches = list(atlas_roi_dir.glob(f"*{roi}*.nii.gz"))
            if matches:
                found += 1
        
        if found == len(key_rois):
            print(f"  ✅ Atlas ROIs found in {atlas_roi_dir}")
        else:
             issues.append("Some Atlas ROI files missing")
             print(f"  ❌ Missing some Atlas ROIs in {atlas_roi_dir}")
             print(f"     Found {found}/{len(key_rois)} required ROIs.")
    else:
         # atlas_rois_dir might not be strictly required if using manual ROIs, but good to check
         print(f"  ⚠️ Atlas ROI dir not found: {atlas_roi_dir}")

    print("\n" + "="*60)
    if issues:
        print(f"VALIDATION FAILED with {len(issues)} issues:")
        for i in issues:
            print(f"  - {i}")
        return False
    else:
        print("VALIDATION SUCCESSFUL")
        return True

if __name__ == "__main__":
    success = validate_pipeline_config()
    sys.exit(0 if success else 1)
