#!/usr/bin/env python3
"""
Atlas ROI Extraction for CTRL subjects (flip1 pipeline)

This script registers the Allen CCF atlas to a subject's T2w image and extracts
ROIs for tractography targets. Used for CTRL subjects (no graft = use atlas ROIs).

Usage:
    python functions/atlas_roi_extraction.py --subject sub-129-921-R1F

Requirements:
    - ANTs (antsRegistrationSyN.sh, antsApplyTransforms)
    - nibabel, numpy
    - Subject must have t2w_highres_flip.mif in dwi/ folder

Output:
    - atlas/annotation_warped.nii.gz (warped Allen annotation)
    - atlas/rois/roi_{region}.nii.gz (bilateral ROIs)
    - atlas/rois/roi_{region}_left.nii.gz (left hemisphere)
    - atlas/rois/roi_{region}_right.nii.gz (right hemisphere)

Author: MBRAIN project
Date: 2026-01-29
"""

import argparse
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
import numpy as np
import nibabel as nib
import json
import yaml

def load_config():
    """Load config1.yaml"""
    config_path = Path(__file__).parent.parent / "config1.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)

def find_structure(structs, target_id):
    """Recursively find structure by ID in Allen hierarchy"""
    for s in structs:
        if s['id'] == target_id:
            return s
        if 'children' in s:
            result = find_structure(s['children'], target_id)
            if result:
                return result
    return None

def collect_ids(struct):
    """Collect all child IDs recursively"""
    ids = [struct['id']]
    if 'children' in struct:
        for child in struct['children']:
            ids.extend(collect_ids(child))
    return ids

def get_all_child_ids(structures, parent_id):
    """Get all descendant IDs for a parent structure"""
    struct = find_structure(structures, parent_id)
    if struct:
        return collect_ids(struct)
    return [parent_id]

def run_atlas_registration(t2w_path, template_path, output_dir):
    """Run ANTs SyN registration of atlas template to subject T2w"""
    print(f"Running ANTs registration (this takes ~5 minutes)...")
    
    # Convert T2w to nifti for ANTs
    t2w_nii = output_dir / "t2w.nii.gz"
    subprocess.run([
        "mrconvert", str(t2w_path), str(t2w_nii), "-force"
    ], check=True)
    
    # Copy template
    template_nii = output_dir / "template.nii.gz"
    shutil.copy(template_path, template_nii)
    
    # Run ANTs registration
    subprocess.run([
        "antsRegistrationSyN.sh",
        "-d", "3",
        "-t", "s",  # SyN
        "-f", str(t2w_nii),
        "-m", str(template_nii),
        "-o", str(output_dir / "atlas_coreg_")
    ], check=True)
    
    return output_dir / "atlas_coreg_1Warp.nii.gz", output_dir / "atlas_coreg_0GenericAffine.mat"

def apply_warp_to_annotation(annotation_path, t2w_path, warp_path, affine_path, output_dir):
    """Apply warp to annotation volume using nearest neighbor interpolation"""
    print("Applying warp to annotation...")
    
    t2w_nii = output_dir / "t2w.nii.gz"
    annotation_nii = output_dir / "annotation.nii.gz"
    shutil.copy(annotation_path, annotation_nii)
    
    warped_annotation = output_dir / "annotation_warped.nii.gz"
    
    subprocess.run([
        "antsApplyTransforms",
        "-d", "3",
        "-i", str(annotation_nii),
        "-r", str(t2w_nii),
        "-o", str(warped_annotation),
        "-n", "NearestNeighbor",
        "-t", str(warp_path),
        "-t", str(affine_path)
    ], check=True)
    
    return warped_annotation

def extract_rois(annotation_path, structures, roi_defs, output_dir):
    """Extract ROIs from warped annotation"""
    print("Extracting ROIs...")
    
    img = nib.load(annotation_path)
    data = img.get_fdata().astype(np.int32)
    affine = img.affine
    
    # FLIP1 ORIENTATION: midline along axis 2
    # CRITICAL: In flip1 space, anatomical LEFT = HIGH Z indices
    midline_z = data.shape[2] // 2
    
    rois_dir = output_dir / "rois"
    rois_dir.mkdir(exist_ok=True)
    
    results = {}
    
    for roi_name, roi_cfg in roi_defs.items():
        parent_ids = roi_cfg['allen_ids']
        
        # Get ALL child IDs
        all_ids = []
        for pid in parent_ids:
            child_ids = get_all_child_ids(structures, pid)
            all_ids.extend(child_ids)
        
        # Create binary mask
        mask = np.isin(data, all_ids).astype(np.uint8)
        
        if mask.sum() == 0:
            print(f"  Warning: No voxels found for {roi_name}")
            continue
        
        # Save bilateral ROI
        nib.save(nib.Nifti1Image(mask, affine), rois_dir / f"roi_{roi_name}.nii.gz")
        
        # Split into hemispheres
        # FLIP1 CONVENTION: Anatomical LEFT = HIGH Z, Anatomical RIGHT = LOW Z
        mask_left = mask.copy()
        mask_left[:, :, :midline_z] = 0  # Keep HIGH Z (anatomical left)
        
        mask_right = mask.copy()
        mask_right[:, :, midline_z:] = 0  # Keep LOW Z (anatomical right)
        
        nib.save(nib.Nifti1Image(mask_left, affine), rois_dir / f"roi_{roi_name}_left.nii.gz")
        nib.save(nib.Nifti1Image(mask_right, affine), rois_dir / f"roi_{roi_name}_right.nii.gz")
        
        results[roi_name] = {
            'total': int(mask.sum()),
            'left': int(mask_left.sum()),
            'right': int(mask_right.sum())
        }
        print(f"  {roi_name}: {mask.sum()} voxels (L:{mask_left.sum()}, R:{mask_right.sum()})")
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract Allen atlas ROIs for CTRL subject")
    parser.add_argument("--subject", required=True, help="Subject ID (e.g., sub-129-921-R1F)")
    parser.add_argument("--session", help="Session ID (auto-detected if not provided)")
    args = parser.parse_args()
    
    # Load config
    config = load_config()
    
    # Paths from config
    deriv_base = Path(config['paths']['derivatives_flip1'])
    atlas_dir = Path(config['paths']['local']) / "atlas" / "allen_ccf"
    
    # Find session if not provided
    subject_dir = deriv_base / args.subject
    if args.session:
        session = args.session
    else:
        sessions = list(subject_dir.glob("ses-*"))
        if not sessions:
            print(f"Error: No sessions found for {args.subject}")
            sys.exit(1)
        session = sessions[0].name
    
    session_dir = subject_dir / session
    t2w_path = session_dir / "dwi" / "t2w_highres_flip.mif"
    
    if not t2w_path.exists():
        print(f"Error: T2w not found: {t2w_path}")
        sys.exit(1)
    
    # Atlas files
    template_path = atlas_dir / "average_template_25.nii.gz"
    annotation_path = atlas_dir / "annotation_25.nii.gz"
    structures_path = atlas_dir / "structures.json"
    
    # Load Allen structures
    with open(structures_path) as f:
        structures_data = json.load(f)
    structures = structures_data.get('msg', [])
    
    # ROI definitions from config
    roi_source = config['atlas'].get('use_roi_definitions_from', 'rois')
    roi_defs = config[roi_source]
    
    # Create output directory
    atlas_output = session_dir / "atlas"
    atlas_output.mkdir(exist_ok=True)
    
    # Use temp directory for processing (avoids Google Drive path issues)
    with tempfile.TemporaryDirectory(prefix="mbrain_atlas_") as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Run registration
        warp_path, affine_path = run_atlas_registration(t2w_path, template_path, tmpdir)
        
        # Apply warp to annotation
        warped_annotation = apply_warp_to_annotation(
            annotation_path, t2w_path, warp_path, affine_path, tmpdir
        )
        
        # Extract ROIs
        results = extract_rois(warped_annotation, structures, roi_defs, tmpdir)
        
        # Copy results to subject folder
        print(f"\nCopying results to {atlas_output}...")
        shutil.copy(tmpdir / "annotation_warped.nii.gz", atlas_output)
        shutil.copytree(tmpdir / "rois", atlas_output / "rois", dirs_exist_ok=True)
        
        # Copy transform files
        for f in tmpdir.glob("atlas_coreg_*"):
            shutil.copy(f, atlas_output)
    
    print(f"\n✓ Atlas ROI extraction complete for {args.subject}")
    print(f"  Output: {atlas_output}")
    
    return results

if __name__ == "__main__":
    main()
