#!/usr/bin/env python3
"""
Generic atlas coregistration script.

Coregisters Allen CCF atlas to subject T2w images using ANTs.
Supports manual shift corrections defined in config1.yaml.

ROI extraction uses Allen CCF structure hierarchy to get all child IDs
for each parent structure ID defined in config1.yaml. This is necessary
because the annotation volume contains leaf node IDs, not parent IDs.

Usage:
    python atlas_coregistration.py                          # All subjects in all groups
    python atlas_coregistration.py --group CTRL             # All subjects in CTRL group
    python atlas_coregistration.py --subject sub-129-921-R1F  # Specific subject
    python atlas_coregistration.py --dry-run                # Preview commands
"""

import argparse
import subprocess
import os
import json
from pathlib import Path
import numpy as np
import nibabel as nib

from config_loader import (
    load_config, get_paths, get_subject_groups,
    get_roi_definitions, get_manual_shift
)


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session) tuples based on filters."""
    groups = get_subject_groups()
    subjects = []
    
    if subject:
        # Find specific subject in any group
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session']))
                    return subjects
        raise ValueError(f"Subject {subject} not found in any group")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found. Available: {list(groups.keys())}")
        for s in groups[group]:
            subjects.append((s['id'], s['session']))
    else:
        # All subjects
        for g, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session']))
    
    return subjects


def apply_manual_shift(nifti_path, shift_mm, output_path):
    """Apply manual shift to a NIfTI file by modifying the affine."""
    img = nib.load(nifti_path)
    new_affine = img.affine.copy()
    
    # shift_mm is [x, y, z] in mm
    new_affine[0, 3] += shift_mm[0]
    new_affine[1, 3] += shift_mm[1]
    new_affine[2, 3] += shift_mm[2]
    
    new_img = nib.Nifti1Image(img.get_fdata(), new_affine, img.header)
    nib.save(new_img, output_path)
    return output_path


def run_atlas_coregistration(subject_id, session, dry_run=False):
    """Run atlas coregistration for a single subject."""
    cfg = load_config()
    paths = get_paths()
    
    # Paths
    deriv_dir = Path(paths['derivatives_dir'])
    atlas_dir = Path(paths['atlas_dir'])
    
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    atlas_out_dir = dwi_dir / "atlas"
    rois_out_dir = atlas_out_dir / "rois"
    
    # Input files
    t2w_file = dwi_dir / "t2w_highres_flip.mif"
    template_file = atlas_dir / "average_template_25.nii.gz"
    annotation_file = atlas_dir / "annotation_25.nii.gz"
    
    print(f"\n{'='*60}")
    print(f"Atlas Coregistration: {subject_id}/{session}")
    print(f"{'='*60}")
    
    # Check inputs
    if not t2w_file.exists():
        print(f"  ERROR: T2w file not found: {t2w_file}")
        return False
    
    if dry_run:
        print(f"  [DRY RUN] Would process {subject_id}")
        return True
    
    # Create output directories
    atlas_out_dir.mkdir(parents=True, exist_ok=True)
    rois_out_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert T2w to NIfTI for ANTs
    t2w_nii = atlas_out_dir / "t2w_for_atlas.nii.gz"
    subprocess.run(["mrconvert", str(t2w_file), str(t2w_nii), "-force"], check=True)
    
    # Get manual shift if defined
    shift = get_manual_shift(subject_id)
    if any(s != 0 for s in shift):
        print(f"  Applying manual shift: {shift} mm")
        template_to_use = atlas_out_dir / f"allen_template_shifted_{shift[1]}mm.nii.gz"
        apply_manual_shift(template_file, shift, template_to_use)
    else:
        template_to_use = template_file
    
    # Run ANTs registration
    output_prefix = atlas_out_dir / "atlas_coreg_"
    
    ants_params = cfg.get('ants', {})
    transform_type = ants_params.get('transform_type', 's')
    
    cmd = [
        "antsRegistrationSyN.sh",
        "-d", "3",
        "-t", transform_type,
        "-f", str(t2w_nii),
        "-m", str(template_to_use),
        "-o", str(output_prefix)
    ]
    
    print(f"  Running ANTs registration...")
    print(f"    {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  ERROR: ANTs failed: {result.stderr[:500]}")
        return False
    
    # Apply warp to annotation (atlas labels)
    warped_annotation = atlas_out_dir / "atlas_annotation_warped.nii.gz"
    warp_file = atlas_out_dir / "atlas_coreg_1Warp.nii.gz"
    affine_file = atlas_out_dir / "atlas_coreg_0GenericAffine.mat"
    
    cmd = [
        "antsApplyTransforms",
        "-d", "3",
        "-i", str(annotation_file),
        "-r", str(t2w_nii),
        "-o", str(warped_annotation),
        "-t", str(warp_file),
        "-t", str(affine_file),
        "-n", "NearestNeighbor"
    ]
    
    print(f"  Applying warp to annotation...")
    subprocess.run(cmd, check=True)
    
    # Also warp template for QC
    warped_template = atlas_out_dir / "atlas_template_warped.nii.gz"
    cmd = [
        "antsApplyTransforms",
        "-d", "3",
        "-i", str(template_to_use),
        "-r", str(t2w_nii),
        "-o", str(warped_template),
        "-t", str(warp_file),
        "-t", str(affine_file),
        "-n", "Linear"
    ]
    subprocess.run(cmd, check=True)
    
    # Extract ROIs
    print(f"  Extracting ROIs...")
    paths = get_paths()
    atlas_dir = Path(paths['atlas_dir'])
    structures_file = atlas_dir / "structures.json"
    extract_rois(warped_annotation, rois_out_dir, subject_id, structures_file)
    
    print(f"  ✓ Atlas coregistration complete")
    return True


def load_allen_structures(structures_file):
    """Load Allen CCF structure hierarchy from JSON."""
    import json
    with open(structures_file) as f:
        data = json.load(f)
    return data.get('msg', [])


def get_all_child_ids(structures, parent_id):
    """Get all descendant structure IDs from Allen CCF hierarchy.
    
    The annotation volume contains LEAF node IDs, not parent IDs.
    This function recursively collects all child IDs under a parent.
    """
    def find_structure(structs, target_id):
        for s in structs:
            if s['id'] == target_id:
                return s
            if 'children' in s:
                result = find_structure(s['children'], target_id)
                if result:
                    return result
        return None
    
    def collect_ids(struct):
        ids = [struct['id']]
        if 'children' in struct:
            for child in struct['children']:
                ids.extend(collect_ids(child))
        return ids
    
    struct = find_structure(structures, parent_id)
    if struct:
        return collect_ids(struct)
    return [parent_id]  # Return parent ID if not found in hierarchy


def extract_rois(annotation_file, output_dir, subject_id, structures_file):
    """Extract ROIs from warped annotation based on config definitions.
    
    Uses Allen CCF hierarchy to get all child structure IDs for each parent ID
    defined in config1.yaml. This is necessary because the annotation volume
    contains leaf node IDs, not parent structure IDs.
    """
    roi_defs = get_roi_definitions()
    structures = load_allen_structures(structures_file)
    
    img = nib.load(annotation_file)
    data = img.get_fdata().astype(np.int32)
    affine = img.affine
    
    # FLIP1 ORIENTATION (after axes 2,1,0 permutation + flip 1):
    # - Axis 0 (X): Anterior-Posterior (A→P)
    # - Axis 1 (Y): Dorsal-Ventral (D→V)  
    # - Axis 2 (Z): Medial-Lateral (Left-Right)
    #
    # IMPORTANT: The affine has NEGATIVE Z direction (affine[2,2] < 0)!
    # This means: High voxel Z index → LOW world Z coordinate
    #
    # Based on validated AP manual ROIs (in WORLD coordinates):
    #   Graft_left world Z = 73.9 (LOW world Z → HIGH voxel Z index)
    #   Graft_right world Z = 107.7 (HIGH world Z → LOW voxel Z index)
    #
    # Therefore, in VOXEL INDEX terms:
    #   LEFT = HIGH voxel Z indices (Z >= midline)
    #   RIGHT = LOW voxel Z indices (Z < midline)
    midline_z = data.shape[2] // 2
    
    for roi_name, roi_cfg in roi_defs.items():
        parent_ids = roi_cfg['allen_ids']
        
        # Get ALL child IDs for each parent ID
        all_ids = []
        for pid in parent_ids:
            child_ids = get_all_child_ids(structures, pid)
            all_ids.extend(child_ids)
        
        # Create binary mask using np.isin for efficient matching
        mask = np.isin(data, all_ids).astype(np.uint8)
        
        if mask.sum() == 0:
            print(f"    Warning: No voxels found for {roi_name} (checked {len(all_ids)} IDs)")
            continue
        
        # Save full ROI
        roi_img = nib.Nifti1Image(mask, affine)
        nib.save(roi_img, output_dir / f"roi_{roi_name}.nii.gz")
        
        # Split into hemispheres along axis 2 (Z = medial-lateral)
        # LEFT = HIGH voxel Z indices (Z >= midline) → LOW world Z
        mask_left = mask.copy()
        mask_left[:, :, :midline_z] = 0  # Zero out low voxel Z
        
        # RIGHT = LOW voxel Z indices (Z < midline) → HIGH world Z
        mask_right = mask.copy()
        mask_right[:, :, midline_z:] = 0  # Zero out high voxel Z
        
        nib.save(nib.Nifti1Image(mask_left, affine), 
                 output_dir / f"roi_{roi_name}_left.nii.gz")
        nib.save(nib.Nifti1Image(mask_right, affine), 
                 output_dir / f"roi_{roi_name}_right.nii.gz")
        
        print(f"    {roi_name}: {mask.sum()} voxels (L:{mask_left.sum()}, R:{mask_right.sum()})")


def main():
    parser = argparse.ArgumentParser(
        description='Atlas coregistration for mouse brain MRI',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group (CTRL, AP, APTR)')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    
    print(f"Processing {len(subjects)} subjects...")
    
    results = []
    for subject_id, session in subjects:
        success = run_atlas_coregistration(subject_id, session, dry_run=args.dry_run)
        results.append((subject_id, success))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for subj, success in results:
        status = "✓" if success else "✗"
        print(f"  {status} {subj}")


if __name__ == "__main__":
    main()
