#!/usr/bin/env python3
"""Check which Allen CCF IDs are present in warped annotation and create ROIs."""

import nibabel as nib
import numpy as np
import json
import os
from pathlib import Path

# Paths
DERIV = Path("/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2")
ATLAS_DIR = Path("/Users/glerma/soft/mbrain/local/atlas/allen_ccf")
SUBJECT = "sub-129-921-R1F"
SESSION = "ses-260115121231"

# Load warped annotation
annot_path = DERIV / SUBJECT / SESSION / "dwi/atlas/atlas_annotation_warped.nii.gz"
print(f"Loading: {annot_path}")
img = nib.load(annot_path)
data = img.get_fdata().astype(int)
present_ids = set(np.unique(data))
print(f"Total unique IDs in warped annotation: {len(present_ids)}")

# Load structures
with open(ATLAS_DIR / "structures.json") as f:
    structures_data = json.load(f)
structures = structures_data.get('msg', [])

def find_structure(structs, target_id):
    for s in structs:
        if s['id'] == target_id:
            return s
        if 'children' in s:
            result = find_structure(s['children'], target_id)
            if result:
                return result
    return None

def get_all_ids(struct):
    ids = [struct['id']]
    if 'children' in struct:
        for child in struct['children']:
            ids.extend(get_all_ids(child))
    return ids

# Define ROIs with their parent Allen IDs
ROIS = {
    'isocortex': {'parent_id': 315, 'name': 'Isocortex'},
    'piriform_complex': {'parent_ids': [961, 909, 788, 566], 'name': 'Piriform_Complex'},  # PIR, ENT, PAA, TR
    'striatum': {'parent_ids': [485, 493], 'name': 'Striatum'},  # dorsal + ventral
    'thalamus': {'parent_id': 549, 'name': 'Thalamus'},
    'cerebellum': {'parent_id': 512, 'name': 'Cerebellum'},
}

# Output directory
out_dir = DERIV / SUBJECT / SESSION / "dwi/atlas/rois"
out_dir.mkdir(parents=True, exist_ok=True)

print("\n" + "=" * 80)
print("Creating ROI masks from Allen CCF annotation")
print("=" * 80)

for roi_key, roi_def in ROIS.items():
    print(f"\n{roi_def['name']}:")
    
    # Get all structure IDs (handle single or multiple parent IDs)
    all_ids = []
    if 'parent_id' in roi_def:
        struct = find_structure(structures, roi_def['parent_id'])
        if struct:
            all_ids = get_all_ids(struct)
    elif 'parent_ids' in roi_def:
        for pid in roi_def['parent_ids']:
            struct = find_structure(structures, pid)
            if struct:
                all_ids.extend(get_all_ids(struct))
    
    print(f"  Total structure IDs: {len(all_ids)}")
    
    # Check which are present
    present = [i for i in all_ids if i in present_ids]
    print(f"  IDs present in annotation: {len(present)}")
    
    # Create mask
    mask = np.isin(data, all_ids).astype(np.uint8)
    n_voxels = np.sum(mask > 0)
    print(f"  Total voxels: {n_voxels:,}")
    
    if n_voxels > 0:
        # Split into hemispheres
        x_center = mask.shape[0] // 2
        
        left_mask = np.zeros_like(mask)
        right_mask = np.zeros_like(mask)
        
        left_mask[:x_center, :, :] = mask[:x_center, :, :]
        right_mask[x_center:, :, :] = mask[x_center:, :, :]
        
        left_voxels = np.sum(left_mask > 0)
        right_voxels = np.sum(right_mask > 0)
        print(f"    Left hemisphere: {left_voxels:,} voxels")
        print(f"    Right hemisphere: {right_voxels:,} voxels")
        
        # Save masks
        roi_name = roi_def['name'].lower()
        
        left_file = out_dir / f"roi_{roi_name}_left.nii.gz"
        right_file = out_dir / f"roi_{roi_name}_right.nii.gz"
        
        left_img = nib.Nifti1Image(left_mask, img.affine, img.header)
        right_img = nib.Nifti1Image(right_mask, img.affine, img.header)
        
        nib.save(left_img, left_file)
        nib.save(right_img, right_file)
        
        print(f"  ✓ Saved: {left_file.name}")
        print(f"  ✓ Saved: {right_file.name}")
    else:
        print(f"  ⚠ WARNING: No voxels found for {roi_def['name']}")

print("\n" + "=" * 80)
print("Done! ROIs saved to:", out_dir)
print("=" * 80)
