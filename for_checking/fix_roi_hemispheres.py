#!/usr/bin/env python3
"""Re-extract ROIs with corrected left/right hemisphere assignment."""

import nibabel as nib
import numpy as np
from pathlib import Path
import json

# Paths
deriv = Path('/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2')
atlas_dir = Path('/Users/glerma/soft/mbrain/local/atlas/allen_ccf')
subject = 'sub-129-921-R1F'
session = 'ses-260115121231'

base = deriv / subject / session / 'dwi' / 'atlas'
output_dir = base / 'rois'
annotation_file = base / 'atlas_annotation_warped.nii.gz'
structures_file = atlas_dir / 'structures.json'

# Load structures
with open(structures_file) as f:
    data = json.load(f)
structures = data['msg'] if 'msg' in data else data

def get_all_child_ids(structures, parent_id):
    """Recursively get all child IDs for a parent structure."""
    all_ids = set()
    def recurse(nodes, pid):
        for node in nodes:
            if node['id'] == pid:
                all_ids.add(node['id'])
                if 'children' in node:
                    for child in node['children']:
                        recurse([child], child['id'])
            elif 'children' in node:
                recurse(node['children'], pid)
    recurse(structures, parent_id)
    if not all_ids:
        all_ids.add(parent_id)
    return list(all_ids)

# ROI definitions
roi_defs = {
    'isocortex': [315],
    'piriform_complex': [961, 909, 788, 566],
    'striatum': [485, 493],
    'thalamus': [549],
    'cerebellum': [512],
    'bla': [295],
    'aca': [31]
}

# Load annotation
img = nib.load(annotation_file)
ann_data = img.get_fdata().astype(np.int32)
affine = img.affine

# FLIP2 ORIENTATION:
# Axis 2 (Z) = medial-lateral (left-right)
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
midline_z = ann_data.shape[2] // 2
print(f'Image shape: {ann_data.shape}, midline_z = {midline_z}')

for roi_name, parent_ids in roi_defs.items():
    all_ids = []
    for pid in parent_ids:
        child_ids = get_all_child_ids(structures, pid)
        all_ids.extend(child_ids)
    
    mask = np.isin(ann_data, all_ids).astype(np.uint8)
    
    if mask.sum() == 0:
        print(f'  Warning: No voxels for {roi_name}')
        continue
    
    # Save full ROI
    nib.save(nib.Nifti1Image(mask, affine), output_dir / f'roi_{roi_name}.nii.gz')
    
    # LEFT = HIGH voxel Z indices (Z >= midline) → LOW world Z
    mask_left = mask.copy()
    mask_left[:, :, :midline_z] = 0  # Zero out low voxel Z (which is high world Z)
    
    # RIGHT = LOW voxel Z indices (Z < midline) → HIGH world Z
    mask_right = mask.copy()
    mask_right[:, :, midline_z:] = 0  # Zero out high voxel Z (which is low world Z)
    
    nib.save(nib.Nifti1Image(mask_left, affine), output_dir / f'roi_{roi_name}_left.nii.gz')
    nib.save(nib.Nifti1Image(mask_right, affine), output_dir / f'roi_{roi_name}_right.nii.gz')
    
    print(f'  {roi_name}: {mask.sum()} voxels (L:{mask_left.sum()}, R:{mask_right.sum()})')

print('Done!')
