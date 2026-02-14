#!/usr/bin/env python3
"""Quick test of FOD loading."""

import numpy as np
import subprocess
import os
import tempfile
import nibabel as nib

# Test loading one subject
deriv = 'data/derivatives/mrtrix_preproc_flip1'
subj = 'sub-Apallial17NmF'
ses = 'ses-240816091633'
base = os.path.join(deriv, subj, ses)

# Load FOD
fod_path = os.path.join(base, 'dwi', 'wm_fod_norm.mif')
print(f'Loading FOD from: {fod_path}')

with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
    tmp_path = tmp.name

subprocess.run(['mrconvert', fod_path, tmp_path, '-force', '-quiet'], check=True)
img = nib.load(tmp_path)
fod_data = img.get_fdata()
print(f'FOD shape: {fod_data.shape}')
os.remove(tmp_path)

# Load graft mask
mask_path = os.path.join(base, 'manual_roi', 'fod_space', 'roi_Graft_left_flip_fod.mif')
print(f'Loading mask from: {mask_path}')

with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
    tmp_path = tmp.name

subprocess.run(['mrconvert', mask_path, tmp_path, '-force', '-quiet'], check=True)
img = nib.load(tmp_path)
mask_data = img.get_fdata() > 0.5
print(f'Mask shape: {mask_data.shape}')
print(f'Mask voxels: {mask_data.sum()}')
os.remove(tmp_path)

# Find A-P slice range
mask_slices = np.any(mask_data, axis=(1, 2))
slice_indices = np.where(mask_slices)[0]
print(f'Graft spans slices: {slice_indices.min()} to {slice_indices.max()} ({len(slice_indices)} slices)')

# Get sample FOD coefficients
slice_idx = slice_indices[len(slice_indices)//2]  # middle slice
slice_mask = mask_data[slice_idx, :, :]
y_coords, z_coords = np.where(slice_mask)
print(f'Slice {slice_idx} has {len(y_coords)} voxels')

# Check one voxel's SH coefficients
y, z = y_coords[0], z_coords[0]
sh_coeffs = fod_data[slice_idx, y, z, :]
print(f'SH coefficients range: {sh_coeffs.min():.4f} to {sh_coeffs.max():.4f}')
print(f'First 5 coefficients: {sh_coeffs[:5]}')
