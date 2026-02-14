#!/usr/bin/env python3
"""Segment isocortex ROI into 3 parts for CTRL mouse: sub-129-921-R1F"""

import subprocess
import numpy as np
import nibabel as nib
import os

# Subject directory
subj_dir = "/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-129-921-R1F/ses-260115121231/dwi"
os.chdir(subj_dir)

# Create output directory
os.makedirs("isocortex_segments", exist_ok=True)

# Load isocortex ROI - using RIGHT isocortex (appears on right in viewer)
isocortex_file = "atlas/rois/roi_isocortex_right.nii.gz"

# Load the nifti
img = nib.load(isocortex_file)
data = img.get_fdata()
affine = img.affine

print(f"Image shape: {data.shape}")
print(f"Voxel size: {np.abs(np.diag(affine)[:3])}")

# Find non-zero voxels
nz = np.where(data > 0.5)
print(f"Non-zero voxels: {len(nz[0])}")

if len(nz[0]) > 0:
    # Get bounding box
    x_min, x_max = nz[0].min(), nz[0].max()
    y_min, y_max = nz[1].min(), nz[1].max()
    z_min, z_max = nz[2].min(), nz[2].max()
    
    print(f"\nBounding box:")
    print(f"X range (index 0): {x_min} - {x_max}")
    print(f"Y range (index 1): {y_min} - {y_max}")
    print(f"Z range (index 2): {z_min} - {z_max}")
    
    # Calculate thirds for anterior/central/posterior (along X dimension)
    x_range = x_max - x_min
    x_third = x_range / 3
    x_ant_end = x_min + x_third
    x_cent_end = x_min + 2 * x_third
    
    print(f"\nSegmentation boundaries:")
    print(f"Anterior: X <= {x_ant_end:.1f}")
    print(f"Central: {x_ant_end:.1f} < X <= {x_cent_end:.1f}")
    print(f"Posterior: X > {x_cent_end:.1f}")
    
    # Create 3 segments (anterior/central/posterior) - no dorsal/ventral for isocortex
    # Using posterior since that matches the "dorsal anterior" convention visually
    segments = {
        'anterior': (lambda x: x <= x_ant_end),
        'central': (lambda x: (x > x_ant_end) & (x <= x_cent_end)),
        'posterior': (lambda x: x > x_cent_end),
    }
    
    print(f"\nCreating segments:")
    for name, mask_func in segments.items():
        segment = np.zeros_like(data)
        mask = mask_func(nz[0])
        
        for i in range(len(nz[0])):
            if mask[i]:
                segment[nz[0][i], nz[1][i], nz[2][i]] = 1
        
        count = int(np.sum(segment > 0))
        print(f"  {name}: {count} voxels")
        
        # Save segment
        out_img = nib.Nifti1Image(segment.astype(np.float32), affine)
        nib.save(out_img, f"isocortex_segments/isocortex_right_{name}.nii.gz")
    
    print("\nAll segments saved to isocortex_segments/")
    print("Main file for viewing: isocortex_segments/isocortex_right_posterior.nii.gz")
else:
    print("ERROR: No voxels found in isocortex ROI")
