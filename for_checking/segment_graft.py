#!/usr/bin/env python3
"""Segment graft ROI into 6 parts: dorsal/ventral × anterior/central/posterior"""

import subprocess
import numpy as np
import nibabel as nib
import os
import sys

# Subject directory
subj_dir = "/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-AP72-862-901-NmF/ses-251210093342/dwi"
os.chdir(subj_dir)

# Create output directory
os.makedirs("graft_segments", exist_ok=True)

# Load graft ROI - using LEFT graft (which appears on RIGHT hemisphere visually)
graft_file = "../manual_roi/flipped/roi_Graft_left_flip.mif"

# Convert to nifti temporarily
subprocess.run(["mrconvert", graft_file, "graft_segments/graft_left_temp.nii.gz", "-force"], 
               capture_output=True)

# Load the nifti
img = nib.load("graft_segments/graft_left_temp.nii.gz")
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
    
    # Calculate half for dorsal/ventral (along Y dimension)
    y_range = y_max - y_min
    y_mid = y_min + y_range / 2
    
    print(f"\nSegmentation boundaries:")
    print(f"Anterior: X <= {x_ant_end:.1f}")
    print(f"Central: {x_ant_end:.1f} < X <= {x_cent_end:.1f}")
    print(f"Posterior: X > {x_cent_end:.1f}")
    print(f"Dorsal: Y >= {y_mid:.1f}")
    print(f"Ventral: Y < {y_mid:.1f}")
    
    # Create all 6 segments
    segments = {
        'anterior_dorsal': (lambda x, y: (x <= x_ant_end) & (y >= y_mid)),
        'anterior_ventral': (lambda x, y: (x <= x_ant_end) & (y < y_mid)),
        'central_dorsal': (lambda x, y: (x > x_ant_end) & (x <= x_cent_end) & (y >= y_mid)),
        'central_ventral': (lambda x, y: (x > x_ant_end) & (x <= x_cent_end) & (y < y_mid)),
        'posterior_dorsal': (lambda x, y: (x > x_cent_end) & (y >= y_mid)),
        'posterior_ventral': (lambda x, y: (x > x_cent_end) & (y < y_mid)),
    }
    
    print(f"\nCreating segments:")
    for name, mask_func in segments.items():
        segment = np.zeros_like(data)
        mask = mask_func(nz[0], nz[1])
        
        for i in range(len(nz[0])):
            if mask[i]:
                segment[nz[0][i], nz[1][i], nz[2][i]] = 1
        
        count = int(np.sum(segment > 0))
        print(f"  {name}: {count} voxels")
        
        # Save segment - using "left" in filename (appears on right hemisphere)
        out_img = nib.Nifti1Image(segment.astype(np.float32), affine)
        nib.save(out_img, f"graft_segments/graft_left_{name}.nii.gz")
    
    print("\nAll segments saved to graft_segments/")
    print("Main file for viewing: graft_segments/graft_left_posterior_ventral.nii.gz")
else:
    print("ERROR: No voxels found in graft ROI")
    sys.exit(1)
