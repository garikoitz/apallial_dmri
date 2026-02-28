#!/usr/bin/env python3
"""Check where CST ROI voxels are located in FOD space and if they overlap with brain mask."""
import nibabel as nib
import numpy as np
import subprocess
import os

DERIV = "data/derivatives/mrtrix_preproc_flip1"

subjects = {
    "sub-769-C-O5-R1M": "ses-240911115733",
    "sub-Apallial20R1F": "ses-240909094232",
}

for subj, ses in subjects.items():
    base = f"{DERIV}/{subj}/{ses}"
    out_dir = f"{base}/dwi/wholebrain/passthrough_cst"
    mask_path = f"{base}/dwi/mask_anat.mif"
    
    print(f"\n{'='*60}")
    print(f"Subject: {subj}")
    print(f"{'='*60}")
    
    # Convert mask to nifti
    mask_nii = f"/tmp/mask_{subj}.nii.gz"
    subprocess.run(["mrconvert", mask_path, mask_nii, "-force", "-quiet"],
                   capture_output=True)
    mask = nib.load(mask_nii).get_fdata()
    print(f"Mask shape: {mask.shape}, non-zero: {np.sum(mask > 0)}")
    
    for hemi in ["left", "right"]:
        roi_fod = f"{out_dir}/roi_CST_{hemi}_fod.mif"
        if not os.path.exists(roi_fod):
            print(f"  {hemi}: ROI not found")
            continue
            
        roi_nii = f"/tmp/cst_{hemi}_{subj}.nii.gz"
        subprocess.run(["mrconvert", roi_fod, roi_nii, "-force", "-quiet"],
                       capture_output=True)
        roi = nib.load(roi_nii).get_fdata()
        coords = np.argwhere(roi > 0)
        
        print(f"\n  {hemi} CST ROI:")
        print(f"    Shape: {roi.shape}, Non-zero voxels: {len(coords)}")
        for c in coords:
            i, j, k = c
            in_mask = mask[i, j, k] > 0 if i < mask.shape[0] and j < mask.shape[1] and k < mask.shape[2] else False
            print(f"    Voxel ({i:3d}, {j:3d}, {k:3d}) → in mask: {in_mask}")
        
        # Also check graft ROI
        graft_path = f"{base}/manual_roi/fod_space/roi_Graft_{hemi}_flip_fod.mif"
        if os.path.exists(graft_path):
            graft_nii = f"/tmp/graft_{hemi}_{subj}.nii.gz"
            subprocess.run(["mrconvert", graft_path, graft_nii, "-force", "-quiet"],
                           capture_output=True)
            graft = nib.load(graft_nii).get_fdata()
            graft_coords = np.argwhere(graft > 0)
            print(f"    Graft ROI: {len(graft_coords)} voxels")
            # Check overlap
            overlap = 0
            for c in coords:
                i, j, k = c
                if i < graft.shape[0] and j < graft.shape[1] and k < graft.shape[2]:
                    if graft[i, j, k] > 0:
                        overlap += 1
            print(f"    CST-Graft overlap: {overlap} voxels")
            # Find nearest graft voxel
            if len(coords) > 0 and len(graft_coords) > 0:
                from scipy.spatial.distance import cdist
                dists = cdist(coords, graft_coords)
                min_dist = dists.min()
                print(f"    Min distance CST→Graft: {min_dist:.1f} voxels")
