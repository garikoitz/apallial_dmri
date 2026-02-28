#!/usr/bin/env python3
"""Check if CST ROI voxels have DWI/FOD signal (even if outside mask)."""
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
    fod_path = f"{base}/dwi/wm_fod_norm_group.mif"
    dwi_path = f"{base}/dwi/dwi_flip_gradcheck_denoise_bias_anat.mif"
    
    print(f"\n{'='*60}")
    print(f"Subject: {subj}")
    print(f"{'='*60}")
    
    # Convert mask, FOD (first volume), DWI (first volume) to nifti
    mask_nii = f"/tmp/mask_{subj}.nii.gz"
    fod_nii = f"/tmp/fod_{subj}.nii.gz"
    dwi_nii = f"/tmp/dwi_{subj}.nii.gz"
    
    subprocess.run(["mrconvert", mask_path, mask_nii, "-force", "-quiet"], capture_output=True)
    # Extract first FOD volume (l=0 component)
    subprocess.run(["mrconvert", fod_path, fod_nii, "-coord", "3", "0", "-force", "-quiet"], capture_output=True)
    # Extract first DWI volume (b0)
    subprocess.run(["mrconvert", dwi_path, dwi_nii, "-coord", "3", "0", "-force", "-quiet"], capture_output=True)
    
    mask = nib.load(mask_nii).get_fdata()
    fod = nib.load(fod_nii).get_fdata().squeeze()
    dwi = nib.load(dwi_nii).get_fdata().squeeze()
    
    print(f"Mask shape: {mask.shape}, non-zero: {np.sum(mask > 0)}")
    print(f"FOD shape: {fod.shape}")
    print(f"DWI (b0) shape: {dwi.shape}")
    
    for hemi in ["left", "right"]:
        roi_fod = f"{out_dir}/roi_CST_{hemi}_fod.mif"
        if not os.path.exists(roi_fod):
            print(f"  {hemi}: ROI not found")
            continue
            
        roi_nii = f"/tmp/cst_{hemi}_{subj}.nii.gz"
        subprocess.run(["mrconvert", roi_fod, roi_nii, "-force", "-quiet"], capture_output=True)
        roi = nib.load(roi_nii).get_fdata()
        coords = np.argwhere(roi > 0)
        
        print(f"\n  {hemi} CST ROI ({len(coords)} voxels):")
        for c in coords:
            i, j, k = c
            m = mask[i, j, k] > 0
            f = fod[i, j, k] if i < fod.shape[0] and j < fod.shape[1] and k < fod.shape[2] else 0
            d = dwi[i, j, k] if i < dwi.shape[0] and j < dwi.shape[1] and k < dwi.shape[2] else 0
            print(f"    ({i:3d},{j:3d},{k:3d})  mask={m}  FOD_l0={f:.4f}  DWI_b0={d:.1f}")
        
        # Also check neighborhood: expand 2 voxels around ROI
        if len(coords) > 0:
            center = coords.mean(axis=0).astype(int)
            i, j, k = center
            # Check mask coverage in neighborhood
            r = 3
            neighborhood = mask[max(0,i-r):i+r+1, max(0,j-r):j+r+1, max(0,k-r):k+r+1]
            fod_neighborhood = fod[max(0,i-r):i+r+1, max(0,j-r):j+r+1, max(0,k-r):k+r+1]
            dwi_neighborhood = dwi[max(0,i-r):i+r+1, max(0,j-r):j+r+1, max(0,k-r):k+r+1]
            print(f"    Center: ({i},{j},{k}), ±{r} neighborhood:")
            print(f"      Mask voxels: {np.sum(neighborhood > 0)} / {neighborhood.size}")
            print(f"      FOD range: [{fod_neighborhood.min():.4f}, {fod_neighborhood.max():.4f}]")
            print(f"      DWI range: [{dwi_neighborhood.min():.1f}, {dwi_neighborhood.max():.1f}]")
