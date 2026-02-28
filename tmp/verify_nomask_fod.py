#!/usr/bin/env python3
"""Verify nomask FOD has signal at CST ROI voxels."""
import nibabel as nib
import numpy as np
import subprocess

DERIV = "data/derivatives/mrtrix_preproc_flip1"
subjects = {
    "sub-769-C-O5-R1M": "ses-240911115733",
    "sub-Apallial20R1F": "ses-240909094232",
}

for subj, ses in subjects.items():
    base = f"{DERIV}/{subj}/{ses}/dwi"
    out_dir = f"{base}/wholebrain/passthrough_cst"
    
    # Convert nomask FOD l=0 to nifti
    fod_nii = f"/tmp/fod_nomask_{subj}.nii.gz"
    subprocess.run(["mrconvert", f"{base}/wm_fod_group_nomask.mif", fod_nii,
                    "-coord", "3", "0", "-force", "-quiet"], capture_output=True)
    fod = nib.load(fod_nii).get_fdata().squeeze()
    
    print(f"\n{subj}:")
    for hemi in ["left", "right"]:
        roi_nii = f"/tmp/cst_{hemi}_{subj}.nii.gz"
        subprocess.run(["mrconvert", f"{out_dir}/roi_CST_{hemi}_fod.mif", roi_nii,
                        "-force", "-quiet"], capture_output=True)
        roi = nib.load(roi_nii).get_fdata()
        coords = np.argwhere(roi > 0)
        
        for c in coords:
            i, j, k = c
            print(f"  {hemi} ({i},{j},{k}): FOD_nomask_l0={fod[i,j,k]:.4f}")
