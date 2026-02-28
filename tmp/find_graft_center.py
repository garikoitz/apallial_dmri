#!/usr/bin/env python3
"""Find central AP slice of the graft for each APTR mouse."""
import subprocess
import numpy as np

DERIV = "data/derivatives/mrtrix_preproc_flip1"
subjects = {
    "sub-AP72-862-901-NmF": "ses-251210093342",
    "sub-AP72-862-901-R1M": "ses-251210100043",
    "sub-Apallial17NmF": "ses-240816091633",
    "sub-Apallial20R1F": "ses-240909094232",
}

import nibabel as nib
import tempfile, os

for subj, ses in subjects.items():
    graft = f"{DERIV}/{subj}/{ses}/manual_roi/fod_space/roi_Graft_left_flip_fod.mif"
    fod = f"{DERIV}/{subj}/{ses}/dwi/wm_fod_norm_group.mif"
    
    # Get FOD dimensions
    fod_size = subprocess.run(
        ["mrinfo", fod, "-size"], capture_output=True, text=True
    ).stdout.strip().split()
    n_ap = int(fod_size[0])
    
    # Convert graft to nifti, load with nibabel
    tmp_nii = tempfile.mktemp(suffix=".nii.gz")
    subprocess.run(["mrconvert", graft, tmp_nii, "-quiet", "-force"], check=True)
    img = nib.load(tmp_nii)
    data = img.get_fdata()
    os.unlink(tmp_nii)
    
    # MRtrix axis 0 = AP corresponds to nifti axis 0 (after mrconvert)
    # Sum along DV (axis 1) and ML (axis 2) to get AP profile
    ap_profile = np.sum(data, axis=(1, 2))
    
    nonzero = np.where(ap_profile > 0)[0]
    if len(nonzero) > 0:
        ap_min, ap_max = nonzero[0], nonzero[-1]
        ap_center = (ap_min + ap_max) // 2
        ap_com = np.average(nonzero, weights=ap_profile[nonzero])
        print(f"{subj}: FOD AP={n_ap}, Graft AP range=[{ap_min}-{ap_max}], "
              f"center={ap_center}, CoM={ap_com:.1f}, extent={ap_max-ap_min+1} slices")
    else:
        print(f"{subj}: No graft voxels found!")
