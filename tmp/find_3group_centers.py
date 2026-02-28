#!/usr/bin/env python3
"""Find central AP slice for 3 representative subjects."""
import subprocess
import numpy as np
import nibabel as nib
import tempfile, os

DERIV = "data/derivatives/mrtrix_preproc_flip1"

subjects = {
    "CTRL": ("sub-138-1010-L1F", "ses-260115101441",
             "atlas/rois/fod_space/roi_isocortex_left.nii.gz"),
    "AP":   ("sub-769-C-O5-R1M", "ses-240911115733",
             "manual_roi/fod_space/roi_Graft_left_flip_fod.mif"),
    "APTR": ("sub-Apallial17NmF", "ses-240816091633",
             "manual_roi/fod_space/roi_Graft_left_flip_fod.mif"),
}

def load_any(path):
    if path.endswith('.nii.gz') or path.endswith('.nii'):
        return nib.load(path).get_fdata()
    tmp = tempfile.mktemp(suffix=".nii.gz")
    subprocess.run(["mrconvert", path, tmp, "-quiet", "-force"], check=True)
    data = nib.load(tmp).get_fdata()
    os.unlink(tmp)
    return data

for group, (subj, ses, roi_rel) in subjects.items():
    roi_path = f"{DERIV}/{subj}/{ses}/{roi_rel}"
    fod_path = f"{DERIV}/{subj}/{ses}/dwi/wm_fod_norm_group.mif"
    
    fod_size = subprocess.run(
        ["mrinfo", fod_path, "-size"], capture_output=True, text=True
    ).stdout.strip().split()
    n_ap = int(fod_size[0])
    
    data = load_any(roi_path)
    ap_profile = np.sum(data, axis=(1, 2))
    nonzero = np.where(ap_profile > 0)[0]
    
    if len(nonzero) > 0:
        ap_min, ap_max = nonzero[0], nonzero[-1]
        ap_center = (ap_min + ap_max) // 2
        roi_name = "isocortex" if group == "CTRL" else "graft"
        print(f"{group} ({subj}): FOD AP={n_ap}, {roi_name} AP=[{ap_min}-{ap_max}], "
              f"center={ap_center}, extent={ap_max-ap_min+1}")
    else:
        print(f"{group} ({subj}): No ROI voxels found!")
