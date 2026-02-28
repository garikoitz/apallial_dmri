#!/usr/bin/env python3
"""Create shrunk Graft ROI files in manual_roi/fod_space/shrink/ for all 4 subjects."""
import os, subprocess, tempfile
import nibabel as nib
import numpy as np
from scipy import ndimage
from pathlib import Path

def load_mif_as_array(path):
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    subprocess.run(['mrconvert', str(path), tmp_path, '-force', '-quiet'], check=True, capture_output=True)
    img = nib.load(tmp_path)
    data = img.get_fdata()
    os.remove(tmp_path)
    return data, img

def save_array_as_mif(data, ref_img, out_path):
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    new_img = nib.Nifti1Image(data.astype(np.float32), ref_img.affine, ref_img.header)
    nib.save(new_img, tmp_path)
    subprocess.run(['mrconvert', tmp_path, str(out_path), '-force', '-quiet'], check=True, capture_output=True)
    os.remove(tmp_path)

struct = ndimage.generate_binary_structure(3, 1)

subjects = [
    ("sub-769-C-O5-R1M", "ses-240911115733"),
    ("sub-APtestmouseNmM", "ses-240911111853"),
    ("sub-Apallial17NmF", "ses-240816091633"),
    ("sub-Apallial20R1F", "ses-240909094232"),
]

BASE = Path("/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1")

for subj, sess in subjects:
    shrink_dir = BASE / subj / sess / "manual_roi/fod_space/shrink"
    shrink_dir.mkdir(parents=True, exist_ok=True)
    fod_dir = BASE / subj / sess / "manual_roi/fod_space"
    
    for hemi in ['left', 'right']:
        roi_path = fod_dir / f"roi_Graft_{hemi}_flip_fod.mif"
        out_path = shrink_dir / f"roi_Graft_{hemi}_flip_fod_shrink.mif"
        
        if not roi_path.exists():
            print(f"  MISSING: {roi_path.name}")
            continue
        
        data, img = load_mif_as_array(roi_path)
        mask = data > 0.5
        eroded = ndimage.binary_erosion(mask, structure=struct)
        save_array_as_mif(eroded, img, out_path)
        print(f"  {subj} {hemi}: {int(mask.sum())} -> {int(eroded.sum())} voxels -> {out_path.name}")

print("Done!")
