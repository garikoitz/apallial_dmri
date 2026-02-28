#!/usr/bin/env python3
"""Convert FOD central voxel coords to T2w voxel coords for mrview."""
import subprocess, numpy as np

DERIV = "data/derivatives/mrtrix_preproc_flip1"
subjects = [
    ("CTRL", "sub-138-1010-L1F", "ses-260115101441", 13),
    ("AP",   "sub-769-C-O5-R1M", "ses-240911115733", 14),
    ("APTR", "sub-Apallial17NmF", "ses-240816091633", 15),
]

def get_transform(path):
    lines = subprocess.run(
        ["mrinfo", path, "-transform"],
        capture_output=True, text=True
    ).stdout.strip().split('\n')
    aff = np.zeros((4, 4))
    for i, line in enumerate(lines):
        aff[i] = [float(x) for x in line.split()]
    aff[3] = [0, 0, 0, 1]
    return aff

def get_size(path):
    return [int(x) for x in subprocess.run(
        ["mrinfo", path, "-size"],
        capture_output=True, text=True
    ).stdout.strip().split()]

for group, subj, ses, fod_slice in subjects:
    base = f"{DERIV}/{subj}/{ses}/dwi"
    fod_path = f"{base}/wm_fod_norm_group.mif"
    t2_path = f"{base}/t2w_highres_flip.mif"
    
    fod_aff = get_transform(fod_path)
    t2_aff = get_transform(t2_path)
    fod_size = get_size(fod_path)
    t2_size = get_size(t2_path)
    
    # FOD center voxel
    fod_vox = np.array([fod_slice, fod_size[1] // 2, fod_size[2] // 2, 1.0])
    # Convert to world coords, then to T2w voxels
    world = fod_aff @ fod_vox
    t2_vox = np.linalg.inv(t2_aff) @ world
    
    print(f"{group}: FOD voxel=({fod_slice},{fod_size[1]//2},{fod_size[2]//2}) "
          f"-> T2w voxel=({int(round(t2_vox[0]))},{int(round(t2_vox[1]))},{int(round(t2_vox[2]))})")
    print(f"  T2w size={t2_size[:3]}, FOD size={fod_size[:3]}")
