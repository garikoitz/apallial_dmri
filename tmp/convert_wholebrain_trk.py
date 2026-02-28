#!/usr/bin/env python3
"""Convert wholebrain_100k.tck to .trk for Freeview."""
import subprocess
from pathlib import Path
import numpy as np
import nibabel as nib
from nibabel.streamlines import load as tck_load
from nibabel.streamlines.trk import TrkFile

BASE = Path("data/derivatives/mrtrix_preproc_flip1/sub-Apallial20R1F/ses-240909094232")
TCK = BASE / "dwi/wholebrain/wholebrain_100k.tck"
T2W_NII = Path("tmp/t2w_highres_flip_ap20.nii.gz")
OUT = BASE / "dwi/wholebrain/passthrough_tracts/trk/wholebrain_100k.trk"
FREEVIEW = Path("tmp/freeview_ap20/wholebrain_100k.trk")

ref_img = nib.load(str(T2W_NII))
aff = ref_img.affine
aff_inv = np.linalg.inv(aff)
zooms = np.array(ref_img.header.get_zooms()[:3], dtype=np.float32)
dims  = np.array(ref_img.shape[:3], dtype=np.int16)

print(f"Loading {TCK} ...")
tck = tck_load(str(TCK))
streamlines = list(tck.streamlines)
print(f"  {len(streamlines)} streamlines loaded")

print("Converting world -> voxel coords ...")
vox_streamlines = []
for sl in streamlines:
    ones = np.ones((sl.shape[0], 1))
    vox = (aff_inv @ np.hstack([sl, ones]).T).T[:, :3]
    vox_streamlines.append(vox.astype(np.float32))

hdr = TrkFile.create_empty_header()
hdr['dimensions']  = dims
hdr['voxel_sizes'] = zooms
hdr['voxel_order'] = b'RAS'
hdr['vox_to_ras']  = aff.astype(np.float32)

tractogram = nib.streamlines.Tractogram(vox_streamlines, affine_to_rasmm=aff)
nib.streamlines.save(TrkFile(tractogram, header=hdr), str(OUT))
print(f"  Saved: {OUT}")

# Symlink to freeview folder
import os
if FREEVIEW.exists() or FREEVIEW.is_symlink():
    FREEVIEW.unlink()
FREEVIEW.symlink_to(OUT.resolve())
print(f"  Symlinked: {FREEVIEW}")
