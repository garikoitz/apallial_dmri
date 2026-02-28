#!/usr/bin/env python3
"""Convert wholebrain_5k.tck to .trk for Freeview."""
import sys, os
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from pathlib import Path
import numpy as np
import nibabel as nib
from nibabel.streamlines import load as tck_load
from nibabel.streamlines.trk import TrkFile

BASE = Path("data/derivatives/mrtrix_preproc_flip1/sub-Apallial20R1F/ses-240909094232")
TCK  = BASE / "dwi/wholebrain/wholebrain_5k.tck"
T2W  = Path("tmp/t2w_highres_flip_ap20.nii.gz")
OUT  = BASE / "dwi/wholebrain/passthrough_tracts/trk/wholebrain_5k.trk"
LINK = Path("tmp/freeview_ap20/wholebrain_5k.trk")

ref = nib.load(str(T2W))
aff = ref.affine
aff_inv = np.linalg.inv(aff)

tck = tck_load(str(TCK))
streamlines = list(tck.streamlines)
print(f"{len(streamlines)} streamlines")

vox_sl = []
for sl in streamlines:
    ones = np.ones((sl.shape[0], 1))
    vox = (aff_inv @ np.hstack([sl, ones]).T).T[:, :3]
    vox_sl.append(vox.astype(np.float32))

hdr = TrkFile.create_empty_header()
hdr['dimensions']  = np.array(ref.shape[:3], dtype=np.int16)
hdr['voxel_sizes'] = np.array(ref.header.get_zooms()[:3], dtype=np.float32)
hdr['voxel_order'] = b'RAS'
hdr['vox_to_ras']  = aff.astype(np.float32)

tractogram = nib.streamlines.Tractogram(vox_sl, affine_to_rasmm=aff)
nib.streamlines.save(TrkFile(tractogram, header=hdr), str(OUT))
print(f"Saved: {OUT}")

if LINK.exists() or LINK.is_symlink():
    LINK.unlink()
LINK.symlink_to(OUT.resolve())
print(f"Symlinked: {LINK}")
