#!/usr/bin/env python3
"""Convert all _cleaned.tck tracts for sub-Apallial20R1F to .trk (Freeview-compatible).

Uses T2w_highres_flip as the reference image for the .trk header.
Output goes in the same folder as the .tck files, with .trk extension.
"""
import subprocess, tempfile, os, glob
from pathlib import Path
import numpy as np
import nibabel as nib
from nibabel.streamlines import load as tck_load, save as trk_save
from nibabel.streamlines.trk import TrkFile

BASE = Path("data/derivatives/mrtrix_preproc_flip1/sub-Apallial20R1F/ses-240909094232")
TCK_DIR = BASE / "dwi/wholebrain/passthrough_tracts"
T2W_MIF = BASE / "dwi/t2w_highres_flip.mif"
OUT_DIR = TCK_DIR / "trk"
OUT_DIR.mkdir(exist_ok=True)

# Step 1: Convert T2w to NIfTI (needed for trk header)
print("Converting T2w to NIfTI...")
t2w_nii = Path("tmp/t2w_highres_flip_ap20.nii.gz")
subprocess.run(
    ["micromamba", "run", "-n", "mbrain",
     "mrconvert", str(T2W_MIF), str(t2w_nii), "-force", "-quiet"],
    check=True
)
ref_img = nib.load(str(t2w_nii))
print(f"  Reference shape: {ref_img.shape}, voxel size: {ref_img.header.get_zooms()}")

# Step 2: Find all _cleaned.tck files
tck_files = sorted(TCK_DIR.glob("*_cleaned.tck"))
print(f"\nFound {len(tck_files)} _cleaned.tck files")

# Step 3: Convert each
success, failed = 0, 0
aff = ref_img.affine
aff_inv = np.linalg.inv(aff)
zooms = np.array(ref_img.header.get_zooms()[:3], dtype=np.float32)
dims = np.array(ref_img.shape[:3], dtype=np.int16)

for tck_path in tck_files:
    out_path = OUT_DIR / tck_path.with_suffix(".trk").name
    try:
        # Load .tck — streamlines are in world (RAS mm) space
        tck = tck_load(str(tck_path))
        streamlines = list(tck.streamlines)

        if len(streamlines) == 0:
            print(f"  ⚠ EMPTY: {tck_path.name}")
            continue

        # Convert world (mm) → voxel coordinates  (required by TrkFile)
        vox_streamlines = []
        for sl in streamlines:
            ones = np.ones((sl.shape[0], 1))
            vox = (aff_inv @ np.hstack([sl, ones]).T).T[:, :3]
            vox_streamlines.append(vox.astype(np.float32))

        # Build TrkFile header
        hdr = TrkFile.create_empty_header()
        hdr['dimensions']  = dims
        hdr['voxel_sizes'] = zooms
        hdr['voxel_order'] = b'RAS'
        hdr['vox_to_ras']  = aff.astype(np.float32)

        tractogram = nib.streamlines.Tractogram(vox_streamlines,
                                                 affine_to_rasmm=aff)
        nib.streamlines.save(TrkFile(tractogram, header=hdr), str(out_path))

        print(f"  ✓ {tck_path.name} → {out_path.name}  ({len(streamlines)} streamlines)")
        success += 1

    except Exception as e:
        print(f"  ✗ {tck_path.name}: {e}")
        failed += 1

print(f"\nDone: {success} converted, {failed} failed")
print(f"Output: {OUT_DIR.resolve()}")
