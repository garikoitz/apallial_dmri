from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np


def build_mrconvert_command(reference_mif_path: Path, tmp_nii_path: Path, modules_cmd: str = "") -> str:
    """Build a command that converts reference .mif image to temporary .nii.gz."""
    convert_cmd = f"mrconvert {reference_mif_path} {tmp_nii_path} -force -quiet"
    if modules_cmd:
        return f"{modules_cmd} && {convert_cmd}"
    return convert_cmd


def convert_streamlines_world_to_voxel(streamlines: list[np.ndarray], affine: np.ndarray) -> list[np.ndarray]:
    """Convert streamline points from world (RAS mm) to voxel coordinates."""
    aff_inv = np.linalg.inv(affine)
    vox_streamlines: list[np.ndarray] = []
    for sl in streamlines:
        ones = np.ones((sl.shape[0], 1))
        vox = (aff_inv @ np.hstack([sl, ones]).T).T[:, :3]
        vox_streamlines.append(vox.astype(np.float32))
    return vox_streamlines


def convert_tck_to_trk(
    *,
    tck_path: Path,
    trk_path: Path,
    reference_mif_path: Path,
    modules_cmd: str = "",
    mrconvert_timeout: int = 30,
) -> int:
    """Convert one .tck file to .trk using a .mif reference image for geometry/header."""
    import nibabel as nib
    from nibabel.streamlines import load as tck_load
    from nibabel.streamlines.trk import TrkFile

    with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
        tmp_nii = Path(tmp.name)

    try:
        subprocess.run(
            build_mrconvert_command(reference_mif_path, tmp_nii, modules_cmd=modules_cmd),
            shell=True,
            check=True,
            capture_output=True,
            timeout=mrconvert_timeout,
        )

        ref_img = nib.load(str(tmp_nii))
        affine = ref_img.affine
        zooms = np.array(ref_img.header.get_zooms()[:3], dtype=np.float32)
        dims = np.array(ref_img.shape[:3], dtype=np.int16)

        tck = tck_load(str(tck_path))
        streamlines = list(tck.streamlines)
        vox_streamlines = convert_streamlines_world_to_voxel(streamlines, affine)

        hdr = TrkFile.create_empty_header()
        hdr["dimensions"] = dims
        hdr["voxel_sizes"] = zooms
        hdr["voxel_order"] = b"RAS"
        hdr["vox_to_ras"] = affine.astype(np.float32)

        tractogram = nib.streamlines.Tractogram(vox_streamlines, affine_to_rasmm=affine)
        nib.streamlines.save(TrkFile(tractogram, header=hdr), str(trk_path))
        return len(streamlines)
    finally:
        if os.path.exists(tmp_nii):
            os.remove(tmp_nii)