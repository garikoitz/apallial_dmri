import numpy as np
from pathlib import Path

from functions.pipeline_new_mask.trk_conversion import build_mrconvert_command, convert_streamlines_world_to_voxel


def test_build_mrconvert_command_without_modules():
    cmd = build_mrconvert_command(Path("/a/ref.mif"), Path("/tmp/ref.nii.gz"))
    assert cmd == "mrconvert /a/ref.mif /tmp/ref.nii.gz -force -quiet"


def test_build_mrconvert_command_with_modules():
    cmd = build_mrconvert_command(
        Path("/a/ref.mif"),
        Path("/tmp/ref.nii.gz"),
        modules_cmd="module load mrtrix3",
    )
    assert cmd.startswith("module load mrtrix3 && mrconvert")


def test_convert_streamlines_world_to_voxel_identity_affine():
    sl = np.array([[1.0, 2.0, 3.0], [4.5, 5.5, 6.5]], dtype=np.float32)
    out = convert_streamlines_world_to_voxel([sl], np.eye(4, dtype=np.float32))
    assert np.allclose(out[0], sl)
    assert out[0].dtype == np.float32


def test_convert_streamlines_world_to_voxel_scaled_affine():
    # World points are in mm; with voxel size 2 mm, voxel coordinates are /2.
    aff = np.array(
        [
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=np.float32,
    )
    sl = np.array([[2.0, 4.0, 6.0], [10.0, 12.0, 14.0]], dtype=np.float32)
    out = convert_streamlines_world_to_voxel([sl], aff)
    expected = np.array([[1.0, 2.0, 3.0], [5.0, 6.0, 7.0]], dtype=np.float32)
    assert np.allclose(out[0], expected)