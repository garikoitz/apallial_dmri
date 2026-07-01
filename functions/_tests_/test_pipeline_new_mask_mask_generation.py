from pathlib import Path

from functions.pipeline_new_mask.mask_generation import (
    DEFAULT_GROUP_ATLAS_IDS,
    build_group_atlas_info,
    build_mask_generation_args,
)


def _subject(sid: str, group: str, base: Path) -> dict:
    dwi = base / sid / "ses-a" / "dwi"
    wb = dwi / "wholebrain"
    return {
        "id": sid,
        "group": group,
        "session": "ses-a",
        "base": base / sid / "ses-a",
        "dwi_dir": dwi,
        "wb_dir": wb,
        "t2w": dwi / "t2w_highres_flip.mif",
        "mask_new_edited": dwi / "mask_new_edited.mif",
        "mask_new": dwi / "mask_new.nii.gz",
        "mask_new_mif": dwi / "mask_new.mif",
    }


def test_build_group_atlas_info_matches_expected_paths(tmp_path: Path):
    subjects = [
        _subject(DEFAULT_GROUP_ATLAS_IDS["CTRL"], "CTRL", tmp_path),
        _subject(DEFAULT_GROUP_ATLAS_IDS["AP"], "AP", tmp_path),
        _subject(DEFAULT_GROUP_ATLAS_IDS["APTR"], "APTR", tmp_path),
    ]

    atlas_info = build_group_atlas_info(subjects, DEFAULT_GROUP_ATLAS_IDS)

    assert set(atlas_info.keys()) == {"CTRL", "AP", "APTR"}
    assert atlas_info["CTRL"]["atlas_id"] == DEFAULT_GROUP_ATLAS_IDS["CTRL"]
    assert atlas_info["AP"]["mask_nii"].name == "mask_new_edited.nii.gz"
    assert atlas_info["APTR"]["t2w_n4"].name == "_atlas_t2w_n4.nii.gz"


def test_build_mask_generation_args_uses_group_atlas_resources(tmp_path: Path):
    ctrl_atlas = _subject(DEFAULT_GROUP_ATLAS_IDS["CTRL"], "CTRL", tmp_path)
    ctrl_target = _subject("sub-ctrl-target", "CTRL", tmp_path)
    ap_atlas = _subject(DEFAULT_GROUP_ATLAS_IDS["AP"], "AP", tmp_path)
    aptr_atlas = _subject(DEFAULT_GROUP_ATLAS_IDS["APTR"], "APTR", tmp_path)

    subjects = [ctrl_atlas, ctrl_target, ap_atlas, aptr_atlas]
    atlas_info = build_group_atlas_info(subjects, DEFAULT_GROUP_ATLAS_IDS)
    args = build_mask_generation_args(subjects, atlas_info, modules="module load ants")

    assert len(args) == 4
    target_tuple = next(t for t in args if t[0] == "sub-ctrl-target")
    assert target_tuple[1] == "CTRL"
    assert target_tuple[6].endswith("_atlas_t2w_n4.nii.gz")
    assert target_tuple[8] == "module load ants"
