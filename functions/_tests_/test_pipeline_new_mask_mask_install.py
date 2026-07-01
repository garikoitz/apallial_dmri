from pathlib import Path

from functions.pipeline_new_mask.mask_install import (
    build_mask_install_plan,
    resolve_dwi_template,
    select_new_mask,
)


def _subject(base: Path) -> dict:
    dwi = base / "sub-01" / "ses-1" / "dwi"
    dwi.mkdir(parents=True, exist_ok=True)
    return {
        "id": "sub-01",
        "group": "CTRL",
        "dwi_dir": dwi,
        "mask_new_edited": dwi / "mask_new_edited.mif",
        "mask_new_mif": dwi / "mask_new.mif",
    }


def test_select_new_mask_prefers_edited(tmp_path: Path):
    s = _subject(tmp_path)
    s["mask_new_edited"].touch()
    s["mask_new_mif"].touch()

    path, source = select_new_mask(s)
    assert path == s["mask_new_edited"]
    assert source == "edited"


def test_select_new_mask_falls_back_to_atlas(tmp_path: Path):
    s = _subject(tmp_path)
    s["mask_new_mif"].touch()

    path, source = select_new_mask(s)
    assert path == s["mask_new_mif"]
    assert source == "atlas-based"


def test_resolve_dwi_template_priority(tmp_path: Path):
    dwi = tmp_path / "dwi"
    dwi.mkdir()
    p = dwi / "dwi_flip_gradcheck_denoise_bias_anat.mif"
    p.touch()
    assert resolve_dwi_template(dwi) == p


def test_resolve_dwi_template_fallback(tmp_path: Path):
    dwi = tmp_path / "dwi"
    dwi.mkdir()
    f = dwi / "wm_fod_norm.mif"
    f.touch()
    assert resolve_dwi_template(dwi) == f


def test_build_mask_install_plan_fields(tmp_path: Path):
    s = _subject(tmp_path)
    s["mask_new_mif"].touch()
    (s["dwi_dir"] / "dwi_preproc.mif").touch()

    plan = build_mask_install_plan([s], backup_suffix="20260701")[0]

    assert plan["subject_id"] == "sub-01"
    assert plan["source"] == "atlas-based"
    assert plan["mask_mif"].name == "mask.mif"
    assert plan["backup_mask"].name == "mask_old_20260701.mif"
    assert plan["dwi_template"].name == "dwi_preproc.mif"
