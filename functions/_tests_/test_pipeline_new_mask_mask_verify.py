from pathlib import Path

from functions.pipeline_new_mask.mask_install import build_mask_grid_report, compute_pipeline_checkpoint


def _subject(base: Path, sid: str = "sub-01") -> dict:
    dwi = base / sid / "ses-1" / "dwi"
    dwi.mkdir(parents=True, exist_ok=True)
    wb = dwi / "wholebrain"
    wb.mkdir(parents=True, exist_ok=True)
    return {
        "id": sid,
        "group": "CTRL",
        "dwi_dir": dwi,
        "wb_tck": wb / "wholebrain_1000k.tck",
    }


def test_build_mask_grid_report_missing_mask(tmp_path: Path):
    s = _subject(tmp_path)
    rows = build_mask_grid_report([s], size_reader=lambda _p: ["28", "100", "107"])
    assert rows[0]["subject_id"] == "sub-01"
    assert rows[0]["missing_mask"] is True


def test_build_mask_grid_report_compare_with_fod(tmp_path: Path):
    s = _subject(tmp_path)
    (s["dwi_dir"] / "mask.mif").touch()
    (s["dwi_dir"] / "wm_fod_norm.mif").touch()

    def size_reader(path: Path):
        if path.name == "mask.mif":
            return ["28", "100", "107"]
        return ["28", "100", "107", "45"]

    row = build_mask_grid_report([s], size_reader=size_reader)[0]
    assert row["missing_mask"] is False
    assert row["has_fod"] is True
    assert row["match"] is True


def test_compute_pipeline_checkpoint_classification(tmp_path: Path):
    s = _subject(tmp_path)
    mask = s["dwi_dir"] / "mask.mif"
    mask_dws = s["dwi_dir"] / "mask_new_dwi_space.mif"
    fod = s["dwi_dir"] / "wm_fod_norm.mif"
    grp = s["dwi_dir"] / "wm_fod_norm_group.mif"
    tck = s["wb_tck"]
    for f in [mask, mask_dws, fod, grp, tck]:
        f.touch()

    # Ensure mtimes satisfy "ok" conditions.
    mask_dws.touch()
    fod.touch()
    grp.touch()
    tck.touch()

    chk = compute_pipeline_checkpoint([s], size_reader=lambda _p: ["28", "100", "107"])
    assert chk["mask_ok"] == ["sub-01"]
    assert chk["fod_ok"] == ["sub-01"]
    assert chk["grp_ok"] == ["sub-01"]
    assert chk["tck_ok"] == ["sub-01"]