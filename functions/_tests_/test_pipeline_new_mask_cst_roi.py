from pathlib import Path

from functions.pipeline_new_mask.cst_roi import (
    build_temp_cst_exclusion_tasks,
    build_cst_final_qc_summary,
    build_cst_regrid_plan,
    collect_cst_cleaned_counts,
    collect_cst_raw_tract_counts,
    collect_temp_exclusion_rois,
    collect_cst_tckgen_counts,
    collect_cst_status,
    resolve_cst_source_rois,
)


def _subject(base: Path, sid: str = "sub-01", group: str = "AP") -> dict:
    subject_base = base / sid / "ses-1"
    return {
        "id": sid,
        "group": group,
        "base": subject_base,
        "wb_dir": subject_base / "dwi" / "wholebrain",
    }


def test_resolve_cst_source_rois_prefers_manual(tmp_path: Path):
    s = _subject(tmp_path)
    manual = s["base"] / "manual_roi"
    legacy = tmp_path / "legacy"
    manual.mkdir(parents=True)
    legacy.mkdir(parents=True)

    left = manual / "roi_CST_left_flip.mif"
    right = manual / "roi_CST_right_flip.mif"
    left.touch()
    right.touch()

    out_left, out_right = resolve_cst_source_rois(s["id"], manual, legacy)
    assert out_left == left
    assert out_right == right


def test_resolve_cst_source_rois_fallback_legacy(tmp_path: Path):
    s = _subject(tmp_path, sid="sub-A")
    manual = s["base"] / "manual_roi"
    legacy = tmp_path / "legacy"
    manual.mkdir(parents=True)
    legacy.mkdir(parents=True)

    left = legacy / f"roi_CST_left_{s['id']}.mif"
    right = legacy / f"roi_CST_right_{s['id']}_CST.mif"
    left.touch()
    right.touch()

    out_left, out_right = resolve_cst_source_rois(s["id"], manual, legacy)
    assert out_left == left
    assert out_right == right


def test_collect_cst_status_classifies_subjects(tmp_path: Path):
    legacy = tmp_path / "legacy"
    legacy.mkdir(parents=True)

    s_complete = _subject(tmp_path, sid="sub-complete")
    s_regrid = _subject(tmp_path, sid="sub-regrid")
    s_missing = _subject(tmp_path, sid="sub-missing")

    complete_dir = s_complete["wb_dir"] / "passthrough_cst"
    complete_dir.mkdir(parents=True)
    (complete_dir / "roi_CST_left_fod.mif").touch()
    (complete_dir / "roi_CST_right_fod.mif").touch()

    regrid_manual = s_regrid["base"] / "manual_roi"
    regrid_manual.mkdir(parents=True)
    (regrid_manual / "roi_CST_left_flip.mif").touch()
    (regrid_manual / "roi_CST_right_flip.mif").touch()

    out = collect_cst_status([s_complete, s_regrid, s_missing], legacy_dir=legacy)
    assert len(out["complete"]) == 1
    assert len(out["needs_regrid"]) == 1
    assert len(out["missing"]) == 1


def test_build_cst_regrid_plan_returns_expected_lists(tmp_path: Path):
    legacy = tmp_path / "legacy"
    legacy.mkdir(parents=True)

    s_done = _subject(tmp_path, sid="sub-done")
    s_plan = _subject(tmp_path, sid="sub-plan")
    s_miss = _subject(tmp_path, sid="sub-miss")

    done_dir = s_done["wb_dir"] / "passthrough_cst"
    done_dir.mkdir(parents=True)
    (done_dir / "roi_CST_left_fod.mif").touch()
    (done_dir / "roi_CST_right_fod.mif").touch()

    plan_manual = s_plan["base"] / "manual_roi"
    plan_manual.mkdir(parents=True)
    (plan_manual / "roi_CST_left_flip.mif").touch()
    (plan_manual / "roi_CST_right_flip.mif").touch()

    out = build_cst_regrid_plan([s_done, s_plan, s_miss], legacy_dir=legacy)
    assert len(out["already_complete"]) == 1
    assert len(out["to_regrid"]) == 1
    assert len(out["missing_sources"]) == 1


def test_build_cst_final_qc_summary_uses_manual_or_legacy_sources(tmp_path: Path):
    legacy = tmp_path / "legacy"
    legacy.mkdir(parents=True)

    s_manual = _subject(tmp_path, sid="sub-manual")
    s_legacy = _subject(tmp_path, sid="sub-legacy")

    mdir = s_manual["base"] / "manual_roi"
    mdir.mkdir(parents=True)
    (mdir / "roi_CST_left_flip.mif").touch()
    (mdir / "roi_CST_right_flip.mif").touch()
    (s_manual["wb_dir"] / "passthrough_cst").mkdir(parents=True)
    (s_manual["wb_dir"] / "passthrough_cst" / "roi_CST_left_fod.mif").touch()
    (s_manual["wb_dir"] / "passthrough_cst" / "roi_CST_right_fod.mif").touch()

    (legacy / f"roi_CST_left_{s_legacy['id']}.mif").touch()
    (legacy / f"roi_CST_right_{s_legacy['id']}_CST.mif").touch()
    (s_legacy["wb_dir"] / "passthrough_cst").mkdir(parents=True)

    out = build_cst_final_qc_summary([s_manual, s_legacy], legacy_dir=legacy)
    assert out["all_ok"] is False
    rows = {r["subject"]["id"]: r for r in out["rows"]}
    assert rows["sub-manual"]["has_src"] is True
    assert rows["sub-manual"]["has_fod_left"] is True
    assert rows["sub-legacy"]["has_src"] is True
    assert rows["sub-legacy"]["has_fod_left"] is False


def test_collect_cst_tckgen_counts_reports_found_and_missing(tmp_path: Path):
    s = _subject(tmp_path, sid="sub-ap")
    cst_dir = s["wb_dir"] / "passthrough_cst"
    cst_dir.mkdir(parents=True)
    found = cst_dir / "tckgen_CST_left_to_graft_left.tck"
    found.touch()

    def _fake_runner(path: Path) -> str:
        assert path == found
        return "123"

    rows = collect_cst_tckgen_counts([s], run_tckinfo=_fake_runner)
    left = next(r for r in rows if r["hemi"] == "left")
    right = next(r for r in rows if r["hemi"] == "right")
    assert left["exists"] is True and left["count"] == "123"
    assert right["exists"] is False and right["count"] is None


def test_collect_temp_exclusion_rois_found_and_missing(tmp_path: Path):
    roi_dir = tmp_path / "rois"
    roi_dir.mkdir(parents=True)
    (roi_dir / "roi_Striatum_left_flip_fod.mif").touch()
    (roi_dir / "roi_Thalamus_right_flip_fod.mif").touch()

    out = collect_temp_exclusion_rois(roi_dir)
    found_names = {p.name for p in out["found"]}
    assert "roi_Striatum_left_flip_fod.mif" in found_names
    assert "roi_Thalamus_right_flip_fod.mif" in found_names
    assert "roi_Striatum_right_flip_fod.mif" in out["missing"]


def test_build_temp_cst_exclusion_tasks_builds_only_existing_sources(tmp_path: Path):
    cst_dir = tmp_path / "cst"
    out_dir = tmp_path / "out"
    out_dir.mkdir(parents=True)
    cst_dir.mkdir(parents=True)
    left_src = cst_dir / "tckgen_CST_left_to_graft_left.tck"
    left_src.touch()

    roi_a = tmp_path / "a.mif"
    roi_b = tmp_path / "b.mif"
    roi_a.touch()
    roi_b.touch()

    out = build_temp_cst_exclusion_tasks(
        cst_dir=cst_dir,
        out_dir=out_dir,
        output_prefix="Apallial17_CST",
        modules="module load mrtrix3",
        exclude_rois=[roi_a, roi_b],
    )

    assert len(out["tasks"]) == 1
    assert out["tasks"][0]["hemi"] == "left"
    assert "-exclude" in out["exclude_flags"]
    assert "tckgen_CST_right_to_graft_right.tck" in out["missing_sources"]


def test_collect_cst_cleaned_counts_ignores_ctrl(tmp_path: Path):
    s_ap = _subject(tmp_path, sid="sub-ap", group="AP")
    s_ctrl = _subject(tmp_path, sid="sub-ctrl", group="CTRL")
    ap_dir = s_ap["wb_dir"] / "passthrough_cst"
    ctrl_dir = s_ctrl["wb_dir"] / "passthrough_cst"
    ap_dir.mkdir(parents=True)
    ctrl_dir.mkdir(parents=True)
    (ap_dir / "a_cleaned.tck").touch()
    (ap_dir / "b_cleaned.tck").touch()
    (ctrl_dir / "c_cleaned.tck").touch()

    rows = collect_cst_cleaned_counts([s_ap, s_ctrl])
    assert len(rows) == 1
    assert rows[0]["subject"]["id"] == "sub-ap"
    assert rows[0]["n_cleaned"] == 2


def test_collect_cst_raw_tract_counts_counts_normal_and_shrink(tmp_path: Path):
    """Test that raw tract counts separate normal and shrink variants correctly."""
    s_ap = _subject(tmp_path, sid="sub-ap", group="AP")
    s_ctrl = _subject(tmp_path, sid="sub-ctrl", group="CTRL")
    ap_dir = s_ap["wb_dir"] / "passthrough_cst"
    ctrl_dir = s_ctrl["wb_dir"] / "passthrough_cst"
    ap_dir.mkdir(parents=True)
    ctrl_dir.mkdir(parents=True)

    # Create normal tract files
    (ap_dir / "normal_1.tck").touch()
    (ap_dir / "normal_2.tck").touch()
    # Create shrink variants
    (ap_dir / "shrink_1_shrink.tck").touch()
    (ap_dir / "shrink_2_shrink.tck").touch()
    # Create cleaned (should be ignored)
    (ap_dir / "cleaned_1_cleaned.tck").touch()
    (ap_dir / "cleaned_2_shrink_cleaned.tck").touch()
    # CTRL baseline
    (ctrl_dir / "ctrl_1.tck").touch()
    (ctrl_dir / "ctrl_1_shrink.tck").touch()

    rows = collect_cst_raw_tract_counts([s_ap, s_ctrl])
    assert len(rows) == 1
    assert rows[0]["subject"]["id"] == "sub-ap"
    assert rows[0]["n_normal"] == 2
    assert rows[0]["n_shrink"] == 2