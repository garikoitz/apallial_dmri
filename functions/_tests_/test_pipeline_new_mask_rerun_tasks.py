from pathlib import Path

from functions.pipeline_new_mask.rerun_tasks import (
    build_batch_process_tasks,
    build_ctrl_cleaning_tasks,
    build_ctrl_passthrough_extraction_tasks,
    build_cst_cleaning_tasks,
    build_cst_passthrough_extraction_tasks,
    build_cst_tckgen_tasks,
    build_passthrough_cleaning_tasks,
    build_passthrough_extraction_tasks,
    build_shrink_segments_command,
    build_sift2_tasks,
    build_tck_downsample_tasks,
    build_wholebrain_tractography_tasks,
    collect_passthrough_tract_counts,
    collect_ctrl_passthrough_counts,
)


def _subjects():
    return [
        {"id": "sub-01", "session": "ses-a"},
        {"id": "sub-02", "session": "ses-b"},
    ]


def test_build_batch_process_tasks_steps_and_force():
    tasks = build_batch_process_tasks(
        _subjects(),
        modules="module load mrtrix",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        steps="6,7",
        timeout=1800,
        force=True,
    )
    assert len(tasks) == 2
    assert "--steps 6,7" in tasks[0][1]
    assert "--force" in tasks[0][1]
    assert str(Path("/repo")) in tasks[0][1]


def test_build_wholebrain_tractography_tasks_content():
    tasks = build_wholebrain_tractography_tasks(
        _subjects(),
        modules="module load mrtrix",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
    )
    assert "generate_wholebrain_tractogram.py" in tasks[0][1]
    assert "--derivatives-dir /deriv" in tasks[0][1]


def test_build_sift2_tasks_content():
    tasks = build_sift2_tasks(
        _subjects(),
        modules="module load mrtrix",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
    )
    assert "run_passthrough_pipeline.py" in tasks[0][1]
    assert "--step sift2" in tasks[0][1]


def test_build_tck_downsample_tasks_paths():
    tasks = build_tck_downsample_tasks(
        _subjects(),
        modules="module load mrtrix",
        derivatives_dir=Path("/deriv"),
        number=10000,
    )
    cmd = tasks[0][1]
    assert "wholebrain_1000k.tck" in cmd
    assert "wholebrain_10k.tck" in cmd
    assert "-number 10000" in cmd


def test_build_cst_tckgen_tasks_builds_for_non_ctrl(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "session": "ses-1",
        "base": tmp_path / "sub-ap" / "ses-1",
        "dwi_dir": tmp_path / "sub-ap" / "ses-1" / "dwi",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
    }
    (subj["dwi_dir"]).mkdir(parents=True)
    (subj["wb_dir"] / "passthrough_cst").mkdir(parents=True)
    (subj["base"] / "manual_roi" / "fod_space").mkdir(parents=True)
    (subj["dwi_dir"] / "wm_fod_norm_group.mif").touch()
    (subj["dwi_dir"] / "mask_anat.mif").touch()
    (subj["wb_dir"] / "passthrough_cst" / "roi_CST_left_fod.mif").touch()
    (subj["wb_dir"] / "passthrough_cst" / "roi_CST_right_fod.mif").touch()
    (subj["base"] / "manual_roi" / "fod_space" / "roi_Graft_left_flip_fod.mif").touch()
    (subj["base"] / "manual_roi" / "fod_space" / "roi_Graft_right_flip_fod.mif").touch()

    out = build_cst_tckgen_tasks([subj], modules="module load mrtrix3")
    assert len(out["tasks"]) == 1
    assert "tckgen" in out["tasks"][0][1]
    assert "-select 10000" in out["tasks"][0][1]


def test_build_cst_tckgen_tasks_reports_missing_inputs(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "session": "ses-1",
        "base": tmp_path / "sub-ap" / "ses-1",
        "dwi_dir": tmp_path / "sub-ap" / "ses-1" / "dwi",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
    }
    (subj["dwi_dir"]).mkdir(parents=True)
    (subj["wb_dir"] / "passthrough_cst").mkdir(parents=True)

    out = build_cst_tckgen_tasks([subj], modules="module load mrtrix3")
    assert out["tasks"] == []
    assert any("no FOD image" in msg for msg in out["skips"])


def test_build_cst_cleaning_tasks_adds_subject_needing_clean(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "base": tmp_path / "sub-ap" / "ses-1",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
    }
    cst_dir = subj["wb_dir"] / "passthrough_cst"
    cst_dir.mkdir(parents=True)
    (cst_dir / "graft_left_to_CST_left.tck").touch()

    out = build_cst_cleaning_tasks(
        [subj],
        modules="module load mrtrix3",
        module_load_matlab="module load matlab",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        timeout=None,
    )
    assert len(out["tasks"]) == 1
    assert "--passthrough-cst" in out["tasks"][0][1]
    assert "clean_tracts.py" in out["tasks"][0][1]


def test_build_cst_cleaning_tasks_skips_when_all_cleaned(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "base": tmp_path / "sub-ap" / "ses-1",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
    }
    cst_dir = subj["wb_dir"] / "passthrough_cst"
    cst_dir.mkdir(parents=True)
    raw = cst_dir / "graft_left_to_CST_left.tck"
    raw.touch()
    (cst_dir / "graft_left_to_CST_left_cleaned.tck").touch()

    out = build_cst_cleaning_tasks(
        [subj],
        modules="module load mrtrix3",
        module_load_matlab="module load matlab",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        timeout=None,
    )
    assert out["tasks"] == []
    assert any("already cleaned" in msg for msg in out["skips"])


def test_build_cst_passthrough_extraction_tasks_builds_normal_and_shrink(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "base": tmp_path / "sub-ap" / "ses-1",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
        "wb_tck": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain" / "wholebrain_1000k.tck",
    }
    cst_dir = subj["wb_dir"] / "passthrough_cst"
    graft_dir = subj["base"] / "manual_roi" / "fod_space"
    (graft_dir / "shrink").mkdir(parents=True)
    cst_dir.mkdir(parents=True)
    subj["wb_tck"].touch()
    (cst_dir / "roi_CST_left_fod.mif").touch()
    (cst_dir / "roi_CST_right_fod.mif").touch()
    (graft_dir / "roi_Graft_left_flip_fod.mif").touch()
    (graft_dir / "roi_Graft_right_flip_fod.mif").touch()
    (graft_dir / "shrink" / "roi_Graft_left_flip_fod_shrink.mif").touch()
    (graft_dir / "shrink" / "roi_Graft_right_flip_fod_shrink.mif").touch()

    out = build_cst_passthrough_extraction_tasks([subj], modules="module load mrtrix3")
    assert len(out["tasks"]) == 1
    cmd = out["tasks"][0][1]
    assert cmd.count("tckedit") == 4
    assert "_shrink.tck" in cmd


def test_build_cst_passthrough_extraction_tasks_reports_missing_tractogram(tmp_path: Path):
    subj = {
        "id": "sub-ap",
        "group": "AP",
        "base": tmp_path / "sub-ap" / "ses-1",
        "wb_dir": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain",
        "wb_tck": tmp_path / "sub-ap" / "ses-1" / "dwi" / "wholebrain" / "wholebrain_1000k.tck",
    }
    out = build_cst_passthrough_extraction_tasks([subj], modules="module load mrtrix3")
    assert out["tasks"] == []
    assert any("no whole-brain tractogram" in msg for msg in out["skips"])


def test_build_shrink_segments_command_none_when_no_ap_group():
    subjects = [{"id": "sub-01", "group": "CTRL"}]
    cmd = build_shrink_segments_command(
        subjects,
        modules="module load mrtrix3",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
    )
    assert cmd is None


def test_build_shrink_segments_command_uses_first_ap_subject():
    subjects = [
        {"id": "sub-01", "group": "CTRL"},
        {"id": "sub-ap", "group": "AP"},
    ]
    cmd = build_shrink_segments_command(
        subjects,
        modules="module load mrtrix3",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
    )
    assert cmd is not None
    assert "--subject sub-ap" in cmd
    assert "--create-shrink-segments" in cmd


def test_build_passthrough_extraction_tasks_includes_shrink_for_ap():
    subjects = [
        {"id": "sub-ctrl", "group": "CTRL"},
        {"id": "sub-ap", "group": "AP"},
    ]
    tasks = build_passthrough_extraction_tasks(
        subjects,
        modules="module load mrtrix3",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        timeout=None,
    )
    ctrl_cmd = next(cmd for sid, cmd, _ in tasks if sid == "sub-ctrl")
    ap_cmd = next(cmd for sid, cmd, _ in tasks if sid == "sub-ap")
    assert "--hemi left" in ctrl_cmd and "--hemi right" in ctrl_cmd
    assert "--shrink" not in ctrl_cmd
    assert ap_cmd.count("--shrink") == 2


def test_build_passthrough_cleaning_tasks_timeouts_and_shrink():
    subjects = [
        {"id": "sub-ctrl", "group": "CTRL"},
        {"id": "sub-ap", "group": "AP"},
    ]
    tasks = build_passthrough_cleaning_tasks(
        subjects,
        modules="module load mrtrix3",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
    )

    ctrl = next(t for t in tasks if t[0] == "sub-ctrl")
    ap = next(t for t in tasks if t[0] == "sub-ap")
    assert ctrl[2] == 3600
    assert ap[2] == 7200
    assert "--shrink" not in ctrl[1]
    assert ap[1].count("--shrink") == 1


def test_build_ctrl_passthrough_extraction_tasks_left_right_only():
    ctrl_subjects = [{"id": "sub-ctrl", "group": "CTRL"}]
    tasks = build_ctrl_passthrough_extraction_tasks(
        ctrl_subjects,
        modules="module load mrtrix3",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        timeout=None,
    )
    assert len(tasks) == 1
    cmd = tasks[0][1]
    assert "--hemi left" in cmd
    assert "--hemi right" in cmd
    assert "--shrink" not in cmd


def test_build_ctrl_cleaning_tasks_includes_matlab_module():
    ctrl_subjects = [{"id": "sub-ctrl", "group": "CTRL"}]
    tasks = build_ctrl_cleaning_tasks(
        ctrl_subjects,
        modules="module load mrtrix3",
        module_load_matlab="module load matlab",
        repo_dir=Path("/repo"),
        derivatives_dir=Path("/deriv"),
        timeout=None,
    )
    assert len(tasks) == 1
    cmd = tasks[0][1]
    assert "module load matlab" in cmd
    assert "--passthrough" in cmd


def test_collect_passthrough_tract_counts_counts_raw_only(tmp_path: Path):
    """Test passthrough tract count collection (raw only)."""
    s1 = {
        "id": "sub-01",
        "pt_dir": tmp_path / "sub-01" / "pt",
    }
    s1["pt_dir"].mkdir(parents=True)
    (s1["pt_dir"] / "tract_1.tck").touch()
    (s1["pt_dir"] / "tract_2.tck").touch()
    (s1["pt_dir"] / "tract_1_cleaned.tck").touch()  # Should be ignored in raw count

    s2 = {
        "id": "sub-02",
        "pt_dir": tmp_path / "sub-02" / "pt",
    }
    # s2 has no pt_dir

    rows = collect_passthrough_tract_counts([s1, s2])
    assert len(rows) == 2
    assert rows[0]["subject"]["id"] == "sub-01"
    assert rows[0]["n_raw"] == 2
    assert rows[0]["has_dir"] is True
    assert rows[1]["subject"]["id"] == "sub-02"
    assert rows[1]["n_raw"] is None
    assert rows[1]["has_dir"] is False


def test_collect_ctrl_passthrough_counts_separates_raw_and_cleaned(tmp_path: Path):
    """Test CTRL passthrough count collection (raw vs cleaned)."""
    ctrl_subjects = [
        {
            "id": "sub-ctrl-01",
            "pt_dir": tmp_path / "sub-ctrl-01" / "pt",
        }
    ]
    ctrl_subjects[0]["pt_dir"].mkdir(parents=True)
    (ctrl_subjects[0]["pt_dir"] / "raw_1.tck").touch()
    (ctrl_subjects[0]["pt_dir"] / "raw_2.tck").touch()
    (ctrl_subjects[0]["pt_dir"] / "raw_1_cleaned.tck").touch()
    (ctrl_subjects[0]["pt_dir"] / "raw_2_cleaned.tck").touch()

    rows = collect_ctrl_passthrough_counts(ctrl_subjects)
    assert len(rows) == 1
    assert rows[0]["subject"]["id"] == "sub-ctrl-01"
    assert rows[0]["n_raw"] == 2
    assert rows[0]["n_cleaned"] == 2
    assert rows[0]["has_dir"] is True
