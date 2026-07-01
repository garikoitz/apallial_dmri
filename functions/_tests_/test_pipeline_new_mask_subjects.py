from pathlib import Path

from functions.pipeline_new_mask.subjects import build_subject_records, summarize_subjects


def _mock_groups():
    return {
        "CTRL": [
            {"id": "sub-ctrl-01", "session": "ses-a"},
            {"id": "sub-ctrl-02", "session": "ses-b"},
        ],
        "AP": [
            {"id": "sub-ap-01", "session": "ses-c"},
        ],
        "APTR": [
            {"id": "sub-aptr-01", "session": "ses-d"},
        ],
    }


def test_build_subject_records_respects_group_order_and_skip():
    groups = _mock_groups()
    records = build_subject_records(
        groups=groups,
        deriv_base=Path("/tmp/deriv"),
        skip_ids={"sub-ctrl-02"},
    )

    assert [r["id"] for r in records] == ["sub-ctrl-01", "sub-ap-01", "sub-aptr-01"]
    assert all("dwi_dir" in r for r in records)
    assert all("wb_tck" in r for r in records)


def test_build_subject_records_path_schema_matches_notebook_usage():
    groups = {"CTRL": [{"id": "sub-1", "session": "ses-1"}]}
    rec = build_subject_records(groups=groups, deriv_base=Path("/x"))[0]

    assert rec["base"] == Path("/x/sub-1/ses-1")
    assert rec["dwi_dir"] == Path("/x/sub-1/ses-1/dwi")
    assert rec["mask_new_mif"] == Path("/x/sub-1/ses-1/dwi/mask_new.mif")
    assert rec["sift2_weights"] == Path("/x/sub-1/ses-1/dwi/wholebrain/sift2_weights.txt")


def test_summarize_subjects_counts_groups_and_missing_t2w(tmp_path: Path):
    groups = _mock_groups()
    records = build_subject_records(groups=groups, deriv_base=tmp_path)

    # Touch one t2w file only, leaving the others missing.
    records[0]["t2w"].parent.mkdir(parents=True, exist_ok=True)
    records[0]["t2w"].touch()

    summary = summarize_subjects(records)

    assert summary["total"] == 4
    assert summary["by_group"] == {"CTRL": 2, "AP": 1, "APTR": 1}
    assert summary["missing_t2w"] == 3
