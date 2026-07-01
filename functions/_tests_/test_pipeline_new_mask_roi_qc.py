from pathlib import Path

from functions.pipeline_new_mask.roi_qc import evaluate_laterality_consistency, get_roi_pairs_for_group, select_sample_tract, select_raw_cleaned_comparison_pair


def test_get_roi_pairs_for_group_ctrl():
    pairs = get_roi_pairs_for_group("CTRL")
    assert "isocortex" in pairs
    assert pairs["thalamus"][0].endswith("left.nii.gz")


def test_get_roi_pairs_for_group_ap_style():
    pairs = get_roi_pairs_for_group("AP")
    assert "Graft" in pairs
    assert pairs["Graft"][0].endswith("_left_flip_fod.mif")


def test_evaluate_laterality_consistency_ok():
    left = [("A", 90.0), ("B", 92.0)]
    right = [("A", 20.0), ("B", 25.0)]
    out = evaluate_laterality_consistency(left, right, axis_mid=50.0)
    assert out["lr_separated"] is True
    assert out["issues"] == []


def test_evaluate_laterality_consistency_detects_mixed():
    left = [("A", 90.0), ("B", 20.0)]
    right = [("A", 10.0), ("B", 15.0)]
    out = evaluate_laterality_consistency(left, right, axis_mid=50.0)
    assert out["lr_separated"] is False
    assert any("LEFT ROIs on mixed sides" in m for m in out["issues"])


def test_select_sample_tract_finds_tract_ctrl(tmp_path: Path):
    """Test sample tract selection for CTRL subject."""
    s = {
        "id": "sub-ctrl",
        "group": "CTRL",
        "base": tmp_path / "sub-ctrl" / "ses-1",
        "pt_dir": tmp_path / "sub-ctrl" / "ses-1" / "pt",
    }
    s["base"].mkdir(parents=True)
    s["pt_dir"].mkdir(parents=True)
    s["pt_dir"].joinpath("atlas").mkdir(parents=True)
    s["pt_dir"].joinpath("atlas/rois/fod_space").mkdir(parents=True)
    
    # Create sample tract and ROI
    sample_tract = s["pt_dir"] / "x_thalamus_left_y.tck"
    sample_tract.touch()
    roi = s["base"] / "atlas" / "rois" / "fod_space" / "roi_thalamus_left.nii.gz"
    roi.touch()
    
    tract, roi_overlay = select_sample_tract(s)
    assert tract == sample_tract
    assert roi_overlay == roi


def test_select_sample_tract_returns_none_when_no_dir(tmp_path: Path):
    """Test that select_sample_tract returns None when pt_dir doesn't exist."""
    s = {
        "id": "sub-01",
        "group": "AP",
        "base": tmp_path / "sub-01" / "ses-1",
        "pt_dir": tmp_path / "sub-01" / "ses-1" / "pt",
    }
    # pt_dir doesn't exist
    tract, roi = select_sample_tract(s)
    assert tract is None
    assert roi is None


def test_select_raw_cleaned_comparison_pair_finds_matching_pair(tmp_path: Path):
    """Test raw/cleaned tract pair selection."""
    s = {
        "id": "sub-01",
        "pt_dir": tmp_path / "sub-01" / "pt",
    }
    s["pt_dir"].mkdir(parents=True)
    
    # Create cleaned tract (this is the discovery point)
    cleaned = s["pt_dir"] / "tract_1_cleaned.tck"
    cleaned.touch()
    # Create corresponding raw tract
    raw = s["pt_dir"] / "tract_1.tck"
    raw.touch()
    # Create another unmatched cleaned tract
    (s["pt_dir"] / "tract_2_cleaned.tck").touch()
    
    raw_found, cleaned_found = select_raw_cleaned_comparison_pair(s)
    assert raw_found == raw
    assert cleaned_found == cleaned


def test_select_raw_cleaned_comparison_pair_returns_none_when_no_cleaned(tmp_path: Path):
    """Test that function returns None when no cleaned tracts exist."""
    s = {
        "id": "sub-01",
        "pt_dir": tmp_path / "sub-01" / "pt",
    }
    s["pt_dir"].mkdir(parents=True)
    (s["pt_dir"] / "tract_1.tck").touch()  # Only raw, no cleaned
    
    raw, cleaned = select_raw_cleaned_comparison_pair(s)
    assert raw is None
    assert cleaned is None
