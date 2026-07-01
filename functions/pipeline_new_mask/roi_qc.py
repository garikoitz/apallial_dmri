from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


def get_roi_pairs_for_group(group: str) -> dict[str, tuple[str, str]]:
    """Return ROI file-name pairs for one subject group."""
    if group == "CTRL":
        return {
            "isocortex": ("roi_isocortex_left.nii.gz", "roi_isocortex_right.nii.gz"),
            "thalamus": ("roi_thalamus_left.nii.gz", "roi_thalamus_right.nii.gz"),
            "striatum": ("roi_striatum_left.nii.gz", "roi_striatum_right.nii.gz"),
            "cerebellum": ("roi_cerebellum_left.nii.gz", "roi_cerebellum_right.nii.gz"),
        }
    return {
        "Graft": ("roi_Graft_left_flip_fod.mif", "roi_Graft_right_flip_fod.mif"),
        "Thalamus": ("roi_Thalamus_left_flip_fod.mif", "roi_Thalamus_right_flip_fod.mif"),
        "Striatum": ("roi_Striatum_left_flip_fod.mif", "roi_Striatum_right_flip_fod.mif"),
        "Palecortex": ("roi_Palecortex_left_flip_fod.mif", "roi_Palecortex_right_flip_fod.mif"),
        "cerebellum": ("roi_cerebellum_left_flip_fod.mif", "roi_cerebellum_right_flip_fod.mif"),
    }


def load_roi_data(roi_path: Path):
    """Load ROI voxel data from .mif or .nii.gz into a numpy array."""
    import nibabel as nib

    p = str(roi_path)
    if p.endswith(".mif"):
        with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp:
            tmp_nii = tmp.name
        try:
            subprocess.run(["mrconvert", p, tmp_nii, "-quiet", "-force"], capture_output=True, timeout=30)
            img = nib.load(tmp_nii)
            return img.get_fdata()
        finally:
            if os.path.exists(tmp_nii):
                os.unlink(tmp_nii)

    img = nib.load(p)
    return img.get_fdata()


def get_grid_info(image_path: Path) -> tuple[list[int], list[float]]:
    """Return spatial dimensions and voxel sizes using mrinfo."""
    result = subprocess.run(
        f"mrinfo {image_path} -size -spacing",
        shell=True,
        capture_output=True,
        text=True,
        timeout=10,
    )
    lines = result.stdout.strip().split("\n")
    dims = [int(x) for x in lines[0].split()][:3]
    vox = [float(x) for x in lines[1].split()][:3]
    return dims, vox


def evaluate_laterality_consistency(
    left_centroids_ax2: list[tuple[str, float]],
    right_centroids_ax2: list[tuple[str, float]],
    axis_mid: float,
) -> dict:
    """Evaluate left/right centroid consistency and separation."""
    left_sides = ["HIGH" if c > axis_mid else "LOW" for _, c in left_centroids_ax2]
    right_sides = ["HIGH" if c > axis_mid else "LOW" for _, c in right_centroids_ax2]

    left_consistent = len(set(left_sides)) == 1
    right_consistent = len(set(right_sides)) == 1
    lr_separated = left_consistent and right_consistent and left_sides and right_sides and (left_sides[0] != right_sides[0])

    issues: list[str] = []
    if not left_consistent:
        issues.append(f"LEFT ROIs on mixed sides: {list(zip([r for r, _ in left_centroids_ax2], left_sides))}")
    if not right_consistent:
        issues.append(f"RIGHT ROIs on mixed sides: {list(zip([r for r, _ in right_centroids_ax2], right_sides))}")
    if left_consistent and right_consistent and left_sides and right_sides and left_sides[0] == right_sides[0]:
        issues.append(f"L and R centroids BOTH on {left_sides[0]} side — SWAP ERROR?")

    return {
        "left_sides": left_sides,
        "right_sides": right_sides,
        "left_consistent": left_consistent,
        "right_consistent": right_consistent,
        "lr_separated": lr_separated,
        "issues": issues,
    }


def select_sample_tract(
    subject_dict: dict,
    tract_pattern: str = "*thalamus_left*.tck",
) -> tuple[Path | None, Path | None]:
    """Select first sample tract matching pattern and resolve overlay ROI.
    
    Returns (sample_tract_path, roi_overlay_path) tuple, or (None, None) if not found.
    """
    pt_d = subject_dict.get("pt_dir")
    if not pt_d or not pt_d.exists():
        return None, None
    
    sample = sorted(pt_d.glob(tract_pattern))
    sample_tract = sample[0] if sample else None
    
    if not sample_tract:
        return None, None
    
    # Determine overlay ROI based on group
    group = subject_dict.get("group", "")
    base = subject_dict.get("base")
    
    if group == "CTRL":
        roi_overlay = base / "atlas" / "rois" / "fod_space" / "roi_thalamus_left.nii.gz"
    else:
        roi_overlay = base / "manual_roi" / "fod_space" / "roi_Thalamus_left_flip_fod.mif"
    
    return sample_tract, roi_overlay if roi_overlay.exists() else None


def select_raw_cleaned_comparison_pair(subject_dict: dict) -> tuple[Path | None, Path | None]:
    """Select first cleaned tract and its corresponding raw tract for comparison.
    
    Returns (raw_tract_path, cleaned_tract_path) tuple, or (None, None) if not found.
    """
    pt_d = subject_dict.get("pt_dir")
    if not pt_d or not pt_d.exists():
        return None, None
    
    import glob
    cleaned_files = sorted(glob.glob(str(pt_d / "*_cleaned.tck")))
    if not cleaned_files:
        return None, None
    
    cleaned = Path(cleaned_files[0])
    raw = cleaned.parent / cleaned.name.replace("_cleaned.tck", ".tck")
    
    return (raw if raw.exists() else None, cleaned if cleaned.exists() else None)
