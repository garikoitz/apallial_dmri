from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Callable


def _run_mrinfo_size(run_cmd: Callable, modules: str, image_path: Path) -> str:
    result = run_cmd(f"{modules} && mrinfo {image_path} -size -quiet", timeout=30, verbose=False)
    return result.stdout.strip()


def verify_subject_fod(subject: dict, run_cmd: Callable, modules: str) -> tuple[dict, list[str]]:
    """Compute one subject row plus issues for Step 5b verification table."""
    dwi = subject["dwi_dir"]
    sid = subject["id"]
    row = {"subject": sid, "group": subject["group"]}
    issues: list[str] = []

    fod_norm_grp = dwi / "wm_fod_norm_group.mif"
    if not fod_norm_grp.exists():
        issues.append(f"❌ {sid}: wm_fod_norm_group.mif MISSING")
        row["status"] = "MISSING"
        return row, issues

    mtime = datetime.fromtimestamp(fod_norm_grp.stat().st_mtime)
    row["mtime"] = mtime.strftime("%Y-%m-%d %H:%M")
    row["size_MB"] = round(fod_norm_grp.stat().st_size / 1e6, 1)

    fod_grp_raw = dwi / "wm_fod_group.mif"
    csf_fod_norm = dwi / "csf_fod_norm_group.mif"
    if fod_grp_raw.exists() and fod_norm_grp.exists():
        sz_raw = fod_grp_raw.stat().st_size
        sz_norm = fod_norm_grp.stat().st_size
        if sz_raw == sz_norm and abs(fod_grp_raw.stat().st_mtime - fod_norm_grp.stat().st_mtime) < 2:
            row["normalised"] = "❌ RAW COPY (not normalised)"
            issues.append(f"❌ {sid}: mtnormalise completely failed — wm_fod_norm_group.mif is a raw copy")
        elif not csf_fod_norm.exists() and sz_raw != sz_norm:
            row["normalised"] = "⚠️  WM-only norm"
        else:
            row["normalised"] = "✅ normalised"
    else:
        row["normalised"] = "? (cannot check)"

    try:
        fod_grid = _run_mrinfo_size(run_cmd, modules, fod_norm_grp)
        mask_grid = _run_mrinfo_size(run_cmd, modules, dwi / "mask_anat.mif")
        row["fod_grid"] = fod_grid
        row["mask_grid"] = mask_grid
        fod_dims = fod_grid.split()[:3]
        mask_dims = mask_grid.split()[:3]
        if fod_dims != mask_dims:
            issues.append(f"❌ {sid}: GRID MISMATCH — FOD {fod_dims} vs mask {mask_dims}")
            row["grid_match"] = "❌ MISMATCH"
        else:
            row["grid_match"] = "✅"
    except Exception as exc:
        row["grid_match"] = f"? error: {exc}"

    missing_resp = [
        f
        for f in ["response_wm.txt", "response_gm.txt", "response_csf.txt"]
        if not (dwi / f).exists()
    ]
    if missing_resp:
        issues.append(f"⚠️  {sid}: missing response files: {missing_resp}")
        row["responses"] = f"⚠️  missing {len(missing_resp)}"
    else:
        row["responses"] = "✅"

    mask_anat = dwi / "mask_anat.mif"
    if mask_anat.exists():
        mask_mt = datetime.fromtimestamp(mask_anat.stat().st_mtime)
        row["mask_date"] = mask_mt.strftime("%Y-%m-%d %H:%M")
    else:
        row["mask_date"] = "MISSING"
        issues.append(f"❌ {sid}: mask_anat.mif MISSING")

    row["status"] = "OK"
    return row, issues


def collect_fod_verification(subjects: list[dict], run_cmd: Callable, modules: str) -> tuple[list[dict], list[str]]:
    """Collect all Step 5b verification rows and issues."""
    rows: list[dict] = []
    issues: list[str] = []

    for subject in subjects:
        row, row_issues = verify_subject_fod(subject, run_cmd=run_cmd, modules=modules)
        rows.append(row)
        issues.extend(row_issues)

    return rows, issues


def collect_group_response_status(derivatives_dir: Path) -> list[dict]:
    """Return status rows for required group response files."""
    group_resp_dir = Path(derivatives_dir) / "group_responses"
    out: list[dict] = []
    for name in ["group_response_wm.txt", "group_response_gm.txt", "group_response_csf.txt"]:
        path = group_resp_dir / name
        if path.exists():
            mtime = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
            out.append({"name": name, "exists": True, "mtime": mtime})
        else:
            out.append({"name": name, "exists": False, "mtime": None})
    return out
