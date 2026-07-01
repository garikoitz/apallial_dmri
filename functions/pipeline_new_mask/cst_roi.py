from __future__ import annotations

from pathlib import Path
import subprocess


def resolve_cst_source_rois(subject_id: str, manual_roi_dir: Path, legacy_dir: Path) -> tuple[Path | None, Path | None]:
    """Resolve CST source ROI files, preferring subject manual_roi and falling back to legacy dir."""
    src_left = manual_roi_dir / "roi_CST_left_flip.mif"
    src_right = manual_roi_dir / "roi_CST_right_flip.mif"

    if src_left.exists() and src_right.exists():
        return src_left, src_right

    if not src_left.exists():
        for cand in [
            legacy_dir / f"roi_CST_left_{subject_id}.mif",
            legacy_dir / f"roi_CST_left_{subject_id}_CST.mif",
        ]:
            if cand.exists():
                src_left = cand
                break

    if not src_right.exists():
        for cand in [
            legacy_dir / f"roi_CST_right_{subject_id}.mif",
            legacy_dir / f"roi_CST_right_{subject_id}_CST.mif",
        ]:
            if cand.exists():
                src_right = cand
                break

    return (src_left if src_left.exists() else None), (src_right if src_right.exists() else None)


def collect_cst_status(subjects: list[dict], *, legacy_dir: Path) -> dict[str, list]:
    """Classify AP/APTR subjects by CST ROI status."""
    out: dict[str, list] = {
        "complete": [],
        "needs_regrid": [],
        "missing": [],
        "rows": [],
    }

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        manual_roi_dir = s["base"] / "manual_roi"
        cst_dir = s["wb_dir"] / "passthrough_cst"

        fod_left = cst_dir / "roi_CST_left_fod.mif"
        fod_right = cst_dir / "roi_CST_right_fod.mif"
        has_fod = fod_left.exists() and fod_right.exists()

        src_left, src_right = resolve_cst_source_rois(sid, manual_roi_dir, legacy_dir)
        has_src = src_left is not None and src_right is not None

        if has_fod:
            status_text = "✅ FOD-space ROIs exist"
            out["complete"].append(s)
        elif has_src:
            status_text = f"⚠️  Source ROIs exist ({src_left.parent.name}/), need regrid to FOD space"
            out["needs_regrid"].append((s, src_left, src_right))
        else:
            status_text = "❌ No CST ROIs — need to create & draw"
            out["missing"].append(s)

        out["rows"].append((s, status_text, src_left, src_right, fod_left, fod_right))

    return out


def build_cst_regrid_plan(subjects: list[dict], *, legacy_dir: Path) -> dict[str, list]:
    """Build regrid plan entries and missing-source reports for AP/APTR subjects."""
    plan: dict[str, list] = {
        "to_regrid": [],
        "missing_sources": [],
        "already_complete": [],
    }

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        cst_dir = s["wb_dir"] / "passthrough_cst"
        fod_left = cst_dir / "roi_CST_left_fod.mif"
        fod_right = cst_dir / "roi_CST_right_fod.mif"
        if fod_left.exists() and fod_right.exists():
            plan["already_complete"].append(s)
            continue

        src_left, src_right = resolve_cst_source_rois(sid, s["base"] / "manual_roi", legacy_dir)
        if src_left is None or src_right is None:
            plan["missing_sources"].append(s)
            continue

        plan["to_regrid"].append((s, src_left, src_right))

    return plan


def build_cst_final_qc_summary(subjects: list[dict], *, legacy_dir: Path) -> dict[str, object]:
    """Build final CST QC rows for AP/APTR subjects and aggregate pass/fail."""
    rows: list[dict] = []
    all_ok = True

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        src_left, src_right = resolve_cst_source_rois(sid, s["base"] / "manual_roi", legacy_dir)

        cst_dir = s["wb_dir"] / "passthrough_cst"
        fod_left = cst_dir / "roi_CST_left_fod.mif"
        fod_right = cst_dir / "roi_CST_right_fod.mif"

        has_src = src_left is not None and src_right is not None
        has_fod_left = fod_left.exists()
        has_fod_right = fod_right.exists()
        has_fod = has_fod_left and has_fod_right
        icon = "✅" if has_fod else ("⚠️ " if has_src else "❌")
        if not has_fod:
            all_ok = False

        rows.append(
            {
                "subject": s,
                "icon": icon,
                "has_src": has_src,
                "has_fod_left": has_fod_left,
                "has_fod_right": has_fod_right,
            }
        )

    return {"rows": rows, "all_ok": all_ok}


def collect_cst_tckgen_counts(
    subjects: list[dict],
    *,
    hemis: tuple[str, ...] = ("left", "right"),
    run_tckinfo=None,
) -> list[dict]:
    """Collect tckgen CST tract counts per subject/hemisphere for reporting."""
    if run_tckinfo is None:
        def run_tckinfo(path: Path) -> str:
            result = subprocess.run(
                f"tckinfo {path} 2>/dev/null | awk '/^ *count:/{{print $2}}'",
                shell=True,
                capture_output=True,
                text=True,
            )
            return result.stdout.strip() or "?"

    rows: list[dict] = []
    for s in subjects:
        if s["group"] == "CTRL":
            continue
        cst_dir = s["wb_dir"] / "passthrough_cst"
        for hemi in hemis:
            out_tck = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"
            if out_tck.exists():
                rows.append({"subject": s, "hemi": hemi, "count": run_tckinfo(out_tck), "exists": True})
            else:
                rows.append({"subject": s, "hemi": hemi, "count": None, "exists": False})
    return rows


def collect_temp_exclusion_rois(
    roi_dir: Path,
    *,
    regions: tuple[str, ...] = ("Striatum", "cerebellum", "Palecortex", "Thalamus"),
    hemis: tuple[str, ...] = ("left", "right"),
) -> dict[str, list]:
    """Collect present and missing exclusion ROIs for TEMP CST filtering block."""
    found: list[Path] = []
    missing: list[str] = []
    for region in regions:
        for hemi in hemis:
            roi = roi_dir / f"roi_{region}_{hemi}_flip_fod.mif"
            if roi.exists():
                found.append(roi)
            else:
                missing.append(roi.name)
    return {"found": found, "missing": missing}


def build_temp_cst_exclusion_tasks(
    *,
    cst_dir: Path,
    out_dir: Path,
    output_prefix: str,
    modules: str,
    exclude_rois: list[Path],
    hemis: tuple[str, ...] = ("left", "right"),
) -> dict[str, list]:
    """Build TEMP CST exclusion commands and missing-source diagnostics."""
    exclude_flags = " ".join(f"-exclude {p}" for p in exclude_rois)
    tasks: list[dict] = []
    missing_sources: list[str] = []

    for hemi in hemis:
        src = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"
        dst = out_dir / f"{output_prefix}_{hemi}_noTargets.tck"
        if not src.exists():
            missing_sources.append(src.name)
            continue
        cmd = f"{modules} && tckedit {src} {dst} {exclude_flags} -force -quiet"
        tasks.append({"hemi": hemi, "src": src, "dst": dst, "cmd": cmd})

    return {"exclude_flags": exclude_flags, "tasks": tasks, "missing_sources": missing_sources}


def collect_cst_cleaned_counts(subjects: list[dict]) -> list[dict]:
    """Collect cleaned CST tract counts for AP/APTR subjects."""
    rows: list[dict] = []
    for s in subjects:
        if s["group"] == "CTRL":
            continue
        cst_dir = s["wb_dir"] / "passthrough_cst"
        n_cleaned = len(list(cst_dir.glob("*_cleaned.tck")))
        rows.append({"subject": s, "n_cleaned": n_cleaned})
    return rows


def collect_cst_raw_tract_counts(subjects: list[dict]) -> list[dict]:
    """Collect raw CST tract counts (normal vs shrink) for AP/APTR subjects."""
    rows: list[dict] = []
    for s in subjects:
        if s["group"] == "CTRL":
            continue
        cst_dir = s["wb_dir"] / "passthrough_cst"
        tcks = list(cst_dir.glob("*.tck"))
        n_normal = len([t for t in tcks if "_shrink" not in t.name and "_cleaned" not in t.name])
        n_shrink = len([t for t in tcks if "_shrink" in t.name and "_cleaned" not in t.name])
        rows.append({"subject": s, "n_normal": n_normal, "n_shrink": n_shrink})
    return rows