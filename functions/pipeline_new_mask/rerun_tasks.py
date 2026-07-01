from __future__ import annotations

from pathlib import Path


def build_batch_process_tasks(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    steps: str,
    timeout: int | None,
    force: bool = True,
) -> list[tuple[str, str, int | None]]:
    """Build task tuples for batch_process_flip1.py reruns."""
    force_flag = " --force" if force else ""
    return [
        (
            s["id"],
            f"{modules} && cd {repo_dir} && python functions/batch_process_flip1.py "
            f"--subject {s['id']} --steps {steps}{force_flag} --derivatives-dir {derivatives_dir}",
            timeout,
        )
        for s in subjects
    ]


def build_wholebrain_tractography_tasks(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int | None = None,
    force: bool = True,
) -> list[tuple[str, str, int | None]]:
    """Build task tuples for wholebrain tractography generation."""
    force_flag = " --force" if force else ""
    return [
        (
            s["id"],
            f"{modules} && cd {repo_dir} && python functions/generate_wholebrain_tractogram.py "
            f"--subject {s['id']}{force_flag} --derivatives-dir {derivatives_dir}",
            timeout,
        )
        for s in subjects
    ]


def build_sift2_tasks(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int = 1800,
) -> list[tuple[str, str, int | None]]:
    """Build task tuples for SIFT2 computation via passthrough pipeline."""
    return [
        (
            s["id"],
            f"{modules} && cd {repo_dir} && python functions/run_passthrough_pipeline.py "
            f"--subject {s['id']} --step sift2 --derivatives-dir {derivatives_dir}",
            timeout,
        )
        for s in subjects
    ]


def build_tck_downsample_tasks(
    subjects: list[dict],
    *,
    modules: str,
    derivatives_dir: Path,
    number: int = 10000,
    timeout: int | None = None,
) -> list[tuple[str, str, int | None]]:
    """Build task tuples to downsample wholebrain tractograms with tckedit."""
    return [
        (
            s["id"],
            f"{modules} && tckedit "
            f"{derivatives_dir}/{s['id']}/{s['session']}/dwi/wholebrain/wholebrain_1000k.tck "
            f"{derivatives_dir}/{s['id']}/{s['session']}/dwi/wholebrain/wholebrain_{number//1000}k.tck "
            f"-number {number} -force",
            timeout,
        )
        for s in subjects
    ]


def build_cst_tckgen_tasks(
    subjects: list[dict],
    *,
    modules: str,
    hemis: tuple[str, ...] = ("left", "right"),
    select: int = 10000,
    timeout: int | None = None,
) -> dict[str, list]:
    """Build CST tckgen task tuples and skip diagnostics for AP/APTR subjects."""
    tasks: list[tuple[str, str, int | None]] = []
    skips: list[str] = []

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        fod = s["dwi_dir"] / "wm_fod_norm_group.mif"
        mask = s["dwi_dir"] / "mask_anat.mif"
        cst_dir = s["wb_dir"] / "passthrough_cst"
        graft_fod_dir = s["base"] / "manual_roi" / "fod_space"

        if not fod.exists():
            skips.append(f"❌ {sid}: no FOD image, skipping")
            continue
        if not mask.exists():
            skips.append(f"❌ {sid}: no mask, skipping")
            continue

        cmds: list[str] = []
        for hemi in hemis:
            cst_roi = cst_dir / f"roi_CST_{hemi}_fod.mif"
            graft_roi = graft_fod_dir / f"roi_Graft_{hemi}_flip_fod.mif"
            out_tck = cst_dir / f"tckgen_CST_{hemi}_to_graft_{hemi}.tck"

            if not cst_roi.exists():
                skips.append(f"⚠️  {sid}: no CST ROI for {hemi}, skipping")
                continue
            if not graft_roi.exists():
                skips.append(f"⚠️  {sid}: no Graft ROI for {hemi}, skipping")
                continue
            if out_tck.exists():
                continue

            cmds.append(
                f"{modules} && tckgen {fod} {out_tck} "
                f"-algorithm iFOD2 "
                f"-seed_image {cst_roi} "
                f"-include {graft_roi} "
                f"-mask {mask} "
                f"-select {select} "
                f"-step 0.02 -angle 60 -cutoff 0.05 "
                f"-minlength 2 -maxlength 30 "
                f"-nthreads 4 -force -quiet"
            )

        if cmds:
            tasks.append((sid, " && ".join(cmds), timeout))
        else:
            skips.append(f"⚪ {sid}: all tckgen CST tracts already exist or missing ROIs")

    return {"tasks": tasks, "skips": skips}


def build_cst_cleaning_tasks(
    subjects: list[dict],
    *,
    modules: str,
    module_load_matlab: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int | None = None,
) -> dict[str, list]:
    """Build CST cleaning task tuples and skip diagnostics for AP/APTR subjects."""
    tasks: list[tuple[str, str, int | None]] = []
    skips: list[str] = []

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        cst_dir = s["wb_dir"] / "passthrough_cst"
        raw_tcks = [t for t in cst_dir.glob("*.tck") if "_cleaned" not in t.name and "_SF" not in t.name]
        needs_clean = [
            t for t in raw_tcks if not (t.parent / t.name.replace(".tck", "_cleaned.tck")).exists()
        ]

        if not needs_clean:
            skips.append(f"⚪ {sid}: all CST tracts already cleaned")
            continue

        base_cmd = (
            f"{modules} && {module_load_matlab} && cd {repo_dir} && python "
            f"functions/clean_tracts.py --subject {sid} --passthrough-cst "
            f"--derivatives-dir {derivatives_dir}"
        )
        tasks.append((sid, base_cmd, timeout))

    return {"tasks": tasks, "skips": skips}


def build_cst_passthrough_extraction_tasks(
    subjects: list[dict],
    *,
    modules: str,
    hemis: tuple[str, ...] = ("left", "right"),
    timeout: int | None = None,
) -> dict[str, list]:
    """Build CST passthrough extraction tasks and skip diagnostics for AP/APTR subjects."""
    tasks: list[tuple[str, str, int | None]] = []
    skips: list[str] = []

    for s in subjects:
        if s["group"] == "CTRL":
            continue

        sid = s["id"]
        wb_tck = s["wb_tck"]
        cst_dir = s["wb_dir"] / "passthrough_cst"
        graft_fod_dir = s["base"] / "manual_roi" / "fod_space"
        graft_shrink_dir = graft_fod_dir / "shrink"

        if not wb_tck.exists():
            skips.append(f"❌ {sid}: no whole-brain tractogram, skipping")
            continue

        cmds: list[str] = []
        for hemi in hemis:
            cst_roi = cst_dir / f"roi_CST_{hemi}_fod.mif"
            if not cst_roi.exists():
                skips.append(f"⚠️  {sid}: no CST ROI for {hemi}, skipping")
                continue

            graft_roi = graft_fod_dir / f"roi_Graft_{hemi}_flip_fod.mif"
            out_tck = cst_dir / f"graft_{hemi}_to_CST_{hemi}.tck"
            if graft_roi.exists() and not out_tck.exists():
                cmds.append(
                    f"{modules} && tckedit {wb_tck} {out_tck} "
                    f"-include {graft_roi} -include {cst_roi} -force -quiet"
                )
            elif not graft_roi.exists():
                skips.append(f"⚠️  {sid}: no whole graft ROI for {hemi}")

            graft_shrink_roi = graft_shrink_dir / f"roi_Graft_{hemi}_flip_fod_shrink.mif"
            out_shrink_tck = cst_dir / f"graft_{hemi}_to_CST_{hemi}_shrink.tck"
            if graft_shrink_roi.exists() and not out_shrink_tck.exists():
                cmds.append(
                    f"{modules} && tckedit {wb_tck} {out_shrink_tck} "
                    f"-include {graft_shrink_roi} -include {cst_roi} -force -quiet"
                )

        if cmds:
            tasks.append((sid, " && ".join(cmds), timeout))
        else:
            skips.append(f"⚪ {sid}: all CST tracts already exist or no graft ROIs")

    return {"tasks": tasks, "skips": skips}


def build_shrink_segments_command(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
) -> str | None:
    """Build one-time command to generate shrunk graft segments for AP/APTR cohorts."""
    first_ap = next((s for s in subjects if s["group"] in ["AP", "APTR"]), None)
    if first_ap is None:
        return None
    return (
        f"{modules} && cd {repo_dir} && python "
        f"functions/extract_passthrough_tracts.py --subject {first_ap['id']} "
        f"--create-shrink-segments --derivatives-dir {derivatives_dir}"
    )


def build_passthrough_extraction_tasks(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int | None = None,
) -> list[tuple[str, str, int | None]]:
    """Build per-subject passthrough extraction tasks (left/right and shrink for AP/APTR)."""
    tasks: list[tuple[str, str, int | None]] = []
    for s in subjects:
        base_cmd = (
            f"{modules} && cd {repo_dir} && python "
            f"functions/extract_passthrough_tracts.py --subject {s['id']}"
        )
        cmds = [
            f"{base_cmd} --hemi left --derivatives-dir {derivatives_dir}",
            f"{base_cmd} --hemi right --derivatives-dir {derivatives_dir}",
        ]
        if s["group"] in ["AP", "APTR"]:
            cmds.append(f"{base_cmd} --hemi left --shrink --derivatives-dir {derivatives_dir}")
            cmds.append(f"{base_cmd} --hemi right --shrink --derivatives-dir {derivatives_dir}")
        tasks.append((s["id"], " && ".join(cmds), timeout))
    return tasks


def build_passthrough_cleaning_tasks(
    subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout_ctrl: int = 3600,
    timeout_ap_aptr: int = 7200,
) -> list[tuple[str, str, int | None]]:
    """Build MATLAB AFQ cleaning tasks for passthrough tracts.

    AP/APTR subjects run normal+shrink cleaning in one compound command.
    """
    tasks: list[tuple[str, str, int | None]] = []
    for s in subjects:
        base_cmd = (
            f"{modules} && cd {repo_dir} && python "
            f"functions/clean_tracts.py --subject {s['id']} --passthrough "
            f"--derivatives-dir {derivatives_dir}"
        )
        if s["group"] in ["AP", "APTR"]:
            tasks.append((s["id"], f"{base_cmd} && {base_cmd} --shrink", timeout_ap_aptr))
        else:
            tasks.append((s["id"], base_cmd, timeout_ctrl))
    return tasks


def build_ctrl_passthrough_extraction_tasks(
    ctrl_subjects: list[dict],
    *,
    modules: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int | None = None,
) -> list[tuple[str, str, int | None]]:
    """Build CTRL-only passthrough extraction tasks (left + right only)."""
    tasks: list[tuple[str, str, int | None]] = []
    for s in ctrl_subjects:
        base_cmd = (
            f"{modules} && cd {repo_dir} && python "
            f"functions/extract_passthrough_tracts.py --subject {s['id']}"
        )
        compound = " && ".join(
            [
                f"{base_cmd} --hemi left --derivatives-dir {derivatives_dir}",
                f"{base_cmd} --hemi right --derivatives-dir {derivatives_dir}",
            ]
        )
        tasks.append((s["id"], compound, timeout))
    return tasks


def build_ctrl_cleaning_tasks(
    ctrl_subjects: list[dict],
    *,
    modules: str,
    module_load_matlab: str,
    repo_dir: Path,
    derivatives_dir: Path,
    timeout: int | None = None,
) -> list[tuple[str, str, int | None]]:
    """Build CTRL-only MATLAB cleaning tasks for passthrough tracts."""
    tasks: list[tuple[str, str, int | None]] = []
    for s in ctrl_subjects:
        base_cmd = (
            f"{modules} && {module_load_matlab} && cd {repo_dir} && python "
            f"functions/clean_tracts.py --subject {s['id']} --passthrough "
            f"--derivatives-dir {derivatives_dir}"
        )
        tasks.append((s["id"], base_cmd, timeout))
    return tasks


def collect_passthrough_tract_counts(subjects: list[dict]) -> list[dict]:
    """Collect passthrough tract counts (raw only) for all subjects."""
    rows: list[dict] = []
    for s in subjects:
        pt_d = s.get("pt_dir")
        if pt_d and pt_d.exists():
            tracts = list(pt_d.glob("*.tck"))
            n_raw = len([t for t in tracts if "_cleaned" not in t.name])
        else:
            n_raw = None
        rows.append({"subject": s, "n_raw": n_raw, "has_dir": pt_d and pt_d.exists() if pt_d else False})
    return rows


def collect_ctrl_passthrough_counts(ctrl_subjects: list[dict]) -> list[dict]:
    """Collect passthrough tract counts (raw vs cleaned) for CTRL subjects."""
    rows: list[dict] = []
    for s in ctrl_subjects:
        pt_d = s.get("pt_dir")
        if pt_d and pt_d.exists():
            tracts = list(pt_d.glob("*.tck"))
            n_raw = len([t for t in tracts if "_cleaned" not in t.name])
            n_cleaned = len([t for t in tracts if "_cleaned" in t.name])
        else:
            n_raw = None
            n_cleaned = None
        rows.append({
            "subject": s,
            "n_raw": n_raw,
            "n_cleaned": n_cleaned,
            "has_dir": pt_d and pt_d.exists() if pt_d else False,
        })
    return rows
