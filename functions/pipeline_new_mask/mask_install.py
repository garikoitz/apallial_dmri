from __future__ import annotations

import subprocess
from pathlib import Path


def select_new_mask(subject: dict) -> tuple[Path | None, str | None]:
    """Select preferred mask source: edited > atlas-based > none."""
    if subject["mask_new_edited"].exists():
        return subject["mask_new_edited"], "edited"
    if subject["mask_new_mif"].exists():
        return subject["mask_new_mif"], "atlas-based"
    return None, None


def resolve_dwi_template(dwi_dir: Path) -> Path | None:
    """Resolve DWI template path for mask regrid with fallback order."""
    primary = dwi_dir / "dwi_flip_gradcheck_denoise_bias_anat.mif"
    if primary.exists():
        return primary

    for fallback in ["wm_fod_norm.mif", "dwi_preproc.mif"]:
        candidate = dwi_dir / fallback
        if candidate.exists():
            return candidate
    return None


def build_mask_install_plan(subjects: list[dict], backup_suffix: str) -> list[dict]:
    """Build per-subject plan for mask replacement step."""
    plans: list[dict] = []
    for idx, subject in enumerate(subjects, start=1):
        label = f"[{idx}/{len(subjects)}] {subject['id']} ({subject['group']})"
        dwi = subject["dwi_dir"]
        new_mask, source = select_new_mask(subject)
        template = resolve_dwi_template(dwi)

        plans.append(
            {
                "subject_id": subject["id"],
                "label": label,
                "dwi": dwi,
                "mask_mif": dwi / "mask.mif",
                "mask_anat": dwi / "mask_anat.mif",
                "mask_resampled": dwi / "mask_new_dwi_space.mif",
                "new_mask": new_mask,
                "source": source,
                "dwi_template": template,
                "backup_mask": dwi / f"mask_old_{backup_suffix}.mif",
                "backup_mask_anat": dwi / f"mask_anat_old_{backup_suffix}.mif",
            }
        )
    return plans


def get_size_tokens(image_path: Path) -> list[str]:
    """Return mrinfo -size tokens for an image path."""
    result = subprocess.run(
        f"mrinfo {image_path} -size",
        shell=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip().split()


def build_mask_grid_report(subjects: list[dict], size_reader=get_size_tokens) -> list[dict]:
    """Build per-subject mask-vs-FOD grid comparison rows for Step 4 verify."""
    rows: list[dict] = []
    for subject in subjects:
        sid = subject["id"]
        dwi = subject["dwi_dir"]
        mask_f = dwi / "mask.mif"
        fod_f = dwi / "wm_fod_norm.mif"

        if not mask_f.exists():
            rows.append({"subject_id": sid, "missing_mask": True})
            continue

        mask_dims = " ".join(size_reader(mask_f))
        row = {
            "subject_id": sid,
            "missing_mask": False,
            "mask_dims": mask_dims,
            "has_fod": fod_f.exists(),
        }

        if fod_f.exists():
            fod_dims = " ".join(size_reader(fod_f)[:3])
            row["fod_dims"] = fod_dims
            row["match"] = mask_dims == fod_dims
        rows.append(row)
    return rows


def compute_pipeline_checkpoint(subjects: list[dict], size_reader=get_size_tokens) -> dict:
    """Compute checkpoint status buckets used in Step 4 readiness summary."""
    chk = {
        "mask_ok": [],
        "mask_bad": [],
        "fod_ok": [],
        "fod_stale": [],
        "fod_missing": [],
        "grp_ok": [],
        "grp_stale": [],
        "grp_missing": [],
        "tck_ok": [],
        "tck_stale": [],
        "tck_missing": [],
    }

    for subject in subjects:
        sid = subject["id"]
        dwi = subject["dwi_dir"]
        mask_f = dwi / "mask.mif"
        mask_dws = dwi / "mask_new_dwi_space.mif"
        fod_f = dwi / "wm_fod_norm.mif"
        grp_f = dwi / "wm_fod_norm_group.mif"
        tck_f = subject["wb_tck"]

        if mask_f.exists():
            dims = size_reader(mask_f)
            if len(dims) == 3 and dims != ["46", "134", "160"]:
                chk["mask_ok"].append(sid)
            else:
                chk["mask_bad"].append(sid)
        else:
            chk["mask_bad"].append(sid)

        if fod_f.exists() and mask_dws.exists():
            if fod_f.stat().st_mtime > mask_dws.stat().st_mtime:
                chk["fod_ok"].append(sid)
            else:
                chk["fod_stale"].append(sid)
        elif fod_f.exists():
            chk["fod_stale"].append(sid)
        else:
            chk["fod_missing"].append(sid)

        if grp_f.exists() and fod_f.exists():
            if grp_f.stat().st_mtime > fod_f.stat().st_mtime:
                chk["grp_ok"].append(sid)
            else:
                chk["grp_stale"].append(sid)
        elif grp_f.exists():
            chk["grp_stale"].append(sid)
        else:
            chk["grp_missing"].append(sid)

        if tck_f.exists() and fod_f.exists():
            if tck_f.stat().st_mtime > fod_f.stat().st_mtime:
                chk["tck_ok"].append(sid)
            else:
                chk["tck_stale"].append(sid)
        elif tck_f.exists():
            chk["tck_stale"].append(sid)
        else:
            chk["tck_missing"].append(sid)

    return chk
