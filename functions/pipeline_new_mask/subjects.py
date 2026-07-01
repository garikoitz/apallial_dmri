from __future__ import annotations

from pathlib import Path
from typing import Iterable


def build_subject_records(
    groups: dict,
    deriv_base: Path,
    skip_ids: Iterable[str] | None = None,
    group_order: Iterable[str] = ("CTRL", "AP", "APTR"),
) -> list[dict]:
    """Build per-subject path records used by the new-mask pipeline.

    The output schema mirrors the dictionary structure currently used in
    notebooks/pipeline_groups_new_mask.py, Step 0c.
    """
    skip = set(skip_ids or [])
    records: list[dict] = []

    for group_name in group_order:
        for subject in groups.get(group_name, []):
            sid = subject["id"]
            if sid in skip:
                continue

            session = subject["session"]
            base_dir = Path(deriv_base) / sid / session
            dwi_dir = base_dir / "dwi"
            wb_dir = dwi_dir / "wholebrain"

            records.append(
                {
                    "id": sid,
                    "session": session,
                    "group": group_name,
                    "base": base_dir,
                    "dwi_dir": dwi_dir,
                    "wb_dir": wb_dir,
                    "t2w": dwi_dir / "t2w_highres_flip.mif",
                    "mask_old": dwi_dir / "mask.mif",
                    "mask_anat_old": dwi_dir / "mask_anat.mif",
                    "mask_new": dwi_dir / "mask_new.nii.gz",
                    "mask_new_mif": dwi_dir / "mask_new.mif",
                    "mask_new_edited": dwi_dir / "mask_new_edited.mif",
                    "wb_tck": wb_dir / "wholebrain_1000k.tck",
                    "sift2_weights": wb_dir / "sift2_weights.txt",
                    "sift2_mu": wb_dir / "sift2_mu.txt",
                    "pt_dir": wb_dir / "passthrough_tracts",
                    "pt_shrink_dir": wb_dir / "passthrough_tracts_shrink",
                    "fixel_dir": dwi_dir / "fixels_group",
                    "afd_voxel": dwi_dir / "afd_voxel_sum.mif",
                }
            )

    return records


def summarize_subjects(records: list[dict]) -> dict:
    """Return summary counts for a subject-record list."""
    by_group: dict[str, int] = {}
    missing_t2w = 0

    for rec in records:
        group = str(rec.get("group", "UNKNOWN"))
        by_group[group] = by_group.get(group, 0) + 1
        t2w = rec.get("t2w")
        if isinstance(t2w, Path) and not t2w.exists():
            missing_t2w += 1

    return {
        "total": len(records),
        "by_group": by_group,
        "missing_t2w": missing_t2w,
    }
