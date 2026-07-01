from __future__ import annotations

from pathlib import Path


DEFAULT_GROUP_ATLAS_IDS = {
    "CTRL": "sub-129-921-R1F",
    "AP": "sub-596-C-OO-NmM",
    "APTR": "sub-AP72-862-901-NmF",
}


def build_group_atlas_info(subjects: list[dict], group_atlas_ids: dict[str, str]) -> dict[str, dict]:
    """Build atlas metadata and atlas-path dictionary per group."""
    atlas_info: dict[str, dict] = {}
    for group, atlas_id in group_atlas_ids.items():
        atlas_subject = next(s for s in subjects if s["id"] == atlas_id)
        atlas_info[group] = {
            "atlas_id": atlas_id,
            "atlas_sub": atlas_subject,
            "t2w": atlas_subject["t2w"],
            "mask_mif": atlas_subject["mask_new_edited"],
            "mask_nii": atlas_subject["dwi_dir"] / "mask_new_edited.nii.gz",
            "t2w_nii": atlas_subject["dwi_dir"] / "_atlas_t2w.nii.gz",
            "t2w_n4": atlas_subject["dwi_dir"] / "_atlas_t2w_n4.nii.gz",
        }
    return atlas_info


def build_mask_generation_args(subjects: list[dict], atlas_info: dict[str, dict], modules: str) -> list[tuple]:
    """Build argument tuples consumed by mask generation workers."""
    tasks: list[tuple] = []
    for subject in subjects:
        group = subject["group"]
        info = atlas_info[group]
        tasks.append(
            (
                subject["id"],
                group,
                str(subject["t2w"]),
                str(subject["dwi_dir"]),
                str(subject["mask_new"]),
                str(subject["mask_new_mif"]),
                str(info["t2w_n4"]),
                str(info["mask_nii"]),
                modules,
            )
        )
    return tasks


def register_and_warp_mask(args):
    """Register atlas T2w to target T2w and warp atlas mask to target space."""
    import glob as _glob
    import os
    import subprocess

    sid, _group, target_t2w_mif, dwi_dir_str, mask_nii_path, mask_mif_path, atlas_t2w_n4, atlas_mask_nii, modules = args

    def _run(cmd, timeout=1800):
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)

    dwi = Path(dwi_dir_str)
    target_t2w = Path(target_t2w_mif)
    mask_nii = Path(mask_nii_path)
    mask_mif = Path(mask_mif_path)

    if not target_t2w.exists():
        return (sid, "skipped", "T2w not found")
    if mask_mif.exists():
        return (sid, "skipped", "mask already exists")

    prefix = dwi / f"_tmp_{sid}_"
    tgt_nii = dwi / f"_tmp_{sid}_t2w.nii.gz"
    tgt_n4 = dwi / f"_tmp_{sid}_t2w_n4.nii.gz"
    warped_nii = dwi / f"_tmp_{sid}_mask_warped.nii.gz"
    mask_clean = dwi / f"_tmp_{sid}_mask_clean.mif"

    warp_file = Path(f"{prefix}1Warp.nii.gz")
    affine_file = Path(f"{prefix}0GenericAffine.mat")

    steps = [
        (f"{modules} && mrconvert {target_t2w} {tgt_nii} -force -quiet", "convert T2w to NIfTI"),
        (
            f"{modules} && N4BiasFieldCorrection -d 3 -i {tgt_nii} -o {tgt_n4} "
            f"-s 2 -b [50] -c [200x200x200,0.0]",
            "N4 bias correction",
        ),
        (
            f"{modules} && antsRegistrationSyN.sh -d 3 -f {tgt_n4} -m {atlas_t2w_n4} -o {prefix} -t s -n 4",
            "ANTs SyN registration",
        ),
        (
            f"{modules} && antsApplyTransforms -d 3 -i {atlas_mask_nii} -r {tgt_n4} -o {warped_nii} "
            f"-t {warp_file} -t {affine_file} -n NearestNeighbor",
            "warp atlas mask",
        ),
        (f"{modules} && mrconvert {warped_nii} {mask_mif} -force -quiet", "convert mask to MIF"),
        (
            f"{modules} && maskfilter {mask_mif} clean {mask_clean} -largest -force -quiet "
            f"&& mv {mask_clean} {mask_mif}",
            "morphological cleanup",
        ),
        (f"{modules} && mrconvert {mask_mif} {mask_nii} -force -quiet", "save NIfTI copy"),
    ]

    for cmd, step_name in steps:
        result = _run(cmd)
        if result.returncode != 0:
            for tmp in _glob.glob(str(prefix) + "*"):
                try:
                    os.remove(tmp)
                except Exception:
                    pass
            for tmp in [tgt_nii, tgt_n4, warped_nii, mask_clean]:
                try:
                    Path(tmp).unlink(missing_ok=True)
                except Exception:
                    pass
            return (sid, "failed", f"{step_name} failed: {result.stderr[:200]}")

    for tmp in _glob.glob(str(prefix) + "*"):
        try:
            os.remove(tmp)
        except Exception:
            pass
    for tmp in [tgt_nii, tgt_n4, warped_nii, mask_clean]:
        try:
            Path(tmp).unlink(missing_ok=True)
        except Exception:
            pass

    frac = ""
    stat = _run(f"{modules} && mrstats {mask_mif}")
    if stat.returncode == 0:
        lines = stat.stdout.strip().split("\n")
        if len(lines) >= 2:
            frac = f", brain fraction={lines[1].split()[1]}"

    return (sid, "success", f"done{frac}")


def copy_atlas_mask(args):
    """Copy mask_new_edited to mask_new for atlas subjects."""
    import shutil
    import subprocess

    sid, _group, _target_t2w_mif, _dwi_dir_str, mask_nii_path, mask_mif_path, _atlas_t2w_n4, _atlas_mask_nii, modules = args
    mask_mif = Path(mask_mif_path)
    mask_nii = Path(mask_nii_path)
    edited = mask_mif.parent / "mask_new_edited.mif"

    if mask_mif.exists():
        return (sid, "skipped", "mask already exists")
    if not edited.exists():
        return (sid, "failed", "mask_new_edited.mif not found")

    shutil.copy2(str(edited), str(mask_mif))
    result = subprocess.run(
        f"{modules} && mrconvert {mask_mif} {mask_nii} -force -quiet",
        shell=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    if result.returncode != 0:
        return (sid, "failed", f"NIfTI conversion failed: {result.stderr[:150]}")

    return (sid, "success", "copied from mask_new_edited.mif (atlas)")
