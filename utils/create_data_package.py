#!/usr/bin/env python
"""Create three MBRAIN component ZIP archives for OSF upload."""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import shutil
import sys
import zipfile
from pathlib import Path

DEFAULT_DATA_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain")
DEFAULT_OUTPUT_DIR = Path("/tmp/mbrain_component_zips")


def log(msg: str) -> None:
    print(msg)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_file(src: Path, dst_root: Path, src_root: Path) -> None:
    rel = src.relative_to(src_root)
    dst = dst_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def matches_any(name: str, patterns: list[str]) -> bool:
    return any(fnmatch.fnmatch(name, p) for p in patterns)


def build_local(stage: Path, data_base: Path) -> None:
    (stage / "local").mkdir(parents=True, exist_ok=True)
    for name in ["atlas", "manual_roi_Slicer", "manual_roi_flip"]:
        shutil.copytree(data_base / "local" / name, stage / "local" / name, dirs_exist_ok=True)
    (stage / ".mbrain_component_local").write_text("local\n", encoding="utf-8")


def build_pre_bids(stage: Path, data_base: Path) -> None:
    shutil.copytree(data_base / "data" / "bids", stage / "data" / "bids", dirs_exist_ok=True)
    (stage / ".mbrain_component_pre_bids").write_text("pre_bids\n", encoding="utf-8")

def build_pre_preproc(stage: Path, data_base: Path) -> None:
    src_root = data_base / "data" / "derivatives" / "mrtrix_newmask"
    dst_root = stage / "data" / "derivatives" / "mrtrix_newmask"
    include = ["t2w_highres_flip*", "dwi*", "mask*", "*.bval", "*.bvec", "*.json"]
    for p in src_root.rglob("*"):
        if p.is_file() and matches_any(p.name, include):
            copy_file(p, dst_root, src_root)
    (stage / ".mbrain_component_pre_preproc").write_text("pre_preproc\n", encoding="utf-8")


def build_post(stage: Path, data_base: Path, exclude_tck: bool = True) -> None:
    src_root = data_base / "data" / "derivatives" / "mrtrix_newmask"
    dst_root = stage / "data" / "derivatives" / "mrtrix_newmask"
    exclude = ["t2w_highres_flip*", "dwi*", "mask*", "*.bval", "*.bvec"]
    if exclude_tck:
        exclude.append("*.tck")
    for p in src_root.rglob("*"):
        if p.is_file() and not matches_any(p.name, exclude):
            copy_file(p, dst_root, src_root)
    (stage / ".mbrain_component_post").write_text("post\n", encoding="utf-8")


def make_zip(stage: Path, zip_path: Path) -> None:
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for p in stage.rglob("*"):
            zf.write(p, p.relative_to(stage))


def verify_required(data_base: Path) -> None:
    required = [
        data_base / "local" / "atlas" / "average_template_25.nii.gz",
        data_base / "local" / "atlas" / "annotation_25.nii.gz",
        data_base / "local" / "manual_roi_Slicer",
        data_base / "local" / "manual_roi_flip",
        data_base / "data" / "bids",
        data_base / "data" / "derivatives" / "mrtrix_newmask",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required paths:\n" + "\n".join(str(p) for p in missing))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create MBRAIN OSF component packages")
    p.add_argument("--data-base", default=str(DEFAULT_DATA_BASE))
    p.add_argument("--stage-root", default=None, help="Defaults to a folder inside --output-dir")
    p.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    p.add_argument("--verify-only", action="store_true")
    p.add_argument("--include-tck", action="store_true", help="Include .tck files (excluded by default)")
    p.add_argument(
        "--components",
        nargs="+",
        choices=["local", "pre", "post", "all"],
        default=["all"],
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    data_base = Path(args.data_base)
    output_dir = Path(args.output_dir)
    
    if args.stage_root is not None:
        stage_root = Path(args.stage_root)
    else:
        stage_root = output_dir / ".mbrain_package_staging"

    comps = ["local", "pre_bids", "pre_preproc", "post"] if "all" in args.components else args.components

    log("MBRAIN component packaging")
    log(f"  Data base: {data_base}")
    log(f"  Components: {' '.join(comps)}")

    verify_required(data_base)
    if args.verify_only:
        log("Verification passed")
        return 0

    for comp in comps:
        stage = stage_root / comp
        if stage.exists():
            shutil.rmtree(stage)
        stage.mkdir(parents=True, exist_ok=True)

        if comp == "local":
            build_local(stage, data_base)
            zip_path = output_dir / "mbrain_local.zip"
        elif comp == "pre_bids":
            build_pre_bids(stage, data_base)
            zip_path = output_dir / "mbrain_pre_bids.zip"
        elif comp == "pre_preproc":
            build_pre_preproc(stage, data_base)
            zip_path = output_dir / "mbrain_pre_preproc.zip"
        else:
            build_post(stage, data_base, exclude_tck=not args.include_tck)
            zip_path = output_dir / "mbrain_post.zip"

        make_zip(stage, zip_path)
        checksum = sha256_file(zip_path)
        zip_path.with_suffix(zip_path.suffix + ".sha256").write_text(
            f"{checksum}  {zip_path.name}\n", encoding="utf-8"
        )
        size_mb = zip_path.stat().st_size / (1024 * 1024)
        log(f"  {zip_path.name}: {size_mb:.1f} MB")

    log(f"Done. Upload ZIPs from {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
