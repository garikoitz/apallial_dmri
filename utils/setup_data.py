#!/usr/bin/env python
"""
setup_data.py — Install MBRAIN test data from component ZIPs on OSF.io

Components:
  - local: atlas + manual ROI files
  - pre: BIDS + pre-mask mrtrix_newmask files
  - post: remaining mrtrix_newmask outputs

Usage:
    python utils/setup_data.py
    python utils/setup_data.py --components local pre post
    python utils/setup_data.py --components local
    python utils/setup_data.py --components pre
    python utils/setup_data.py --components post
    python utils/setup_data.py --check-only

Environment variables:
    MBRAIN_URL_LOCAL: URL for the local ZIP
    MBRAIN_URL_PRE: URL for the pre-mask ZIP
    MBRAIN_URL_POST: URL for the post-mask ZIP
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path

REPO_DIR = Path(__file__).parent.parent.absolute()
TEMP_DIR = Path.home() / ".mbrain_temp"

DEFAULT_COMPONENT_URLS = {
    "local": "https://osf.io/download/6a4715f2dea535d897daa197/",
    "pre_bids": "https://osf.io/download/6a4716dfed812e7fa2daa33b/",
    "pre_preproc": "https://osf.io/download/6a4716e394a8ab0b00df9688/",
    "post": "https://osf.io/download/6a4716f594a8ab0b00df968a/",
}

MARKERS = {
    "local": REPO_DIR / ".mbrain_component_local",
    "pre_bids": REPO_DIR / ".mbrain_component_pre_bids",
    "pre_preproc": REPO_DIR / ".mbrain_component_pre_preproc",
    "post": REPO_DIR / ".mbrain_component_post",
}

SAMPLE_CHECKS = {
    "local": [
        REPO_DIR / "local" / "atlas" / "average_template_25.nii.gz",
        REPO_DIR / "local" / "manual_roi_Slicer",
        REPO_DIR / "local" / "manual_roi_flip",
    ],
    "pre_bids": [
        REPO_DIR / "data" / "bids" / "dataset_description.json",
        REPO_DIR / "data" / "bids" / "participants.tsv",
    ],
    "pre_preproc": [
        REPO_DIR / "data" / "derivatives" / "mrtrix_newmask",
    ],
    "post": [
        REPO_DIR / "data" / "derivatives" / "mrtrix_newmask",
    ],
}


def log(message: str, level: str = "INFO") -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    prefix = {
        "INFO": "[INFO]",
        "OK": "[OK]",
        "WARN": "[WARN]",
        "ERROR": "[ERROR]",
    }.get(level, "[LOG]")
    print(f"[{ts}] {prefix} {message}")


def component_url(component: str) -> str:
    env_name = f"MBRAIN_URL_{component.upper()}"
    return os.environ.get(env_name, DEFAULT_COMPONENT_URLS[component])


def component_installed(component: str) -> bool:
    return MARKERS[component].exists()


def sample_missing(component: str) -> list[Path]:
    missing = []
    if not MARKERS[component].exists():
        missing.append(MARKERS[component])
    for path in SAMPLE_CHECKS[component]:
        if not path.exists():
            missing.append(path)
    return missing


def download_file(url: str, destination: Path, chunk_size: int = 8192) -> bool:
    try:
        log(f"Downloading from: {url}")
        temp_path = destination.with_suffix(destination.suffix + ".tmp")

        with urllib.request.urlopen(url, timeout=300) as response:
            total_size = int(response.headers.get("content-length", 0))
            downloaded = 0
            with open(temp_path, "wb") as handle:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    handle.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(
                            f"\r  Progress: {percent:.1f}% ({downloaded / (1024 * 1024):.1f} MB / {total_size / (1024 * 1024):.1f} MB)",
                            end="",
                            flush=True,
                        )
            print()

        shutil.move(temp_path, destination)
        log(f"Downloaded to: {destination}", "OK")
        return True
    except urllib.error.URLError as exc:
        log(f"Download failed: {exc}", "ERROR")
        if temp_path.exists():
            temp_path.unlink()
        return False
    except Exception as exc:
        log(f"Unexpected error during download: {exc}", "ERROR")
        if temp_path.exists():
            temp_path.unlink()
        return False


def unzip_file(zip_path: Path, extract_to: Path) -> bool:
    try:
        log(f"Extracting to: {extract_to}")
        with zipfile.ZipFile(zip_path, "r") as archive:
            archive.extractall(extract_to)
        log("Extracted successfully", "OK")
        return True
    except zipfile.BadZipFile:
        log(f"Invalid ZIP file: {zip_path}", "ERROR")
        return False
    except Exception as exc:
        log(f"Extraction failed: {exc}", "ERROR")
        return False


def selected_components(args_components: list[str]) -> list[str]:
    return ["local", "pre_bids", "pre_preproc", "post"] if "all" in args_components else args_components


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install MBRAIN test data from OSF component ZIPs")
    parser.add_argument(
        "--components",
        nargs="+",
        choices=["local", "pre_bids", "pre_preproc", "post", "all"],
        default=["all"],
        help="Which components to install",
    )


def install_component(component: str, force: bool = False) -> bool:
    if component_installed(component) and not force:
        log(f"{component}: already installed", "OK")
        return True

    url = component_url(component)
    if "XXXXX" in url:
        log(f"OSF URL for '{component}' is not configured", "ERROR")
        log(f"Set MBRAIN_URL_{component.upper()} or edit DEFAULT_COMPONENT_URLS", "ERROR")
        return False

    TEMP_DIR.mkdir(exist_ok=True)
    zip_path = TEMP_DIR / f"mbrain_{component}.zip"

    if not download_file(url, zip_path):
        return False

    if not unzip_file(zip_path, REPO_DIR):
        return False

    # Ensure marker exists even if archive was produced without marker files.
    MARKERS[component].write_text(f"{component}\n", encoding="utf-8")
    zip_path.unlink(missing_ok=True)
    return True


def main() -> int:
    args = parse_args()
    components = selected_components(args.components)

    log("=" * 70)
    log("MBRAIN Data Setup", "INFO")
    log("=" * 70)
    log(f"Requested components: {' '.join(components)}")

    if args.check_only:
        any_missing = False
        for component in components:
            missing = sample_missing(component)
            if missing:
                any_missing = True
                log(f"{component}: missing {len(missing)} item(s)", "WARN")
                for item in missing[:10]:
                    print(f"  - {item.relative_to(REPO_DIR)}")
            else:
                log(f"{component}: already present", "OK")
        return 1 if any_missing else 0

    ok = True
    for component in components:
        if not install_component(component, force=args.force):
            ok = False

    if not ok:
        return 1

    log("Validating installed components...")
    all_ok = True
    for component in components:
        missing = sample_missing(component)
        if missing:
            all_ok = False
            log(f"{component}: validation failed", "WARN")
        else:
            log(f"{component}: validation passed", "OK")

    if all_ok:
        log("Setup complete! You can now run the notebook.", "OK")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
        
        log(f"Extracted successfully", "OK")
        return True
        
    except zipfile.BadZipFile:
        log(f"Invalid ZIP file: {zip_path}", "ERROR")
        return False
    except Exception as e:
        log(f"Extraction failed: {e}", "ERROR")
        return False


def validate_after_extract():
    """Verify critical files after extraction."""
    checks = [
        ("Allen template", REPO_DIR / "local" / "atlas" / "average_template_25.nii.gz"),
        ("Allen annotation", REPO_DIR / "local" / "atlas" / "annotation_25.nii.gz"),
        ("BIDS subjects", REPO_DIR / "data" / "bids" / "sub-129-921-R1F"),  # Sample
        ("Newmask derivatives", REPO_DIR / "data" / "derivatives" / "mrtrix_newmask" / "sub-129-921-R1F"),  # Sample
    ]
    
    all_ok = True
    for name, path in checks:
        if path.exists():
            log(f"{name}: {path.name}", "OK")
        else:
            log(f"{name}: NOT FOUND at {path}", "WARN")
            all_ok = False
    
    return all_ok


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Auto-download and unzip mbrain test data from OSF.io"
    )
    parser.add_argument("--force", action="store_true", help="Re-download even if local/ exists")
    parser.add_argument("--check-only", action="store_true", help="Check without downloading")
    parser.add_argument("--url", default=None, help="Override OSF.io URL")
    
    args = parser.parse_args()
    
    osf_url = args.url or os.environ.get("MBRAIN_DATA_URL", DEFAULT_OSF_URL)
    
    # ========================================================================
    # STEP 1: Check if data already exists
    # ========================================================================
    log("=" * 70)
    log("MBRAIN Data Setup", level="INFO")
    log("=" * 70)
    
    missing = check_data_structure()
    
    if not missing and not args.force:
        log("All data directories found — no download needed!", "OK")
        return 0
    
    if args.check_only:
        log("Check-only mode. Missing files:", "WARN")
        for m in missing:
            print(f"  - {m}")
        return 1
    
    # ========================================================================
    # STEP 2: Download from OSF
    # ========================================================================
    if osf_url == DEFAULT_OSF_URL or "XXXXX" in osf_url:
        log("ERROR: OSF URL not configured!", "ERROR")
        log("Please set MBRAIN_DATA_URL environment variable or update DEFAULT_OSF_URL in this script", "ERROR")
        log("Example: export MBRAIN_DATA_URL='https://osf.io/ABC12/download'")
        return 1
    
    temp_dir = Path.home() / ".mbrain_temp"
    temp_dir.mkdir(exist_ok=True)
    zip_path = temp_dir / ZIP_FILENAME
    
    if not download_file(osf_url, zip_path):
        return 1
    
    # ========================================================================
    # STEP 3: Extract to repo root
    # ========================================================================
    if not unzip_file(zip_path, REPO_DIR):
        return 1
    
    zip_path.unlink()  # Clean up temp ZIP
    
    # ========================================================================
    # STEP 4: Validate
    # ========================================================================
    log("Validating extracted data...")
    if validate_after_extract():
        log("=" * 70)
        log("Setup complete! You can now run the notebook.", "OK")
        log("=" * 70)
        return 0
    else:
        log("Validation found missing files — check above", "WARN")
        return 1


if __name__ == "__main__":
    sys.exit(main())
