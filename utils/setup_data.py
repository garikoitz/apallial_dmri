#!/usr/bin/env python
"""
setup_data.py — Install MBRAIN test data from component ZIPs on OSF.io

Components:
  - local: atlas + manual ROI files
  - pre_bids: BIDS dataset (dataset_description.json, participants.tsv, etc.)
  - pre_preproc: pre-mask mrtrix_newmask derivatives
  - post: remaining mrtrix_newmask outputs

Usage:
    python utils/setup_data.py
    python utils/setup_data.py --components local pre_bids pre_preproc post
    python utils/setup_data.py --components local
    python utils/setup_data.py --components pre_bids
    python utils/setup_data.py --components pre_preproc
    python utils/setup_data.py --components post
    python utils/setup_data.py --check-only

Environment variables:
    MBRAIN_URL_LOCAL: URL for the local ZIP
    MBRAIN_URL_PRE_BIDS: URL for the pre_bids ZIP
    MBRAIN_URL_PRE_PREPROC: URL for the pre_preproc ZIP
    MBRAIN_URL_POST: URL for the post ZIP
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from datetime import datetime
from pathlib import Path

REPO_DIR = Path(__file__).parent.parent.absolute()
TEMP_DIR = Path.home() / ".mbrain_temp"
MIN_PYTHON = (3, 11)

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

# Top-level directories each component extracts into. Used by the preflight
# check to detect symlinks / pre-existing content that would block extraction.
COMPONENT_TARGET_DIRS = {
    "local": [REPO_DIR / "local"],
    "pre_bids": [REPO_DIR / "data" / "bids"],
    "pre_preproc": [REPO_DIR / "data" / "derivatives" / "mrtrix_newmask"],
    "post": [REPO_DIR / "data" / "derivatives" / "mrtrix_newmask"],
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


def check_python_version() -> bool:
    if sys.version_info >= MIN_PYTHON:
        return True

    candidate_env_python = None
    conda_roots = []
    if os.environ.get("CONDA_ROOT"):
        conda_roots.append(Path(os.environ["CONDA_ROOT"]))
    if os.environ.get("CONDA_EXE"):
        conda_roots.append(Path(os.environ["CONDA_EXE"]).resolve().parent.parent)
    if os.environ.get("CONDA_PREFIX"):
        conda_roots.append(Path(os.environ["CONDA_PREFIX"]).resolve())

    for conda_root in conda_roots:
        env_python = conda_root / "envs" / "mbrain" / "bin" / "python"
        if env_python.exists():
            candidate_env_python = env_python
            break

    required = ".".join(str(part) for part in MIN_PYTHON)
    current = ".".join(str(part) for part in sys.version_info[:3])
    log(f"Python {required}+ is required; found Python {current}.", "ERROR")
    log(f"Current interpreter: {sys.executable}", "ERROR")

    if candidate_env_python and Path(sys.executable).resolve() != candidate_env_python.resolve():
        print()
        print("A Conda mbrain Python was found, but plain python resolves elsewhere.")
        print("For this shell, run setup with the environment Python directly:")
        print(f"  {candidate_env_python} utils/setup_data.py")

    print()
    print("The project environment file pins python=3.11. Recreate and activate it with one of:")
    print("  conda env remove -n mbrain")
    print("  conda env create -f mbrain_env.yml")
    print("  conda activate mbrain")
    print()
    print("or:")
    print("  micromamba env remove -n mbrain")
    print("  micromamba env create -f mbrain_env.yml")
    print("  micromamba activate mbrain")
    return False


def component_url(component: str) -> str:
    env_name = f"MBRAIN_URL_{component.upper()}"
    return os.environ.get(env_name, DEFAULT_COMPONENT_URLS[component])


def component_installed(component: str) -> bool:
    """A component is installed only when the marker AND all sample paths exist.

    Checking the marker alone is insufficient: if the user removes the extracted
    data (or a symlink breaks) the marker becomes stale, and we must re-extract.
    """
    if not MARKERS[component].exists():
        return False
    return all(path.exists() for path in SAMPLE_CHECKS[component])


def sample_missing(component: str) -> list[Path]:
    missing = []
    if not MARKERS[component].exists():
        missing.append(MARKERS[component])
    for path in SAMPLE_CHECKS[component]:
        if not path.exists():
            missing.append(path)
    return missing


def remove_if_exists(path: Path) -> None:
    """Remove a file if present without noisy cleanup errors."""
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def urllib_can_open(url: str) -> bool:
    """Return False for HTTPS URLs when this Python was built without SSL."""
    if not url.lower().startswith("https://"):
        return True
    return hasattr(urllib.request, "HTTPSHandler")


def download_with_urllib(url: str, temp_path: Path, chunk_size: int = 8192) -> bool:
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
    return True


def download_with_curl(url: str, temp_path: Path) -> bool:
    curl = shutil.which("curl")
    if not curl:
        log("Python cannot open HTTPS URLs and curl is not available.", "ERROR")
        return False

    cmd = [
        curl,
        "--fail",
        "--location",
        "--retry",
        "3",
        "--retry-delay",
        "5",
        "--output",
        str(temp_path),
        url,
    ]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        log(f"curl download failed with exit code {result.returncode}", "ERROR")
        return False
    return True


def download_file(url: str, destination: Path) -> bool:
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    remove_if_exists(temp_path)

    try:
        log(f"Downloading from: {url}")
        if urllib_can_open(url):
            download_with_urllib(url, temp_path)
        else:
            log("This Python build lacks HTTPS support; using curl instead.", "WARN")
            if not download_with_curl(url, temp_path):
                remove_if_exists(temp_path)
                return False

        shutil.move(temp_path, destination)
        log(f"Downloaded to: {destination}", "OK")
        return True
    except urllib.error.URLError as exc:
        log(f"Download failed: {exc}", "ERROR")
        remove_if_exists(temp_path)
        return False
    except Exception as exc:
        log(f"Unexpected error during download: {exc}", "ERROR")
        remove_if_exists(temp_path)
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


def _describe_blocker(path: Path) -> str:
    """Return a short human-readable description of why `path` blocks extraction."""
    if path.is_symlink():
        target = path.resolve()
        return f"symlink -> {target}"
    if path.is_dir():
        # Count entries to give the user a sense of what is there.
        try:
            entries = list(path.iterdir())
        except OSError as exc:
            return f"directory (unreadable: {exc})"
        return f"directory with {len(entries)} entr{'y' if len(entries) == 1 else 'ies'}"
    if path.exists():
        return f"file ({path.stat().st_size} bytes)"
    return "exists"


def preflight_check(components: list[str], force: bool) -> tuple[bool, list[tuple[str, Path, str]]]:
    """Detect symlinks / pre-existing content that would block extraction.

    Returns (ok, blockers) where `blockers` is a list of (component, path, reason)
    tuples. When `force` is True the check is skipped (caller has accepted risk).
    """
    if force:
        return True, []

    blockers: list[tuple[str, Path, str]] = []
    seen: set[Path] = set()
    for component in components:
        for target_dir in COMPONENT_TARGET_DIRS[component]:
            if target_dir in seen:
                continue
            seen.add(target_dir)
            if target_dir.is_symlink() or target_dir.exists():
                reason = _describe_blocker(target_dir)
                blockers.append((component, target_dir, reason))

    if not blockers:
        return True, []

    log("Pre-flight check found existing content that would block extraction:", "WARN")
    for component, path, reason in blockers:
        rel = path.relative_to(REPO_DIR)
        print(f"  [{component}] {rel}  ({reason})")
    print()
    print("Extraction will fail or overwrite these paths.")
    print("To proceed, remove them first, e.g.:")
    print(f"  rm -rf {REPO_DIR / 'local'} {REPO_DIR / 'data'}")
    print("Or re-run with --force to attempt extraction anyway (may fail on symlinks).")

    if not sys.stdin.isatty():
        log("Non-interactive shell detected — aborting. Re-run with --force to override.", "ERROR")
        return False, blockers

    try:
        answer = input("Remove the listed paths now? [y/N]: ").strip().lower()
    except EOFError:
        log("No input available — aborting.", "ERROR")
        return False, blockers

    if answer not in {"y", "yes"}:
        log("Aborted by user.", "WARN")
        return False, blockers

    for _, path, _ in blockers:
        try:
            if path.is_symlink() or path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
            log(f"Removed: {path.relative_to(REPO_DIR)}", "OK")
        except OSError as exc:
            log(f"Failed to remove {path}: {exc}", "ERROR")
            return False, blockers

    return True, []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install MBRAIN test data from OSF component ZIPs")
    parser.add_argument(
        "--components",
        nargs="+",
        choices=["local", "pre_bids", "pre_preproc", "post", "all"],
        default=["all"],
        help="Which components to install",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the component is already installed",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Check which components are installed without downloading",
    )
    return parser.parse_args()


def install_component(component: str, force: bool = False) -> bool:
    if component_installed(component) and not force:
        log(f"{component}: already installed", "OK")
        return True

    if MARKERS[component].exists() and not force:
        log(f"{component}: marker present but data missing — re-extracting", "WARN")
        MARKERS[component].unlink()

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
    if zip_path.exists():
        zip_path.unlink()
    return True


def main() -> int:
    if not check_python_version():
        return 1

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

    ok_pre, _ = preflight_check(components, force=args.force)
    if not ok_pre:
        return 1

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
