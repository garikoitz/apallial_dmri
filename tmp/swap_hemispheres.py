#!/usr/bin/env python3
"""Atomically swap _left ↔ _right in filenames for AP and APTR subjects.

Strategy (3-pass to avoid collisions):
  Pass 1: _left  → _SWAPTMP_
  Pass 2: _right → _left
  Pass 3: _SWAPTMP_ → _right

Run with DRY_RUN=True first to verify, then set to False.
"""
import os
from pathlib import Path

DRY_RUN = False   # Executing after dry run verified

DERIV = Path("data/derivatives/mrtrix_preproc_flip1")

# Subjects to process — sub-769 only (Apallial17 and 20 already done)
SUBJECTS = [
    "sub-769-C-O5-R1M/ses-240911115733",
]

# Passthrough tract directories — need left↔right filename swap
PASSTHROUGH_SUBDIRS = [
    "dwi/wholebrain/passthrough_tracts",
    "dwi/wholebrain/passthrough_tracts_shrink",
]

# Also rename .trk files inside trk/ subfolder (Freesurfer viz copies)
TRK_SUBDIRS = [
    "dwi/wholebrain/passthrough_tracts/trk",
]


def process_dir(d: Path, dry_run: bool) -> int:
    """Rename files in d using 3 separate passes to avoid collisions.

    Pass 1: _left  → _SWAPTMP_
    Pass 2: _right → _left
    Pass 3: _SWAPTMP_ → _right

    Each pass collects ALL matching files first, then renames them.
    This prevents any overwrite collision.
    """
    if not d.exists():
        return 0
    count = 0

    # Pass 1: _left → _SWAPTMP_
    targets = [(f, d / f.name.replace("_left", "_SWAPTMP_"))
               for f in d.iterdir()
               if f.is_file() and "_left" in f.name]
    for src, dst in sorted(targets):
        print(f"  {'DRY ' if dry_run else ''}P1: {src.name}")
        print(f"      → {dst.name}")
        if not dry_run:
            src.rename(dst)
        count += 1

    # Pass 2: _right → _left
    targets = [(f, d / f.name.replace("_right", "_left"))
               for f in d.iterdir()
               if f.is_file() and "_right" in f.name and "_SWAPTMP_" not in f.name]
    for src, dst in sorted(targets):
        print(f"  {'DRY ' if dry_run else ''}P2: {src.name}")
        print(f"      → {dst.name}")
        if not dry_run:
            src.rename(dst)
        count += 1

    # Pass 3: _SWAPTMP_ → _right
    targets = [(f, d / f.name.replace("_SWAPTMP_", "_right"))
               for f in d.iterdir()
               if f.is_file() and "_SWAPTMP_" in f.name]
    for src, dst in sorted(targets):
        print(f"  {'DRY ' if dry_run else ''}P3: {src.name}")
        print(f"      → {dst.name}")
        if not dry_run:
            src.rename(dst)
        count += 1

    return count


total = 0
print("\n--- Renaming PASSTHROUGH tract directories ---")
for subj in SUBJECTS:
    print(f"\n{'='*60}")
    print(f"Subject: {subj}")
    for subdir in PASSTHROUGH_SUBDIRS + TRK_SUBDIRS:
        d = DERIV / subj / subdir
        n = process_dir(d, DRY_RUN)
        if n > 0:
            print(f"    [{subdir}]: {n} files")
        total += n

print(f"\n{'='*60}")
print(f"Total files {'to rename' if DRY_RUN else 'renamed'}: {total}")
if DRY_RUN:
    print("\n*** DRY RUN — no files were changed ***")
    print("Set DRY_RUN = False to execute.")
