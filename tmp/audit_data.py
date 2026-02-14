#!/usr/bin/env python3
"""Audit all subjects' data on disk for status report."""
import os
from pathlib import Path

base = Path("/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/"
            "Shared drives/DTI Collab/mbrain_data/data/derivatives/mrtrix_preproc_flip1")

subjects = {
    "CTRL": [
        ("sub-129-921-R1F", "ses-260115121231"),
        ("sub-129-921-R1L1F", "ses-260116085433"),
        ("sub-129-1014-L1M", "ses-260115113814"),
        ("sub-138-1010-L1F", "ses-260115101441"),
    ],
    "AP": [
        ("sub-596-C-OO-NmM", "ses-240911103349"),
        ("sub-769-C-O5-R1M", "ses-240911115733"),
        ("sub-769-C-OO-NmM", "ses-240911111853"),
        ("sub-875-C-OO-NM", "ses-240911100253"),
    ],
    "APTR": [
        ("sub-AP72-862-901-NmF", "ses-251210093342"),
        ("sub-AP72-862-901-R1M", "ses-251210100043"),
        ("sub-Apallial17NmF", "ses-240816091633"),
        ("sub-Apallial20R1F", "ses-240909094232"),
    ],
}

for grp, subs in subjects.items():
    print(f"\n{'='*60}")
    print(f"  GROUP: {grp}")
    print(f"{'='*60}")
    for sid, ses in subs:
        wb = base / sid / ses / "dwi" / "wholebrain"
        pt = wb / "passthrough_tracts"
        pt_sh = wb / "passthrough_tracts_shrink"
        has_wb = (wb / "wholebrain_1000k.tck").exists()

        # Passthrough tracts
        raw = cleaned = sf = 0
        if pt.exists():
            for f in pt.glob("*.tck"):
                if "_cleaned" in f.name:
                    cleaned += 1
                elif "_SF" in f.name:
                    sf += 1
                else:
                    raw += 1

        # Passthrough shrink tracts
        sh_raw = sh_cleaned = sh_sf = 0
        if pt_sh.exists():
            for f in pt_sh.glob("*.tck"):
                if "_cleaned" in f.name:
                    sh_cleaned += 1
                elif "_SF" in f.name:
                    sh_sf += 1
                else:
                    sh_raw += 1

        # Segment tracts
        seg_dir = base / sid / ses / "dwi" / "segment_tracts"
        seg_raw = seg_cleaned = seg_sift2w = 0
        if seg_dir.exists():
            for f in seg_dir.iterdir():
                if f.name.endswith("_sift2_weights.txt"):
                    seg_sift2w += 1
                elif f.name.endswith("_cleaned.tck"):
                    seg_cleaned += 1
                elif f.name.endswith(".tck"):
                    seg_raw += 1

        # Shrunk graft segments
        shrink_seg = base / sid / ses / "dwi" / "graft_segments_fod" / "shrink"
        n_shrink_seg = 0
        if shrink_seg.exists():
            n_shrink_seg = len(list(shrink_seg.glob("*.mif")))

        # Check which passthrough targets exist (sample raw files)
        pt_targets = set()
        if pt.exists():
            for f in pt.glob("*.tck"):
                if "_cleaned" not in f.name and "_SF" not in f.name:
                    name = f.stem
                    # Extract target from name like: ..._to_TARGET_HEMI
                    if "_to_" in name:
                        parts = name.split("_to_")
                        pt_targets.add(parts[-1])

        print(f"  {sid}:")
        print(f"    wholebrain: {'YES' if has_wb else 'NO'}")
        print(f"    segment_tracts: {seg_raw} raw, {seg_cleaned} cleaned, {seg_sift2w} sift2w")
        print(f"    passthrough: {raw} raw, {cleaned} cleaned, {sf} SF  targets={sorted(pt_targets)}")
        if grp in ("AP", "APTR"):
            print(f"    shrink_segments: {n_shrink_seg} files")
            print(f"    pt_shrink: {sh_raw} raw, {sh_cleaned} cleaned, {sh_sf} SF")

# Also list some sample filenames from one subject per group
print(f"\n{'='*60}")
print("  SAMPLE FILENAMES")
print(f"{'='*60}")

# CTRL sample
ctrl_pt = base / "sub-129-921-R1F" / "ses-260115121231" / "dwi" / "wholebrain" / "passthrough_tracts"
if ctrl_pt.exists():
    print("\n  CTRL passthrough raw files:")
    for f in sorted(ctrl_pt.glob("*.tck")):
        if "_cleaned" not in f.name and "_SF" not in f.name:
            print(f"    {f.name}")

# AP sample - check for cerebellum with include
ap_pt = base / "sub-596-C-OO-NmM" / "ses-240911103349" / "dwi" / "wholebrain" / "passthrough_tracts"
if ap_pt.exists():
    print("\n  AP passthrough raw files (first 20):")
    count = 0
    for f in sorted(ap_pt.glob("*.tck")):
        if "_cleaned" not in f.name and "_SF" not in f.name:
            print(f"    {f.name}")
            count += 1
            if count >= 20:
                print("    ...")
                break

# AP shrink sample
ap_sh = base / "sub-596-C-OO-NmM" / "ses-240911103349" / "dwi" / "wholebrain" / "passthrough_tracts_shrink"
if ap_sh.exists():
    print("\n  AP shrink raw files (first 20):")
    count = 0
    for f in sorted(ap_sh.glob("*.tck")):
        if "_cleaned" not in f.name and "_SF" not in f.name:
            print(f"    {f.name}")
            count += 1
            if count >= 20:
                print("    ...")
                break

# Check cerebellum files specifically
print(f"\n{'='*60}")
print("  CEREBELLUM TRACTS CHECK")
print(f"{'='*60}")
for grp, subs in subjects.items():
    sid, ses = subs[0]  # First subject per group
    pt = base / sid / ses / "dwi" / "wholebrain" / "passthrough_tracts"
    if pt.exists():
        cereb = sorted([f.name for f in pt.glob("*cerebellum*") if "_cleaned" not in f.name and "_SF" not in f.name])
        print(f"\n  {grp} ({sid}) cerebellum raw tracts: {len(cereb)}")
        for c in cereb[:10]:
            print(f"    {c}")

# Check which raw passthrough tracts DON'T have cleaned versions
print(f"\n{'='*60}")
print("  UNCLEANED PASSTHROUGH TRACTS")
print(f"{'='*60}")
for grp, subs in subjects.items():
    for sid, ses in subs:
        pt = base / sid / ses / "dwi" / "wholebrain" / "passthrough_tracts"
        if not pt.exists():
            continue
        uncleaned = []
        for f in sorted(pt.glob("*.tck")):
            if "_cleaned" not in f.name and "_SF" not in f.name:
                cleaned = pt / f.name.replace(".tck", "_cleaned.tck")
                if not cleaned.exists():
                    uncleaned.append(f.name)
        if uncleaned:
            print(f"\n  [{grp}] {sid}: {len(uncleaned)} uncleaned passthrough")
            for u in uncleaned[:5]:
                print(f"    {u}")
            if len(uncleaned) > 5:
                print(f"    ... and {len(uncleaned)-5} more")

        # Also check shrink
        pt_sh = base / sid / ses / "dwi" / "wholebrain" / "passthrough_tracts_shrink"
        if pt_sh.exists():
            uncleaned_sh = []
            for f in sorted(pt_sh.glob("*.tck")):
                if "_cleaned" not in f.name and "_SF" not in f.name:
                    cleaned = pt_sh / f.name.replace(".tck", "_cleaned.tck")
                    if not cleaned.exists():
                        uncleaned_sh.append(f.name)
            if uncleaned_sh:
                print(f"  [{grp}] {sid}: {len(uncleaned_sh)} uncleaned shrink")
                for u in uncleaned_sh[:5]:
                    print(f"    {u}")
                if len(uncleaned_sh) > 5:
                    print(f"    ... and {len(uncleaned_sh)-5} more")
