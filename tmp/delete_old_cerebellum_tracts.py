#!/usr/bin/env python3
"""
Delete old cerebellum passthrough tracts (extracted WITHOUT thalamus include)
so they can be re-extracted with the new thalamus waypoint constraint.

This removes *cerebellum* tracts from:
  - passthrough_tracts/ (all subjects)
  - passthrough_tracts_shrink/ (AP/APTR subjects)

After running this, run extract_passthrough_tracts.py to regenerate them.
"""
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

total_deleted = 0

for grp, subs in subjects.items():
    for sid, ses in subs:
        # Passthrough tracts (all groups)
        pt = base / sid / ses / "dwi" / "wholebrain" / "passthrough_tracts"
        if pt.exists():
            for f in sorted(pt.glob("*cerebellum*")):
                print(f"  DELETE: {sid}/{f.name}")
                f.unlink()
                total_deleted += 1

        # Shrink passthrough tracts (AP/APTR only)
        if grp in ("AP", "APTR"):
            pt_sh = base / sid / ses / "dwi" / "wholebrain" / "passthrough_tracts_shrink"
            if pt_sh.exists():
                for f in sorted(pt_sh.glob("*cerebellum*")):
                    print(f"  DELETE shrink: {sid}/{f.name}")
                    f.unlink()
                    total_deleted += 1

print(f"\nTotal deleted: {total_deleted} files")
