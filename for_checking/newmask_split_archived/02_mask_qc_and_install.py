# %% [markdown]
# # Stage 2: Mask QC, Manual Edits, and DWI Install
#
# Scope:
# - Visual QC of generated masks
# - Manual correction loop
# - Regrid/install masks on DWI space
# - Checkpoint readiness report
#
# End-of-stage QC:
# - Installed masks match DWI grid
# - Stage 3 inputs are present

# %% [Stage 2] Shared setup
import sys
from pathlib import Path
from datetime import datetime

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import (
    build_mask_grid_report,
    build_mask_install_plan,
    build_subject_records,
    compute_pipeline_checkpoint,
)

cfg = load_config()
paths = cfg["paths"]
groups = get_subject_groups()
DERIV_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask")
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids={"sub-APtestmouseNmM"})

# %% [Stage 2] Run block
# Run the corresponding Stage 2, 3 and 4 cells migrated from:
# notebooks/pipeline_groups_new_mask.py
#
# Suggested order:
# 1) Step 2a/2b/2c visual QC
# 2) Step 3 manual correction
# 3) Step 4 install mask to DWI grid

# %% [Stage 2 - QC Gate]
rows = build_mask_grid_report(subjects)
mismatches = []
missing = []
for r in rows:
    if r.get("missing_mask"):
        missing.append(r["subject_id"])
    elif r.get("has_fod") and not r.get("match"):
        mismatches.append(r["subject_id"])

chk = compute_pipeline_checkpoint(subjects)

print("=" * 70)
print("STAGE 2 QC")
print("=" * 70)
if missing:
    print(f"FAIL: missing mask.mif in {len(missing)} subjects")
if mismatches:
    print(f"FAIL: mask/FOD grid mismatch in {len(mismatches)} subjects")
if not missing and not mismatches:
    print("PASS: installed masks are valid on DWI grid")

print(f"Next-stage readiness: fod_missing={len(chk['fod_missing'])} fod_stale={len(chk['fod_stale'])}")
