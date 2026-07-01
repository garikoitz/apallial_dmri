# %% [markdown]
# # Stage 1: Setup and Atlas Mask Generation
#
# Scope:
# - Environment setup and config loading
# - Subject manifest construction
# - Derivatives copy to newmask root
# - Group atlas preparation and mask generation
#
# End-of-stage QC:
# - `mask_new.mif` exists for all subjects
# - No unresolved registration failures

# %% [Stage 1] Shared setup
import os, sys, subprocess
from pathlib import Path
from datetime import datetime

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
DATA_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain")
MODULE_LOAD_FS = "module load freesurfer/7.4.1-Rocky8"
MODULE_LOAD_MRTRIX = "module load mrtrix3/3.0.3"
MODULE_LOAD_MATLAB = "module load matlab/R2024A"
MODULE_LOAD_ANTS = "module load ants/latest"
MODULES = f"{MODULE_LOAD_FS} && {MODULE_LOAD_MRTRIX} && {MODULE_LOAD_ANTS}"

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import (
    DEFAULT_GROUP_ATLAS_IDS,
    build_group_atlas_info,
    build_mask_generation_args,
    build_subject_records,
    copy_atlas_mask,
    register_and_warp_mask,
    summarize_subjects,
)

cfg = load_config()
paths = cfg["paths"]
groups = get_subject_groups()

DERIV_SRC = Path(paths["derivatives_flip1"])
DERIV_BASE = DATA_BASE / "data" / "derivatives" / "mrtrix_newmask"
SKIP_IDS = {"sub-APtestmouseNmM"}

# %% [Stage 1] Build subject list
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids=SKIP_IDS)
summary = summarize_subjects(subjects)
print(f"Subjects: {summary['total']} | by_group={summary['by_group']}")

# %% [Stage 1] Run block
# Run the corresponding Stage 0 + Step 1 cells migrated from:
# notebooks/pipeline_groups_new_mask.py
#
# Suggested order:
# 1) Copy derivatives tree (Step 0b)
# 2) Build subject list (Step 0c)
# 3) Prepare group atlases (Step 1a)
# 4) Register/warp masks (Step 1b)

# %% [Stage 1 - QC Gate]
missing_mask_new = []
for s in subjects:
    cand = s["dwi_dir"] / "mask_new.mif"
    if not cand.exists():
        missing_mask_new.append(s["id"])

print("=" * 70)
print("STAGE 1 QC")
print("=" * 70)
if missing_mask_new:
    print(f"FAIL: missing mask_new.mif for {len(missing_mask_new)} subjects")
    for sid in missing_mask_new:
        print(f"  - {sid}")
else:
    print("PASS: mask_new.mif exists for all subjects")
