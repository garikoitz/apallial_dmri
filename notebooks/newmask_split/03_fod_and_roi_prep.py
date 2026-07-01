# %% [markdown]
# # Stage 3: FOD and ROI Preparation
#
# Scope:
# - Step 5b FOD rerun
# - Step 5b-group group response + normalization
# - Step 5c ROI resampling and segmentation
# - Step 5c-QC ROI space/laterality checks
#
# End-of-stage QC:
# - `wm_fod_norm_group.mif` present for all subjects
# - ROI QC reviewed and accepted

# %% [Stage 3] Shared setup
import sys
from pathlib import Path

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import build_subject_records

cfg = load_config()
groups = get_subject_groups()
DERIV_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask")
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids={"sub-APtestmouseNmM"})

# %% [Stage 3] Run block
# Run migrated cells from monolith:
# - Step 5b
# - Step 5b-group
# - Step 5b-verify
# - Step 5c
# - Step 5c-QC

# %% [Stage 3 - QC Gate]
missing_group_fod = []
for s in subjects:
    p = s["dwi_dir"] / "wm_fod_norm_group.mif"
    if not p.exists():
        missing_group_fod.append(s["id"])

print("=" * 70)
print("STAGE 3 QC")
print("=" * 70)
if missing_group_fod:
    print(f"FAIL: missing wm_fod_norm_group.mif in {len(missing_group_fod)} subjects")
    for sid in missing_group_fod:
        print(f"  - {sid}")
else:
    print("PASS: wm_fod_norm_group.mif exists for all subjects")

print("Manual check: review ROI QC issues table before continuing.")
