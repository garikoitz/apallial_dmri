# %% [markdown]
# # Stage 5: CST Pipeline (AP/APTR)
#
# Scope:
# - CST ROI status/create/edit/regrid
# - CST extraction and targeted tckgen
# - CST cleaning
#
# End-of-stage QC:
# - AP/APTR CST ROIs complete in FOD space
# - CST raw/tckgen/cleaned outputs available

# %% [Stage 5] Shared setup
import sys
from pathlib import Path

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import (
    build_cst_final_qc_summary,
    build_subject_records,
)

cfg = load_config()
groups = get_subject_groups()
DERIV_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask")
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids={"sub-APtestmouseNmM"})
CST_LEGACY_DIR = Path("/home/glerma/glerma/PROYECTOS/mbrain/local/manual_roi_flip")

# %% [Stage 5] Run block
# Run migrated cells from monolith:
# - Step 5e-CST-check/create/edit/regrid/QC
# - Step 5e-CST-extract
# - Step 5e-CST-tckgen
# - Step 5e-CST-clean

# %% [Stage 5 - QC Gate]
qc = build_cst_final_qc_summary(subjects, legacy_dir=CST_LEGACY_DIR)
print("=" * 70)
print("STAGE 5 QC")
print("=" * 70)
for row in qc["rows"]:
    s = row["subject"]
    print(
        f"{row['icon']} {s['id']:<30} "
        f"src={'YES' if row['has_src'] else 'NO':>3} "
        f"fod_L={'YES' if row['has_fod_left'] else 'NO':>3} "
        f"fod_R={'YES' if row['has_fod_right'] else 'NO':>3}"
    )

if qc["all_ok"]:
    print("PASS: CST ROI FOD-space completeness is OK")
else:
    print("FAIL: some AP/APTR subjects are still missing CST ROI outputs")
