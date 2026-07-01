# %% [markdown]
# # Stage 6: Passthrough and Cleaning
#
# Scope:
# - Step 5f passthrough extraction
# - Step 5g MATLAB AFQ cleaning
# - Step 5fg CTRL-only rerun branch
#
# End-of-stage QC:
# - raw and cleaned tract counts validated by subject/group

# %% [Stage 6] Shared setup
import sys
from pathlib import Path

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import (
    build_subject_records,
    collect_ctrl_passthrough_counts,
    collect_passthrough_tract_counts,
)

cfg = load_config()
groups = get_subject_groups()
DERIV_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask")
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids={"sub-APtestmouseNmM"})
ctrl_subjects = [s for s in subjects if s["group"] == "CTRL"]

# %% [Stage 6] Run block
# Run migrated cells from monolith:
# - Step 5f and 5f-QC
# - Step 5g and 5g-QC
# - Step 5fg CTRL rerun

# %% [Stage 6 - QC Gate]
print("=" * 70)
print("STAGE 6 QC")
print("=" * 70)
print("All subjects raw passthrough counts:")
for row in collect_passthrough_tract_counts(subjects):
    s = row["subject"]
    if row["has_dir"]:
        print(f"  {s['id']:<30} ({s['group']:<4}) raw={row['n_raw']}")
    else:
        print(f"  {s['id']:<30} ({s['group']:<4}) NO passthrough dir")

print("\nCTRL raw/cleaned counts:")
for row in collect_ctrl_passthrough_counts(ctrl_subjects):
    s = row["subject"]
    if row["has_dir"]:
        print(f"  {s['id']:<30} raw={row['n_raw']} cleaned={row['n_cleaned']}")
    else:
        print(f"  {s['id']:<30} NO passthrough dir")
