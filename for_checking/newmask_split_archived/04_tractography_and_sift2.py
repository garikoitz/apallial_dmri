# %% [markdown]
# # Stage 4: Whole-Brain Tractography and SIFT2
#
# Scope:
# - Step 5d tractography rerun
# - Step 5d-post downsample + TRK conversion
# - Step 5e SIFT2 recomputation
#
# End-of-stage QC:
# - wholebrain outputs present
# - downsample/TRK outputs present
# - SIFT2 outputs present

# %% [Stage 4] Shared setup
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

# %% [Stage 4] Run block
# Run migrated cells from monolith:
# - Step 5d
# - Step 5d-post
# - Step 5e

# %% [Stage 4 - QC Gate]
missing = {"wb_tck": [], "wb_10k_tck": [], "wb_10k_trk": [], "sift2": []}
for s in subjects:
    wb = s["wb_dir"]
    if not (wb / "wholebrain_1000k.tck").exists():
        missing["wb_tck"].append(s["id"])
    if not (wb / "wholebrain_10k.tck").exists():
        missing["wb_10k_tck"].append(s["id"])
    if not (wb / "wholebrain_10k.trk").exists():
        missing["wb_10k_trk"].append(s["id"])
    if not s["sift2_weights"].exists():
        missing["sift2"].append(s["id"])

print("=" * 70)
print("STAGE 4 QC")
print("=" * 70)
for key, ids in missing.items():
    if ids:
        print(f"FAIL: {key} missing in {len(ids)} subjects")
if not any(missing.values()):
    print("PASS: tractography/downsample/TRK/SIFT2 outputs are complete")
