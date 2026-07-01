# %% [markdown]
# # Stage 7: Reports, Fixel, and Final QC
#
# Scope:
# - Step 5h tract metrics reports
# - Step 5i and 5j fixel analyses
# - Step 6 and 7 final QC and summaries
#
# End-of-stage QC:
# - report artifacts exist
# - key fixel outputs exist
# - final QC summary complete

# %% [Stage 7] Shared setup
import sys
from pathlib import Path
from datetime import datetime as _dt

REPO_DIR = Path("/export/home/glerma/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import build_subject_records

cfg = load_config()
groups = get_subject_groups()
DERIV_BASE = Path("/export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask")
subjects = build_subject_records(groups=groups, deriv_base=DERIV_BASE, skip_ids={"sub-APtestmouseNmM"})
report_version = f"v{_dt.now().strftime('%Y%m%d')}_newmask"

# %% [Stage 7] Run block
# Run migrated cells from monolith:
# - Step 5h report generation + preview
# - Step 5i fixel blocks
# - Step 5j square ROI fixel blocks
# - Step 6 and Step 7 final QC blocks

# %% [Stage 7 - QC Gate]
reports_dir = REPO_DIR / "reports"
required = [
    reports_dir / f"tract_metrics_{report_version}.csv",
    reports_dir / f"tract_metrics_count_{report_version}.csv",
    reports_dir / f"tract_metrics_sift2_{report_version}.csv",
    reports_dir / f"tract_metrics_density_{report_version}.csv",
    reports_dir / f"tract_metrics_report_{report_version}.html",
]

missing = [p for p in required if not p.exists()]

print("=" * 70)
print("STAGE 7 QC")
print("=" * 70)
if missing:
    print("FAIL: missing report artifacts")
    for p in missing:
        print(f"  - {p}")
else:
    print("PASS: all report artifacts exist")

print("Manual check: confirm fixel figures are generated and visually valid.")
