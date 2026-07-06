# %% [markdown] Apallial_dMRI Pipeline: Version with manual Brain Mask for All Groups
# # Apallial_dMRI Pipeline: New Brain Mask (Atlas-based) + Re-run Pipeline
#
# This notebook generates **new brain masks** using single-atlas segmentation
# (ANTs SyN registration) for ALL subjects across ALL groups (CTRL, AP, APTR),
# then re-runs the downstream pipeline steps that depend on the mask.
#
# **Workflow:**
# 1. Setup — load config, enumerate all subjects
# 2. Atlas-based mask — register atlas T2w+mask to each target
# 3. Visual QC — scroll through each mask overlaid on T2w
# 4. Manual Correction — fix masks for grafted mice (AP/APTR) where needed
# 5. Re-run downstream — bias correction, registration, FOD, tractography, etc.
#
# Requires modules: mrtrix3, ants, freesurfer, matlab.
#
# ---

# %% [markdown] Execution Guide + Proposed Notebook Split
# ## Execution Guide + Proposed Notebook Split
#
# This section documents each processing stage in detail and proposes a split
# into smaller notebooks, each ending with a QC gate. The goal is to reduce
# cognitive load, isolate failures, and make reruns predictable.
#
# ### Stage-by-stage detail (current monolithic flow)
#
# 1. Stage 0 - Environment, configuration, subject manifest
#    - Inputs: config1.yaml, derivatives_flip1 root, subject-group map.
#    - Outputs: `subjects` list with canonical paths, helper wrappers, tool checks.
#    - Critical assumptions:
#      - Paths resolve on lmc02 exactly as configured.
#      - Module strings are valid for the current shell environment.
#    - Common failure modes:
#      - Missing module versions, missing t2w files, wrong derivatives root.
#
# 2. Stage 1 - Atlas-driven mask generation
#    - Inputs: one atlas subject per group + target subject t2w files.
#    - Core operations: N4 correction, ANTs registration, mask warping, largest CC filter.
#    - Outputs: `mask_new.mif` and `mask_new.nii.gz` per subject.
#    - Common failure modes:
#      - Registration drift, bad transform chain, stale atlas mask paths.
#
# 3. Stage 2 and 3 - Visual QC and manual correction loop
#    - Inputs: generated masks + t2w volumes.
#    - Core operations: interactive overlay review, ROI edits in mrview.
#    - Outputs: `mask_new_edited.mif` for failed subjects.
#    - Common failure modes:
#      - Unsaved edits, wrong file edited, inconsistent pass/fail tracking.
#
# 4. Stage 4 - Mask installation on DWI grid
#    - Inputs: selected new/edited mask and DWI template.
#    - Core operations: regrid to DWI, backup old masks, install mask/mask_anat.
#    - Outputs: installed DWI-space masks and backup snapshots.
#    - Common failure modes:
#      - Grid mismatch, template missing, backup overwrite confusion.
#
# 5. Stage 5 (core downstream reruns)
#    - 5b/5b-group: response functions + FOD + group normalization.
#    - 5c: ROI regrid/segmentation and laterality checks.
#    - 5d/5e: wholebrain tractography, downsample/trk, SIFT2.
#    - 5e-CST: CST ROI creation/regrid/extraction/tckgen/cleaning.
#    - 5f/5g: passthrough extraction and AFQ cleaning, including CTRL rerun branch.
#    - 5h/5i/5j: metrics reports and fixel analyses.
#    - Common failure modes:
#      - Stale outputs mixed with fresh outputs, missing intermediate files,
#        heavy steps rerun without checkpoint gating.
#
# 6. Stage 6 and 7 - Final QC and summary views
#    - Inputs: all final masks/tracts/reports.
#    - Outputs: operator validation and publication screenshots.
#
# ### Proposed split into smaller notebooks (recommended)
#
# Notebook 1: `01_setup_and_masks.ipynb`
# - Contains: Stage 0, Stage 1.
# - Final QC gate:
#   - Confirm all subjects have `mask_new.mif`.
#   - Confirm atlas registrations completed (no failed list).
#   - Quick grid sanity for a sample of subjects.
#
# Notebook 2: `02_mask_qc_and_install.ipynb`
# - Contains: Stage 2, Stage 3, Stage 4.
# - Final QC gate:
#   - Confirm manual edits saved where required.
#   - Confirm installed masks match DWI grid for all subjects.
#   - Confirm checkpoint logic marks Stage 5 as runnable.
#
# Notebook 3: `03_fod_and_roi_prep.ipynb`
# - Contains: Stage 5b, 5b-group, 5c, 5c-QC.
# - Final QC gate:
#   - `wm_fod_norm_group.mif` present for all subjects.
#   - ROI grid/laterality checks pass or have reviewed exceptions.
#
# Notebook 4: `04_tractography_and_sift2.ipynb`
# - Contains: Stage 5d, 5d-post, 5e.
# - Final QC gate:
#   - Wholebrain tractograms present.
#   - 10k downsample + trk conversion complete.
#   - SIFT2 weights present and newer than tractograms.
#
# Notebook 5: `05_cst_pipeline.ipynb`
# - Contains: Stage 5e-CST family (check/create/edit/regrid/extract/tckgen/clean).
# - Final QC gate:
#   - AP/APTR subjects have CST ROIs in FOD space.
#   - CST extraction and cleaned outputs present with expected counts.
#
# Notebook 6: `06_passthrough_and_cleaning.ipynb`
# - Contains: Stage 5f, 5f-QC, 5g, 5g-QC, 5fg CTRL rerun branch.
# - Final QC gate:
#   - Raw and cleaned passthrough tracts present by group.
#   - Spot-check raw vs cleaned streamline counts.
#
# Notebook 7: `07_reports_and_fixel.ipynb`
# - Contains: Stage 5h, 5i, 5j, Stage 6, Stage 7.
# - Final QC gate:
#   - Metrics CSV/HTML generated.
#   - Fixel outputs and comparison figures generated.
#   - Final summary table complete with no missing critical artifacts.
#
# ### Minimal handoff contract between notebooks
#
# Each notebook should start with a short preflight cell that validates its
# required input artifacts and stop immediately if requirements are missing.
# This prevents accidental partial reruns and keeps provenance clear.
#
# Suggested convention:
# - Write one JSON checkpoint per notebook under `DERIV_BASE/logs/checkpoints/`
# - Include: timestamp, user, notebook name, pass/fail, required artifact stats.

# %% [markdown] Steps 1-4: BIDS→Newmask Pipeline (Optional)
# # Steps 1-4: BIDS→Newmask Preprocessing Pipeline
#
# **IMPORTANT: Read this section FIRST**
#
# This notebook normally assumes data is already in `mrtrix_newmask/` (preprocessed DWI, FODs, etc.).
# If you **only have BIDS raw data + atlas + manual ROI delineations**, you must run Steps 1-4 first.
#
# **If you ALREADY have `mrtrix_newmask/` populated:** Skip Steps 1-4 → go to **Step 0**
#
# **If starting from BIDS raw data:** Run Steps 1 → 4 in order, then proceed to Step 0
#
# **Data required for Steps 1-4:**
# ```
# data/bids/sub-*/ses-*/                    ← Raw BIDS input (anat + dwi)
# local/atlas/                              ← Allen CCF atlas
# local/manual_roi_Slicer/                  ← Manual ROI delineations (.seg.nrrd)
# local/manual_roi_flip/                    ← Corrected manual ROI masks (.mif)
# ```

# %% [markdown] Step 1: Verify BIDS Data & Configuration
# ## Step 1: Verify BIDS Data Exists

# %% [Step 1] Auto-download Missing Data (if needed)
# If you're running this on a fresh machine, this cell will automatically
# download the test data from OSF.io if the 'local/' directory is missing.

import os, sys, subprocess, json, yaml
from pathlib import Path
from datetime import datetime

# Resolve repository root by searching upward for config1.yaml.
_cwd = Path.cwd().resolve()
REPO_DIR = next((p for p in [_cwd, *_cwd.parents] if (p / "config1.yaml").exists()), None)
if REPO_DIR is None:
    raise RuntimeError("Could not find repository root containing config1.yaml")

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

from functions.config_loader import load_config, get_subject_groups
cfg = load_config(config_path=REPO_DIR / "config1.yaml")
groups = get_subject_groups()
paths = cfg["paths"]

DATA_BASE = Path(paths["data_base"])
BIDS_DIR = Path(paths["bids_dir"])
DERIV_NEWMASK = Path(paths["derivatives_newmask"])
LOCAL_DIR = Path(paths["local"])

# Check if data needs to be downloaded
local_dir = LOCAL_DIR
if not local_dir.exists():
    print("=" * 70)
    print(" LOCAL DATA NOT FOUND — AUTO-DOWNLOADING FROM OSF.io")
    print("=" * 70)
    print()
    
    setup_script = REPO_DIR / "utils" / "setup_data.py"
    if setup_script.exists():
        result = subprocess.run(
            f"cd {REPO_DIR} && python utils/setup_data.py",
            shell=True,
            capture_output=False
        )
        if result.returncode != 0:
            print("\n⚠️  Auto-download failed.")
            print("   Retry from repo root with setup_data.py:")
            print("     python utils/setup_data.py --check-only")
            print("     python utils/setup_data.py --components local pre_bids pre_preproc post")
            print("   You can also install only missing components, e.g.:")
            print("     python utils/setup_data.py --components local")
            print("     python utils/setup_data.py --components pre_bids")
            print("     python utils/setup_data.py --components pre_preproc")
            print("     python utils/setup_data.py --components post")
    else:
        print(f"❌ Setup script not found at {setup_script}")
        print("   Run from repository root where utils/setup_data.py exists.")
    print()
else:
    print(f"✅ Local data directory found at {local_dir}")

# %% [Step 1b] Verify BIDS Structure & Required Files

print(f"{'=' * 70}")
print(f" STEP 1: VERIFY BIDS DATA")
print(f"{'=' * 70}")
print(f"REPO_DIR:  {REPO_DIR}")
print(f"DATA_BASE: {DATA_BASE}")

print(f"\nData locations:")
print(f"  BIDS:          {BIDS_DIR}")
print(f"  Newmask:       {DERIV_NEWMASK}")
print(f"  Atlas+ROIs:    {LOCAL_DIR}")

# Check for BIDS subjects
print(f"\n{'─' * 70}")
print(f" BIDS Subjects")
print(f"{'─' * 70}")
bids_subjects = set()
if BIDS_DIR.exists():
    bids_subjects = {d.name for d in BIDS_DIR.iterdir() if d.is_dir() and d.name.startswith("sub-")}
    print(f"  Found: {len(bids_subjects)} subjects in BIDS")
else:
    print(f"  ⚠️  BIDS directory not found at {BIDS_DIR}")

# Check for required subjects from config
required = set()
for group, subjs in groups.items():
    for s in subjs:
        required.add(f"sub-{s['id']}")

in_bids = required & bids_subjects
missing_bids = required - bids_subjects

for s in sorted(in_bids):
    print(f"    ✅ {s}")
for s in sorted(missing_bids):
    in_deriv = (DERIV_NEWMASK / s).exists()
    status = "  (in mrtrix_newmask — skip preprocessing)" if in_deriv else "  ❌ MISSING"
    print(f"    ⚠️  {s}{status}")

# Check atlas and ROI files
print(f"\n{'─' * 70}")
print(f" Atlas & ROI Files")
print(f"{'─' * 70}")

atlas_files = {
    "Allen CCF Template": LOCAL_DIR / "atlas" / "average_template_25.nii.gz",
    "Allen Annotation": LOCAL_DIR / "atlas" / "annotation_25.nii.gz",
}
for name, path in atlas_files.items():
    status = "✅" if path.exists() else "❌"
    print(f"  {status} {name}")

slicer_dir = LOCAL_DIR / "manual_roi_Slicer"
n_slicer = len(list(slicer_dir.glob("*.nrrd"))) if slicer_dir.exists() else 0
print(f"  {'✅' if n_slicer > 0 else '❌'} Slicer manual ROIs: {n_slicer} files")

flip_dir = LOCAL_DIR / "manual_roi_flip"
n_flip = len(list(flip_dir.glob("*.mif"))) if flip_dir.exists() else 0
print(f"  {'✅' if n_flip > 0 else '❌'} Flipped manual ROIs: {n_flip} files")

print(f"\n{'=' * 70}")
print(f"If all items above show ✅, proceed to Step 2")
print(f"If anything shows ❌ and you need it, check your data directory structure")
print(f"{'=' * 70}\n")

# %% [markdown] Step 2: Run Batch Preprocessing
# ## Step 2: Batch Preprocessing (BIDS→Flip1)
#
# This runs the complete preprocessing pipeline for all subjects with BIDS data:
# - T2w flip1 transformation
# - DWI formatting + preprocessing (denoise, bias correct, register to T2w)
#
# **Runtime:** ~30-60 min per subject | **Parallelization:** Uses all CPU cores

# %% [Step 2] Run Batch T2w Flip + DWI Preprocessing
print(f"{'=' * 70}")
print(f" STEP 2: BATCH PREPROCESSING (T2w Flip + DWI)")
print(f"{'=' * 70}\n")

cmd_preproc = (
    f"cd {REPO_DIR} && "
    f"python functions/batch_process_flip1.py --group CTRL AP APTR --steps 2,4 --force"
)
print(f"Running:\n  {cmd_preproc}\n")
print(f"This may take 30-60 minutes per subject. Please wait...\n")

result = subprocess.run(cmd_preproc, shell=True, capture_output=False, timeout=14400)
if result.returncode == 0:
    print(f"\n✅ Preprocessing complete!")
else:
    print(f"\n❌ Preprocessing failed (exit code {result.returncode})")

# %% [markdown] Step 3: FOD Estimation + ROI Resampling
# ## Step 3: FOD Estimation + ROI Resampling
#
# Computes FODs and resamples ROIs from T2w→FOD space

# %% [Step 3] FOD + ROI Resampling
print(f"{'=' * 70}")
print(f" STEP 3: FOD ESTIMATION + ROI RESAMPLING")
print(f"{'=' * 70}\n")

cmd_fod = (
    f"cd {REPO_DIR} && "
    f"python functions/batch_process_flip1.py --group CTRL AP APTR --steps 5,6 --force"
)
print(f"Running:\n  {cmd_fod}\n")

result = subprocess.run(cmd_fod, shell=True, capture_output=False, timeout=7200)
if result.returncode == 0:
    print(f"\n✅ FOD + ROI resampling complete!")
else:
    print(f"\n❌ FOD/ROI step failed (exit code {result.returncode})")

# %% [markdown] Step 4: Segmentation + Tractography
# ## Step 4: Segmentation + Whole-Brain Tractography
#
# Creates graft/isocortex segments (FOD space) and generates 1M streamline whole-brain tractograms

# %% [Step 4] Segmentation + Tractography (+ Sync to newmask)
print(f"{'=' * 70}")
print(f" STEP 4: SEGMENTATION + TRACTOGRAPHY")
print(f"{'=' * 70}\n")

cmd_tck = (
    f"cd {REPO_DIR} && "
    f"python functions/batch_process_flip1.py --group CTRL AP APTR --steps 7,8"
)
print(f"Running:\n  {cmd_tck}\n")
print(f"tckgen can take several hours (generates 1M streamlines/subject). Please wait...\n")

result = subprocess.run(cmd_tck, shell=True, capture_output=False, timeout=43200)
if result.returncode == 0:
    print(f"\n✅ Segmentation + tractography complete!")

    # Prepare mrtrix_newmask for downstream steps (Step 0+ uses only newmask)
    deriv_flip1 = Path(paths["derivatives_flip1"])
    deriv_newmask = Path(paths["derivatives_newmask"])
    deriv_newmask.mkdir(parents=True, exist_ok=True)
    cmd_sync = f"rsync -av --progress {deriv_flip1}/ {deriv_newmask}/"
    print(f"\nSyncing flip1 → newmask for downstream notebook steps...")
    sync_res = subprocess.run(cmd_sync, shell=True, capture_output=False, timeout=3600)

    if sync_res.returncode == 0:
        print(f"\n✅ mrtrix_newmask is ready.")
        print(f"Proceed to Step 0 to start the new-mask analysis.")
    else:
        print(f"\n❌ Sync to mrtrix_newmask failed (exit code {sync_res.returncode}).")
        print(f"Please sync flip1→newmask manually before Step 0.")
else:
    print(f"\n❌ Tractography step failed (exit code {result.returncode})")

print(f"\n{'=' * 70}\n")

# %% [markdown] Step 0: Setup & Configuration
# ## Step 0: Setup & Configuration
#
# Load config, enumerate all subjects from all groups.

# %% [Step 0] Setup — Load Config & Build Subject List
import os, sys, subprocess, glob, json
from pathlib import Path
from datetime import datetime

# Resolve repository root by searching upward for config1.yaml.
_cwd = Path.cwd().resolve()
REPO_DIR = next((p for p in [_cwd, *_cwd.parents] if (p / "config1.yaml").exists()), None)
if REPO_DIR is None:
    raise RuntimeError("Could not find repository root containing config1.yaml")

MODULE_LOAD_FS = "module load freesurfer/7.4.1-Rocky8"
MODULE_LOAD_MRTRIX = "module load mrtrix3/3.0.3"
MODULE_LOAD_MATLAB = "module load matlab/R2024A"
MODULE_LOAD_ANTS = "module load ants/latest"
MODULES = f"{MODULE_LOAD_FS} && {MODULE_LOAD_MRTRIX} && {MODULE_LOAD_ANTS}"

sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

# Load config (now resolves {variable} placeholders automatically)
from functions.config_loader import load_config, get_subject_groups
from functions.pipeline_new_mask import build_subject_records, summarize_subjects
from functions.pipeline_new_mask import check_file, launch_gui, run_command, run_parallel as _run_parallel
from functions.pipeline_new_mask import (
    build_batch_process_tasks,
    build_ctrl_cleaning_tasks,
    build_ctrl_passthrough_extraction_tasks,
    build_temp_cst_exclusion_tasks,
    build_cst_final_qc_summary,
    build_cst_cleaning_tasks,
    build_cst_passthrough_extraction_tasks,
    collect_cst_cleaned_counts,
    collect_cst_tckgen_counts,
    collect_temp_exclusion_rois,
    build_cst_tckgen_tasks,
    build_passthrough_cleaning_tasks,
    build_passthrough_extraction_tasks,
    build_shrink_segments_command,
    collect_fod_verification,
    collect_group_response_status,
    build_mask_grid_report,
    DEFAULT_GROUP_ATLAS_IDS,
    build_group_atlas_info,
    build_mask_generation_args,
    build_mask_install_plan,
    build_sift2_tasks,
    build_tck_downsample_tasks,
    build_wholebrain_tractography_tasks,
    compute_pipeline_checkpoint,
    copy_atlas_mask,
    convert_tck_to_trk,
    build_cst_regrid_plan,
    collect_cst_status,
    evaluate_laterality_consistency,
    get_grid_info,
    get_roi_pairs_for_group,
    load_roi_data,
    draw_single_seg_panel,
    resolve_cst_source_rois,
    register_and_warp_mask,
)
cfg = load_config()
paths = cfg['paths']
groups = get_subject_groups()
DATA_BASE = Path(paths['data_base'])

# Config sub-sections used in downstream steps
tck = cfg['tractography']
wb  = cfg.get('wholebrain_tractography', {})
pt  = cfg.get('passthrough_extraction', {})
mat = cfg['matlab_cleaning']
sk  = cfg.get('shrink_analysis', {})

# For AP/APTR: segmentation mode (9 = 3 A-P x 3 D-V — recommended)
SEG_MODE = 9

# Destination derivatives directory (single source of truth for this notebook)
DERIV_BASE = Path(paths['derivatives_newmask'])

print(f"New mask derivatives: {DERIV_BASE}")

# %% [Step 0b] Validate mrtrix_newmask is present
# ============================================================
# This notebook now uses mrtrix_newmask directly as the derivative dataset.
# No copy from mrtrix_preproc_flip1 is performed in this workflow.
# ============================================================

print(f"{'=' * 70}")
print(f" CHECKING mrtrix_newmask")
print(f"{'=' * 70}")
print(f"  Path: {DERIV_BASE}")

if not DERIV_BASE.exists() or not any(DERIV_BASE.iterdir()):
    raise FileNotFoundError(
        f"mrtrix_newmask not found or empty: {DERIV_BASE}\n"
        "Please download the OSF package (fast-track) or prepare mrtrix_newmask before running Step 0c+."
    )

available_subs = [d for d in DERIV_BASE.iterdir() if d.is_dir() and d.name.startswith("sub-")]
print(f"\n✅ Found {len(available_subs)} subject folders in {DERIV_BASE.name}")

# %% [Step 0c] Build subject list from mrtrix_newmask
# Build flat list of all subjects with their group, session, and paths.
# All paths now point to the mrtrix_newmask directory.
# Exclude test mice — only real subjects (4 per group = 12 total).
SKIP_IDS = {'sub-APtestmouseNmM'}

subjects = build_subject_records(
    groups=groups,
    deriv_base=DERIV_BASE,
    skip_ids=SKIP_IDS,
)
subject_summary = summarize_subjects(subjects)

print(f"{'=' * 70}")
print(f" SUBJECTS TO PROCESS: {subject_summary['total']}")
print(f"{'=' * 70}")
for s in subjects:
    exists = "✅" if s['t2w'].exists() else "❌"
    print(f"  {exists} [{s['group']:>4}] {s['id']:<30} {s['session']}")
print(f"  By group: {subject_summary['by_group']}")
print(f"{'=' * 70}")
print(f"DERIV_BASE: {DERIV_BASE}")


# %% [Step 0] Helper Configuration
import multiprocessing

MAX_WORKERS = min(12, multiprocessing.cpu_count())  # 12 subjects = 12 workers max
LOG_DIR = DERIV_BASE / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Backward-compatible aliases so existing notebook cells remain unchanged.
run = run_command


def run_parallel(tasks, label="Processing"):
    return _run_parallel(
        tasks,
        label=label,
        max_workers=MAX_WORKERS,
        log_dir=LOG_DIR,
    )

print("Helpers loaded ✅")

# %% [Step 0] Verify tools: MRtrix3, ANTs, FreeSurfer, MATLAB
result_fs = run_command(f"{MODULE_LOAD_FS} && mri_info --version", check=False, verbose=False)
if result_fs.returncode == 0:
    print(f"✅ FreeSurfer found (module: {MODULE_LOAD_FS})")
else:
    print(f"❌ FreeSurfer NOT found")
    print(f"   Module load: {MODULE_LOAD_FS}")
    print(f"   stderr: {result_fs.stderr[:200]}")

result_mr = run_command(f"{MODULE_LOAD_MRTRIX} && mrconvert --version", check=False, verbose=False)
if result_mr.returncode == 0:
    print(f"✅ MRtrix3 found (module: {MODULE_LOAD_MRTRIX})")
else:
    print(f"❌ MRtrix3 NOT found")
    print(f"   Module load: {MODULE_LOAD_MRTRIX}")
    print(f"   stderr: {result_mr.stderr[:200]}")

result_ants = run_command(f"{MODULE_LOAD_ANTS} && N4BiasFieldCorrection --version", check=False, verbose=False)
if result_ants.returncode == 0 or 'ANTs' in (result_ants.stdout + result_ants.stderr):
    print(f"✅ ANTs found (module: {MODULE_LOAD_ANTS})")
else:
    print(f"❌ ANTs NOT found")
    print(f"   Module load: {MODULE_LOAD_ANTS}")
    print(f"   stderr: {result_ants.stderr[:200]}")

result_mat = run_command(f"{MODULE_LOAD_MATLAB} && matlab -batch \"disp('MATLAB OK')\"", check=False, verbose=False)
if result_mat.returncode == 0 and 'MATLAB OK' in result_mat.stdout:
    print(f"✅ MATLAB found (module: {MODULE_LOAD_MATLAB})")
else:
    print(f"❌ MATLAB NOT found")
    print(f"   Module load: {MODULE_LOAD_MATLAB}")
    print(f"   stderr: {result_mat.stderr[:200]}")


# %% [markdown] Step 0.9: Delete existing mask_new files (re-run safety)
# ---
# ## Step 0.9: Delete existing mask_new files
#
# If re-running Step 1 with new parameters, first delete existing masks
# so they aren't skipped. Only deletes mask_new.* — NOT mask_new_edited.mif.

# %% [Step 0.9] Delete existing mask_new (but NOT mask_new_edited) to allow re-run
deleted = 0
skipped_locked = 0
PROTECT_EDITED_MASKS = True
for s in subjects:
    # If a manual edited mask exists, keep all mask variants untouched for this subject.
    if PROTECT_EDITED_MASKS and s['mask_new_edited'].exists():
        skipped_locked += 1
        print(f"  🔒 Locked (has mask_new_edited): {s['id']}")
        continue

    for f in [s['mask_new'], s['mask_new_mif']]:
        if f.exists():
            f.unlink()
            deleted += 1
            print(f"  🗑️ Deleted: {f.name} ({s['id']})")
print(f"\nDeleted {deleted} mask files.")
print(f"Locked subjects skipped: {skipped_locked}")
print("mask_new_edited.mif files preserved.")

# %% [markdown] Step 1: Generate Brain Masks (Group-specific atlas registration)
# ---
# ## Step 1: Generate Brain Masks — Group-Specific Atlas Segmentation
#
# Each group has its own **manually corrected atlas** subject whose mask is
# registered to the remaining subjects in that group via ANTs SyN:
#
# | Group | Atlas subject           | # targets |
# |-------|-------------------------|-----------|
# | CTRL  | sub-129-921-R1F         | 3         |
# | AP    | sub-596-C-OO-NmM       | 3         |
# | APTR  | sub-AP72-862-901-NmF   | 3         |
#
# Pipeline per target:
#   1. **N4BiasFieldCorrection** (ANTs) — corrects intensity inhomogeneity
#   2. **antsRegistrationSyN.sh** (ANTs) — rigid + affine + SyN nonlinear
#   3. **antsApplyTransforms** (ANTs) — warp atlas mask (NearestNeighbor)
#   4. **maskfilter** (MRtrix) — keep largest connected component
#
# Atlas subjects keep their own `mask_new_edited.mif` as-is.
#
# **Output:** `mask_new.mif` and `mask_new.nii.gz` in each subject's dwi/ folder.

# %% [Step 1a] Prepare group-specific atlases: N4-correct each atlas T2w
import time

# ============================================================
# GROUP-SPECIFIC ATLAS DEFINITIONS — one manually corrected subject per group
# ============================================================
GROUP_ATLAS_IDS = dict(DEFAULT_GROUP_ATLAS_IDS)

# Build atlas info dict: {group: {t2w_n4, mask_nii, atlas_id, atlas_sub}}
atlas_info = build_group_atlas_info(subjects, GROUP_ATLAS_IDS)

print(f"{'=' * 70}")
print(f" GROUP-SPECIFIC ATLASES")
print(f"{'=' * 70}")
for group, info in atlas_info.items():
    print(f"  {group:>4}: {info['atlas_id']}")
    assert info['t2w'].exists(),     f"Atlas T2w not found: {info['t2w']}"
    assert info['mask_mif'].exists(), f"Atlas mask not found: {info['mask_mif']}"
    print(f"         T2w: ✅  Mask: ✅")

# N4-correct each atlas T2w and convert mask to NIfTI
for group, info in atlas_info.items():
    print(f"\nPreparing {group} atlas ({info['atlas_id']})...")
    run(f"{MODULES} && mrconvert {info['t2w']} {info['t2w_nii']} -force -quiet")
    run(f"{MODULES} && N4BiasFieldCorrection -d 3 -i {info['t2w_nii']} -o {info['t2w_n4']} "
        f"-s 2 -b [50] -c [200x200x200,0.0]")
    run(f"{MODULES} && mrconvert {info['mask_mif']} {info['mask_nii']} -force -quiet")
    check_file(info['t2w_n4'], f"{group} atlas T2w N4")
    check_file(info['mask_nii'], f"{group} atlas mask NIfTI")

print("\n✅ All 3 group atlases prepared.")

# %% [Step 1b] Register group atlas → each target and warp mask (PARALLEL)
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Dispatch parallel jobs ---
# ANTs SyN uses ~4 threads per job (-n 4). 9 targets total (3 per group).
# lmc02: 40 cores → 40/4 = 10 concurrent jobs (we have 9 targets + 3 atlas copies)
N_WORKERS = 9
ATLAS_IDS = set(GROUP_ATLAS_IDS.values())

print(f"{'=' * 70}")
print(f" GENERATING BRAIN MASKS (group-specific atlas registration)")
print(f" FOR {len(subjects)} SUBJECTS ({len(ATLAS_IDS)} atlases + {len(subjects) - len(ATLAS_IDS)} targets)")
print(f"{'=' * 70}")
for group, atlas_id in GROUP_ATLAS_IDS.items():
    targets = [s['id'] for s in subjects if s['group'] == group and s['id'] != atlas_id]
    print(f"  {group:>4} atlas: {atlas_id} → targets: {targets}")
print(f"\n Using {N_WORKERS} parallel workers on lmc02 (40 cores, 188 GB RAM)")
print(f" Registration: ANTs SyN (rigid + affine + nonlinear), ~3-5 min/subject")
print(f"{'=' * 70}")

t0 = time.time()
results = {'success': [], 'failed': [], 'skipped': []}
mask_tasks = build_mask_generation_args(subjects, atlas_info, MODULES)

with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
    futures = {}
    for args in mask_tasks:
        sid = args[0]

        if sid in ATLAS_IDS:
            # Atlas subject: just copy its edited mask
            futures[executor.submit(copy_atlas_mask, args)] = sid
        else:
            # Non-atlas: register group atlas → target and warp mask
            futures[executor.submit(register_and_warp_mask, args)] = sid

    for future in as_completed(futures):
        sid = futures[future]
        try:
            sid_ret, status, msg = future.result(timeout=1800)
            results[status].append(sid_ret)
            icon = {"success": "✅", "failed": "❌", "skipped": "⏭️"}[status]
            print(f"  {icon} {sid_ret}: {msg}")
        except Exception as e:
            results['failed'].append(sid)
            print(f"  ❌ {sid}: exception — {e}")

elapsed = time.time() - t0
print(f"\n{'=' * 70}")
print(f" MASK GENERATION RESULTS (group-specific atlas) — {elapsed:.0f}s total")
print(f"   ✅ Success: {len(results['success'])}")
print(f"   ❌ Failed:  {len(results['failed'])}")
print(f"   ⏭️  Skipped: {len(results['skipped'])}")
if results['failed']:
    print(f"   Failed subjects: {results['failed']}")
print(f"{'=' * 70}")

# %% [markdown] Step 2: Visual QC — Compare Old vs New Masks
# ---
# ## Step 2: Visual QC — Inspect Every Mask
#
# This step opens mrview for **each subject** with the atlas-based mask overlaid
# on the T2w image. Scroll through **every slice** and note any problems.
#
# **What to look for:**
# - Mask covers the **entire brain** without gaps
# - Mask **excludes skull/background**
# - For AP/APTR: check the **graft-host border** — the automatic mask may
#   include/exclude the graft incorrectly
#
# After reviewing, note which subjects need manual correction in the next step.

# %% [Step 2a] QC: View ONE subject to test (pick first CTRL)
# Quick test with the first subject to verify mrview works
#
# Copy-paste the command printed below into any terminal that has mrview.
#
s = subjects[0]
print(f"Opening QC for: {s['id']} ({s['group']})")
print(
    "Example command:\n  "
    f"mrview {s['t2w']} -overlay.load {s['mask_new_mif']} "
    "-overlay.opacity 0.3 -overlay.colour 0,1,0"
)
if s['mask_new_mif'].exists():
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 "
        f"-overlay.colour 0,1,0"
    )
    print("Green overlay = atlas-based mask. Scroll through all slices.")
else:
    print(f"❌ Mask not found: {s['mask_new_mif']}")

# %% [Step 2b] QC: Compare OLD vs NEW mask side-by-side (first subject)
s = subjects[0]
print(f"Comparing OLD (red) vs NEW (green) mask for: {s['id']}")
if s['mask_new_mif'].exists() and s['mask_old'].exists():
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_old']} -overlay.opacity 0.3 -overlay.colour 1,0,0 "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 -overlay.colour 0,1,0"
    )
    print("RED = old dwi2mask | GREEN = new atlas-based")
else:
    print("❌ One or both masks missing")

# %% [Step 2c] QC: Systematic review — ALL subjects, one by one
# Run this cell, review one subject, then press Enter to continue to the next.
# Record any subjects that need manual editing below.

SUBJECTS_TO_FIX = []  # ← Will be populated as you review

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    
    if not s['mask_new_mif'].exists():
        print(f"  ⏭️ {label}: no mask found, skipping")
        continue
    
    print(f"\n{'=' * 60}")
    print(f"  QC: {label}")
    print(f"{'=' * 60}")
    
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {s['mask_new_mif']} -overlay.opacity 0.3 "
        f"-overlay.colour 0,1,0"
    )
    
    answer = input(f"  {s['id']}: PASS (Enter) or FAIL (type 'f'): ").strip().lower()
    if answer == 'f':
        SUBJECTS_TO_FIX.append(s['id'])
        print(f"  → Marked for manual correction: {s['id']}")
    else:
        print(f"  → PASS ✅")

print(f"\n{'=' * 60}")
print(f" QC COMPLETE")
print(f"   Subjects needing manual correction: {len(SUBJECTS_TO_FIX)}")
for sid in SUBJECTS_TO_FIX:
    print(f"     - {sid}")
print(f"{'=' * 60}")

# %% [markdown] Step 3: Manual Mask Correction
# ---
# ## Step 3: Manual Mask Correction (for failed QC subjects)
#
# For subjects where the automatic mask got the graft-host border wrong:
# 1. Open the mask in mrview's ROI editor
# 2. Paint/erase voxels to fix the boundary
# 3. Save as `mask_new_edited.mif`
#
# **mrview ROI editor usage:**
# - Tool → ROI editor
# - Load the atlas-based mask as the ROI
# - Left-click = add voxels, Right-click = remove voxels
# - Adjust brush size with scroll wheel
# - File → Save to save the edited mask

# %% [Step 3a] Manual Correction — Set subjects to fix
# ============================================================
# EDIT THIS: List subjects that failed QC in Step 2c.
# If Step 2c was run interactively, SUBJECTS_TO_FIX is already set.
# Otherwise, manually specify them here:
# ============================================================

# Uncomment and edit if needed:
# SUBJECTS_TO_FIX = ["sub-AP72-862-901-NmF", "sub-Apallial17NmF"]

print(f"Subjects to manually correct: {SUBJECTS_TO_FIX}")
if not SUBJECTS_TO_FIX:
    print("None! All masks passed QC. Skip to Step 4.")

# %% [Step 3b] Open mrview ROI editor for each subject needing correction
for sid in SUBJECTS_TO_FIX:
    s = next((x for x in subjects if x['id'] == sid), None)
    if s is None:
        print(f"❌ Subject {sid} not found in subject list!")
        continue
    
    src = s['mask_new_mif']
    dst = s['mask_new_edited']
    
    if not src.exists():
        print(f"❌ {sid}: mask not found")
        continue
    
    # Copy mask to edited version (so we don't lose the original)
    if not dst.exists():
        run(f"cp {src} {dst}", verbose=False)
    
    print(f"\n{'=' * 60}")
    print(f"  EDITING: {sid} ({s['group']})")
    print(f"  Mask file: {dst}")
    print(f"{'=' * 60}")
    print("  Instructions:")
    print("    1. In mrview: Tool → ROI editor")
    print("    2. Load ROI: the mask_new_edited.mif")
    print("    3. Paint (left-click) / Erase (right-click) to fix boundaries")
    print("    4. Save (Ctrl+S or File → Save) when done")
    print("    5. Close mrview and press Enter to continue")
    
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {dst} -overlay.opacity 0.4 -overlay.colour 0,1,0"
    )
    
    input(f"  Press Enter when done editing {sid}...")
    
    if dst.exists():
        print(f"  ✅ Saved: {dst}")
    else:
        print(f"  ⚠️ Edited mask not found at {dst}")

print("\n✅ Manual correction complete.")

# %% [markdown] Step 4: Apply New Masks — Replace Old Masks
# ---
# ## Step 4: Apply New Masks
#
# For each subject:
# - Back up old `mask.mif` and `mask_anat.mif`
# - **Resample** the new mask from T2w resolution (46×134×160, 0.3mm) to DWI
#   resolution (28×100×107, 0.5mm) using `mrgrid regrid -interp nearest`
# - Install the resampled mask as both `mask.mif` and `mask_anat.mif`
#
# **Why resample?** The new mask was drawn on T2w_highres_flip (0.3×0.1×0.1mm),
# but all downstream processing (FOD, tractography, fixels) operates at DWI
# resolution (0.5×0.15×0.15mm). MRtrix commands like dwi2fod, tckgen, and
# mtnormalise expect the mask to match the DWI grid exactly.
#
# **Why NOT re-run DWI preprocessing (step 4)?** Running batch_process_flip1.py
# --steps 4 would overwrite our manually edited mask with a fresh dwi2mask
# output. The DWI denoising, bias correction, and registration are already
# fine from the original pipeline — only the mask needs replacing.

# %% [Step 4] Back up old masks and install new ones (resampled to DWI space)
BACKUP_SUFFIX = datetime.now().strftime("%Y%m%d_%H%M%S")
install_plans = build_mask_install_plan(subjects, BACKUP_SUFFIX)

for plan in install_plans:
    label = plan["label"]
    new_mask = plan["new_mask"]
    source = plan["source"]
    dwi_template = plan["dwi_template"]
    mask_mif = plan["mask_mif"]
    mask_anat = plan["mask_anat"]
    mask_resampled = plan["mask_resampled"]

    if new_mask is None:
        print(f"  ⏭️ {label}: no new mask available, skipping")
        continue
    if dwi_template is None:
        print(f"  ❌ {label}: no DWI template found for resampling, skipping")
        continue

    # Backup old masks
    if mask_mif.exists():
        run(f"mv {mask_mif} {plan['backup_mask']}", verbose=False)
    if mask_anat.exists():
        run(f"mv {mask_anat} {plan['backup_mask_anat']}", verbose=False)

    # Resample new mask from T2w resolution to DWI resolution
    result = run(
        f"{MODULES} && mrgrid {new_mask} regrid {mask_resampled} "
        f"-template {dwi_template} -interp nearest -force -quiet",
        verbose=False,
    )

    if result.returncode != 0:
        print(f"  ❌ {label}: mrgrid resample FAILED")
        continue

    # Install resampled mask as both mask.mif and mask_anat.mif
    run(f"cp {mask_resampled} {mask_mif}", verbose=False)
    run(f"cp {mask_resampled} {mask_anat}", verbose=False)

    print(f"  ✅ {label}: installed {source} mask (resampled to DWI grid)")

print(f"\n✅ All masks updated. Old masks backed up with suffix: {BACKUP_SUFFIX}")

# %% [Step 4-verify] Verify installed masks match DWI grid
print(f"{'=' * 70}")
print(f" VERIFYING INSTALLED MASKS MATCH DWI GRID")
print(f"{'=' * 70}")

grid_rows = build_mask_grid_report(subjects)
for row in grid_rows:
    sid = row["subject_id"]
    if row.get("missing_mask"):
        print(f"  ❌ {sid}: mask.mif not found!")
        continue
    mask_dims = row.get("mask_dims", "")
    if row.get("has_fod"):
        fod_dims = row.get("fod_dims", "")
        match = "✅" if row.get("match") else "❌ MISMATCH"
        print(f"  {match} {sid}: mask={mask_dims}  FOD={fod_dims}")
    else:
        print(f"  ⚠️  {sid}: mask={mask_dims}  (no FOD to compare)")

# %% [Step 4-checkpoint] Check pipeline readiness — what needs re-running?
# Run this BEFORE any Step 5 cell to see what still needs processing.
print(f"{'=' * 70}")
print(f" PIPELINE CHECKPOINT — WHAT NEEDS RE-RUNNING?")
print(f"{'=' * 70}")

chk = compute_pipeline_checkpoint(subjects)

print(f"\n  MASKS (DWI resolution):  ✅ {len(chk['mask_ok'])}/12")
if chk['mask_bad']:
    print(f"    ❌ BAD: {chk['mask_bad']}  → Run Step 4 first!")

needs_fod = chk['fod_stale'] + chk['fod_missing']
print(f"\n  FODs (wm_fod_norm.mif):  ✅ {len(chk['fod_ok'])}/12  stale/missing={len(needs_fod)}")
if needs_fod:
    print(f"    → Subjects needing FOD: {needs_fod}")

needs_grp = chk['grp_stale'] + chk['grp_missing']
print(f"\n  GROUP FODs:              ✅ {len(chk['grp_ok'])}/12  stale/missing={len(needs_grp)}")

needs_tck = chk['tck_stale'] + chk['tck_missing']
print(f"\n  TRACTOGRAMS:             ✅ {len(chk['tck_ok'])}/12  stale/missing={len(needs_tck)}")

print(f"\n{'=' * 70}")
print(f" RECOMMENDED NEXT STEP:")
if chk['mask_bad']:
    print(f"   → Step 4: Install masks (some missing/wrong resolution)")
elif needs_fod:
    print(f"   → Step 5b: FOD estimation ({len(needs_fod)} subjects)")
elif needs_grp:
    print(f"   → Step 5b-group: Group response & normalization")
elif needs_tck:
    print(f"   → Step 5d: Tractography ({len(needs_tck)} subjects)")
else:
    print(f"   → Everything up to date! Check Steps 5e+ for downstream.")
print(f"{'=' * 70}")

# %% [markdown] Step 5: Re-run Downstream Pipeline
# ---
# ## Step 5: Re-run Downstream Pipeline
#
# With new masks installed (resampled to DWI resolution), we re-run everything
# that depends on the mask. We skip DWI preprocessing (step 4 in
# batch_process_flip1.py) because:
# - It would overwrite our mask with a fresh dwi2mask output
# - The DWI denoising, bias correction, and registration are unchanged
# - The N4 bias correction quality difference from using a slightly different
#   mask boundary is negligible
#
# **Parallelization:** All per-subject steps use `run_parallel()` with up to
# 12 workers (one per subject) on lmc02's 40 cores.
#
# | Step | What | Why re-run |
# |------|------|------------|
# | 5b | FOD estimation | Uses mask for response estimation + CSD |
# | 5b-group | Group response + normalization | Uses mask for mtnormalise |
# | 5c | ROI resampling | FOD space reference |
# | 5d | Whole-brain tractography | Uses mask as seed/constraint |
# | 5e | SIFT2 | Depends on tractogram |
# | 5f | Passthrough extraction | Depends on tractogram |
# | 5g | MATLAB cleaning | Depends on tracts |
# | 5h | Metrics report | Depends on cleaned tracts |
# | 5i | Fixel analysis | Uses group-normalized FOD + mask |

# %% [Step 5b] Re-run FOD estimation (step 5) — PARALLEL (12 subjects)
print(f"{'=' * 70}")
print(f" RE-RUNNING FOD ESTIMATION WITH NEW MASKS (PARALLEL)")
print(f" Workers: {MAX_WORKERS} on {multiprocessing.cpu_count()} cores")
print(f"{'=' * 70}")

fod_tasks = build_batch_process_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    steps="5",
    timeout=1800,
)
fod_results = run_parallel(fod_tasks, label="FOD estimation")

# %% [markdown] Step 5b-group: Group Response & Normalized FOD
# ---
# ## Step 5b-group: Recompute Group Response Functions & Group-Normalized FODs
#
# **Why this step is critical:**
# The fixel analysis (Step 5i) and any cross-subject AFD comparison use
# `wm_fod_norm_group.mif` — FODs estimated with a **group-averaged** response
# function, then normalized with `mtnormalise`. This ensures that AFD/FD values
# are comparable between subjects and between groups (CTRL vs AP vs APTR).
#
# Changing masks affects ALL of:
# 1. `dwi2response dhollander` — the mask constrains which voxels are used to
#    estimate per-subject response functions (WM, GM, CSF)
# 2. `dwi2fod msmt_csd` — the mask constrains where FODs are estimated
# 3. `mtnormalise` — the mask defines the brain volume for intensity normalization
#
# **Pipeline:**
# 1. Collect updated per-subject response functions (from Step 5b)
# 2. Average them across ALL subjects → `group_responses/group_response_{wm,gm,csf}.txt`
# 3. Re-estimate FODs per subject using the group response → `wm_fod_group.mif`
# 4. Re-normalize per subject → `wm_fod_norm_group.mif`
#
# If `mtnormalise` fails (known issue for some mouse data), the script falls
# back to using the unnormalized group FOD.

# %% [Step 5b-group] Recompute Group Response & Normalized FOD
print(f"{'=' * 70}")
print(f" RECOMPUTING GROUP RESPONSE FUNCTIONS & GROUP-NORMALIZED FODS")
print(f" (This uses ALL {len(subjects)} subjects for response averaging)")
print(f"{'=' * 70}")

result_grp = run(
    f"{MODULES} && cd {REPO_DIR} && python functions/compute_group_response.py "
    f"--force --derivatives-dir {DERIV_BASE}",
    timeout=7200
)

if result_grp.returncode == 0:
    print(f"\n  ✅ Group response & FOD normalization complete")
else:
    print(f"\n  ❌ Group response computation FAILED — check output above")

# Verify wm_fod_norm_group.mif exists for each subject
grp_fod_ok = 0
grp_fod_missing = 0
for s in subjects:
    fod_grp = s['dwi_dir'] / 'wm_fod_norm_group.mif'
    if fod_grp.exists():
        grp_fod_ok += 1
    else:
        grp_fod_missing += 1
        print(f"  ❌ MISSING: {s['id']} wm_fod_norm_group.mif")

print(f"\n  wm_fod_norm_group.mif: ✅ {grp_fod_ok} present, ❌ {grp_fod_missing} missing")

# %% [Step 5b-verify] Verify FOD estimation — dates, sizes, grids, normalisation
import datetime, json as _json
print(f"{'=' * 70}")
print(f" STEP 5b VERIFICATION — FOD Estimation & Group Normalisation")
print(f"{'=' * 70}\n")

_fod_files = [
    'wm_response.txt', 'gm_response.txt', 'csf_response.txt',  # per-subject responses
    'wm_fod.mif', 'gm_fod.mif', 'csf_fod.mif',                # per-subject FODs
    'wm_fod_norm.mif',                                          # per-subject normalised FOD
    'wm_fod_group.mif', 'gm_fod_group.mif', 'csf_fod_group.mif',  # group FODs
    'wm_fod_norm_group.mif',                                    # group-normalised FOD (key output)
]
_mask_files = ['mask.mif', 'mask_anat.mif', 'mask_new_dwi_space.mif']

_table_rows, _issues = collect_fod_verification(subjects, run_cmd=run, modules=MODULES)

# -- Group response files --
print("  Group response files:")
for entry in collect_group_response_status(DERIV_BASE):
    if entry['exists']:
        print(f"    ✅ {entry['name']}  (modified: {entry['mtime']})")
    else:
        print(f"    ❌ {entry['name']}  MISSING")
        _issues.append(f"❌ Group response {entry['name']} MISSING")

# -- Print per-subject table --
print(f"\n  {'Subject':<28} {'Group':<6} {'Status':<8} {'Norm?':<30} {'FOD MB':<8} "
      f"{'FOD Grid':<20} {'Grid?':<4} {'Resp':<4} {'FOD Date':<16} {'Mask Date':<16}")
print(f"  {'─'*28} {'─'*6} {'─'*8} {'─'*30} {'─'*8} {'─'*20} {'─'*4} {'─'*4} {'─'*16} {'─'*16}")
for r in _table_rows:
    print(f"  {r.get('subject','?'):<28} {r.get('group','?'):<6} {r.get('status','?'):<8} "
          f"{r.get('normalised','?'):<30} {r.get('size_MB','?'):<8} "
          f"{r.get('fod_grid','?'):<20} {r.get('grid_match','?'):<4} "
          f"{r.get('responses','?'):<4} {r.get('mtime','?'):<16} {r.get('mask_date','?'):<16}")

# -- Summary --
print(f"\n  {'=' * 60}")
if _issues:
    print(f"  ⚠️  {len(_issues)} issue(s) found:")
    for iss in _issues:
        print(f"    {iss}")
else:
    print(f"  ✅ ALL {len(subjects)} subjects passed Step 5b verification")
print(f"  {'=' * 60}")

# %% [Step 5c] Re-run ROI resampling + segmentation (steps 6,7) — PARALLEL
print(f"{'=' * 70}")
print(f" RE-RUNNING ROI RESAMPLING + SEGMENTATION (PARALLEL)")
print(f"{'=' * 70}")

seg_tasks = build_batch_process_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    steps="6,7",
    timeout=1800,
)
seg_results = run_parallel(seg_tasks, label="ROI resampling")

# %% [markdown] Step 5c-QC: ROI Space & Laterality Verification
# ---
# ## Step 5c-QC: ROI Space & Laterality Verification
#
# After ROI resampling we must verify:
# 1. **Correct space** — fod_space ROI grid matches the FOD (dims & voxel size)
# 2. **Non-empty** — every ROI has >0 voxels after resampling
# 3. **Left/Right separation** — left and right ROIs have centroids on opposite
#    sides of axis 2 (L/R axis after flip1), with no overlap
# 4. **Within-subject consistency** — all "left" ROIs are on the same side,
#    all "right" ROIs on the opposite side
#
# **Axis convention (flip1):**
# - axis 0 = AP, axis 1 = DV, axis 2 = LR
# - CTRL (atlas):  left → LOW axis-2 indices,  right → HIGH
# - AP/APTR (manual+flip): left → HIGH axis-2 indices, right → LOW

# %% [Step 5c-QC] ROI Space & Laterality Check
import numpy as np

print(f"{'=' * 80}")
print(f" ROI SPACE & LATERALITY VERIFICATION — {len(subjects)} SUBJECTS")
print(f"{'=' * 80}")

roi_qc_issues = []  # collect (subject, issue_str)

for si, s in enumerate(subjects):
    label = f"[{si+1}/{len(subjects)}] {s['id']} ({s['group']})"
    print(f"\n{'─' * 70}")
    print(f"  {label}")
    print(f"{'─' * 70}")

    DWI = s['dwi_dir']
    fod_file = DWI / 'wm_fod_norm.mif'
    mask_file = DWI / 'mask.mif'

    # --- 1. Reference grid from FOD ---
    if not fod_file.exists():
        msg = "FOD not found — skipping"
        print(f"  ❌ {msg}")
        roi_qc_issues.append((s['id'], msg))
        continue

    fod_dims, fod_vox = get_grid_info(fod_file)
    print(f"  FOD grid:  {fod_dims}  voxel: {[f'{v:.4f}' for v in fod_vox]}")

    # --- 2. Enumerate fod_space ROIs ---
    if s['group'] == 'CTRL':
        roi_base = s['base'] / 'atlas' / 'rois' / 'fod_space'
    else:
        roi_base = s['base'] / 'manual_roi' / 'fod_space'
    roi_pairs = get_roi_pairs_for_group(s['group'])

    left_centroids_ax2 = []
    right_centroids_ax2 = []

    for region, (lf, rf) in roi_pairs.items():
        for hemi_tag, fname in [('L', lf), ('R', rf)]:
            roi_path = roi_base / fname
            if not roi_path.exists():
                msg = f"{fname}: NOT FOUND"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
                continue

            # Check grid match
            roi_dims, roi_vox = get_grid_info(roi_path)
            grid_ok = (roi_dims == fod_dims)
            vox_ok = all(abs(rv - fv) < 0.01 for rv, fv in zip(roi_vox, fod_vox))

            if not grid_ok:
                msg = f"{fname}: GRID MISMATCH dims={roi_dims} vs FOD={fod_dims}"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
            if not vox_ok:
                msg = f"{fname}: VOXEL SIZE MISMATCH vox={roi_vox} vs FOD={fod_vox}"
                print(f"    ⚠️  {msg}")
                roi_qc_issues.append((s['id'], msg))

            # Load data and check non-empty + centroid
            data = load_roi_data(roi_path)
            coords = np.argwhere(data > 0.5)
            nvox = len(coords)

            if nvox == 0:
                msg = f"{fname}: EMPTY (0 voxels)"
                print(f"    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))
                continue

            centroid = coords.mean(axis=0)
            mid2 = data.shape[2] / 2.0
            side_str = "HIGH" if centroid[2] > mid2 else "LOW"

            status = "✅" if grid_ok else "⚠️ "
            print(f"    {status} {region:>12}_{hemi_tag}: {nvox:>5} vox  "
                  f"centroid_ax2={centroid[2]:5.1f}  mid={mid2:.1f}  side={side_str}")

            if hemi_tag == 'L':
                left_centroids_ax2.append((region, centroid[2]))
            else:
                right_centroids_ax2.append((region, centroid[2]))

    # --- 3. Left/Right consistency check ---
    if left_centroids_ax2 and right_centroids_ax2:
        laterality = evaluate_laterality_consistency(
            left_centroids_ax2,
            right_centroids_ax2,
            axis_mid=fod_dims[2] / 2.0,
        )

        if laterality['lr_separated']:
            print(f"\n    ✅ L/R laterality: all L={laterality['left_sides'][0]}, all R={laterality['right_sides'][0]} — CONSISTENT")
        else:
            for msg in laterality['issues']:
                print(f"\n    ❌ {msg}")
                roi_qc_issues.append((s['id'], msg))

        # Check L/R overlap for each paired region
        for region, (lf, rf) in roi_pairs.items():
            lp = roi_base / lf
            rp = roi_base / rf
            if lp.exists() and rp.exists():
                ld = load_roi_data(lp)
                rd = load_roi_data(rp)
                overlap = np.sum((ld > 0.5) & (rd > 0.5))
                if overlap > 0:
                    msg = f"{region}: L/R OVERLAP = {overlap} voxels"
                    print(f"    ⚠️  {msg}")
                    roi_qc_issues.append((s['id'], msg))

# --- Summary ---
print(f"\n{'=' * 80}")
print(f" ROI QC SUMMARY")
print(f"{'=' * 80}")
if roi_qc_issues:
    print(f"  ⚠️  {len(roi_qc_issues)} issues found:")
    for subj, issue in roi_qc_issues:
        print(f"    {subj:<30} {issue}")
else:
    print(f"  ✅ All ROIs passed: correct space, non-empty, L/R consistent")
print(f"{'=' * 80}")

# %% [Step 5d] Re-run whole-brain tractography — PARALLEL
# NOTE: tckgen is CPU-intensive. With 12 subjects × ~3 threads each = 36 cores used.
# lmc02 has 40 cores so this saturates the machine — that's fine, it's the fastest way.
# Uses generate_wholebrain_tractogram.py (1M streamlines, iFOD2, mask-seeded)
print(f"{'=' * 70}")
print(f" RE-RUNNING WHOLE-BRAIN TRACTOGRAPHY (PARALLEL)")
print(f" 12 subjects × tckgen 1M streamlines on {multiprocessing.cpu_count()} cores")
print(f"{'=' * 70}")

tck_tasks = build_wholebrain_tractography_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=None,
)
tck_results = run_parallel(tck_tasks, label="Tractography")

# %% [Step 5d-post] Downsample to 10k + convert to .trk
# Downsample wholebrain_1000k.tck → wholebrain_10k.tck using tckedit -number
# Then convert wholebrain_10k.tck → wholebrain_10k.trk using nibabel
print(f"{'=' * 70}")
print(f" DOWNSAMPLING TO 10k + CONVERTING TO .trk — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# --- Part 1: Downsample with tckedit (parallel) ---
ds_tasks = build_tck_downsample_tasks(
    subjects,
    modules=MODULES,
    derivatives_dir=DERIV_BASE,
    number=10000,
    timeout=None,
)
ds_results = run_parallel(ds_tasks, label="Downsample 10k")

# --- Part 2: Convert .tck → .trk using nibabel ---
print(f"\n  Converting .tck → .trk ...")
trk_ok, trk_fail = 0, 0
for s in subjects:
    sid = s['id']
    ses = s['session']
    dwi_dir = DERIV_BASE / sid / ses / "dwi"
    wb_dir = dwi_dir / "wholebrain"
    tck_10k = wb_dir / "wholebrain_10k.tck"
    trk_10k = wb_dir / "wholebrain_10k.trk"
    t2w_mif = dwi_dir / "t2w_highres_flip.mif"

    if not tck_10k.exists():
        print(f"  ❌ {sid}: wholebrain_10k.tck not found")
        trk_fail += 1
        continue

    try:
        n_streamlines = convert_tck_to_trk(
            tck_path=tck_10k,
            trk_path=trk_10k,
            reference_mif_path=t2w_mif,
            modules_cmd=MODULES,
        )
        trk_ok += 1
        print(f"  ✅ {sid}: {trk_10k.name} ({n_streamlines} streamlines)")
    except Exception as e:
        trk_fail += 1
        print(f"  ❌ {sid}: {e}")

print(f"\n  .trk conversion: ✅ {trk_ok} / ❌ {trk_fail}")

# --- QC: visualize in mrview ---
# Run these on BRODMANN126 to visually check the 10k tractograms:
#
# D={DERIV_BASE}
# for d in $D/sub-*/ses-*/dwi/wholebrain; do
#   sub=$(basename $(dirname $(dirname $(dirname "$d"))))
#   mrview $(dirname "$d")/t2w_highres_flip.mif \
#     -tractography.load "$d/wholebrain_10k.tck" &
#   sleep 0.5
# done

# %% [markdown] Step 5e: SIFT2 Weights
# ---
# ## Step 5e: SIFT2 Weights for Whole-Brain Tractogram
#
# Compute SIFT2 cross-section multipliers so each streamline gets a weight
# proportional to how well the tractogram matches the underlying FOD.

# %% [Step 5e] Compute SIFT2 Weights — PARALLEL
print(f"{'=' * 70}")
print(f" COMPUTING SIFT2 WEIGHTS (PARALLEL) — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# Check timestamps: if sift2_weights.txt is older than the tractogram, remove it
# so the pipeline script will recompute it instead of skipping.
import os
stale_count = 0
for s in subjects:
    wb_tck = s['wb_tck']
    wf = s['sift2_weights']
    if wf.exists() and wb_tck.exists():
        tck_mtime = os.path.getmtime(wb_tck)
        wf_mtime  = os.path.getmtime(wf)
        if wf_mtime < tck_mtime:
            print(f"  STALE: {s['id']} — weights ({datetime.datetime.fromtimestamp(wf_mtime):%Y-%m-%d %H:%M}) "
                  f"older than tractogram ({datetime.datetime.fromtimestamp(tck_mtime):%Y-%m-%d %H:%M}), removing")
            wf.unlink()
            mu = s.get('sift2_mu') or (wf.parent / "sift2_mu.txt")
            if mu.exists():
                mu.unlink()
            stale_count += 1
        else:
            print(f"  OK:    {s['id']} — weights ({datetime.datetime.fromtimestamp(wf_mtime):%Y-%m-%d %H:%M}) "
                  f"up-to-date vs tractogram ({datetime.datetime.fromtimestamp(tck_mtime):%Y-%m-%d %H:%M})")
    elif not wf.exists() and wb_tck.exists():
        print(f"  NEW:   {s['id']} — no weights file, will compute")
        stale_count += 1
    elif not wb_tck.exists():
        print(f"  SKIP:  {s['id']} — no tractogram found")
print(f"\n  → {stale_count} subject(s) need SIFT2 (re)computation\n")

sift2_tasks = build_sift2_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=1800,
)
sift2_results = run_parallel(sift2_tasks, label="SIFT2")

# %% [markdown] Step 5e-CST: Create & Edit CST ROIs for AP/APTR Subjects
# ---
# ## Step 5e-CST: CST (Corticospinal Tract) ROI Preparation
#
# Before passthrough extraction, AP and APTR subjects need **CST ROIs** (left & right)
# for the corticospinal tract analysis. These are manually drawn on the t2w_highres_flip
# image using mrview's ROI editor, saved to `manual_roi_flip/`, then regridded to FOD space.
#
# Source dir:  `PROYECTOS/mbrain/local/manual_roi_flip/`
# Target dir:  `<subject>/dwi/wholebrain/passthrough_cst/`
#
# Workflow:
# 1. Check which AP/APTR subjects already have CST ROIs
# 2. Create empty template ROIs for subjects that are missing them
# 3. Open mrview for manual drawing (on BROADMANN126)
# 4. Regrid finished ROIs to FOD space

# %% [Step 5e-CST-check] Check CST ROI status for AP/APTR subjects
import os
from pathlib import Path

# Fix MRtrix tmp dir (system default points to non-existent /flywheel/v0/output/tmp)
os.environ['MRTRIX_TMPDIR'] = '/tmp'

# CST ROIs live in each subject's manual_roi/ dir as roi_CST_{hemi}_flip.mif
# (same location & naming scheme as roi_Graft_left_flip.mif, roi_Thalamus_left_flip.mif, etc.)
# Legacy source dir is config-driven (paths.cst_legacy_manual_roi_flip_dir).
CST_LEGACY_DIR = Path(paths["cst_legacy_manual_roi_flip_dir"])

print(f"{'=' * 70}")
print(f" CST ROI STATUS CHECK — AP & APTR SUBJECTS")
print(f"{'=' * 70}")
print(f"  ROI naming: manual_roi/roi_CST_{{left,right}}_flip.mif")
print(f"  FOD-space:  wholebrain/passthrough_cst/roi_CST_{{left,right}}_fod.mif")
print()

cst_missing = []   # subjects that need CST ROIs created
cst_exists  = []   # subjects that already have CST ROIs in FOD space
cst_needs_regrid = []  # subjects that have source ROIs but no FOD-space ROIs

status_snapshot = collect_cst_status(subjects, legacy_dir=CST_LEGACY_DIR)
cst_exists = status_snapshot['complete']
cst_needs_regrid = status_snapshot['needs_regrid']
cst_missing = status_snapshot['missing']

for s, status, _src_left, _src_right, _fod_left, _fod_right in status_snapshot['rows']:
    print(f"  [{s['group']:>4}] {s['id']:<30} {status}")

print(f"\n  Summary:")
print(f"    ✅ Complete (FOD-space):  {len(cst_exists)}")
print(f"    ⚠️  Need regrid only:     {len(cst_needs_regrid)}")
print(f"    ❌ Need creation & draw:  {len(cst_missing)}")

# %% [Step 5e-CST-create] Create empty CST ROI templates for subjects that need them
# Creates empty bitwise ROIs in t2w_highres_flip space in each subject's manual_roi/
# for manual drawing in mrview. Naming: roi_CST_{hemi}_flip.mif
os.environ['MRTRIX_TMPDIR'] = '/tmp'

print(f"{'=' * 70}")
print(f" CREATING EMPTY CST ROI TEMPLATES FOR MANUAL DRAWING")
print(f"{'=' * 70}")

created_files = []
for s in cst_missing:
    sid = s['id']
    t2w = s['t2w']
    manual_roi_dir = s['base'] / 'manual_roi'

    if not t2w.exists():
        print(f"  ❌ {sid}: t2w not found at {t2w}")
        continue

    manual_roi_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n  {sid} ({s['group']}):")

    for hemi in ['left', 'right']:
        out_file = manual_roi_dir / f"roi_CST_{hemi}_flip.mif"
        if out_file.exists():
            print(f"    SKIP: {out_file.name} already exists")
            continue
        # Create an empty (all-zero) bitwise ROI matching t2w grid
        cmd = f"mrcalc {t2w} 0 -mult {out_file} -datatype bit -force -quiet"
        print(f"    Creating: {out_file.name} ... ", end="", flush=True)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        if result.returncode == 0 and out_file.exists():
            created_files.append(out_file)
            print(f"✅ ({out_file.stat().st_size} bytes)")
        else:
            print(f"❌ {result.stderr[:200]}")

print(f"\n  Created {len(created_files)} empty ROI templates")
print(f"  Now open them in mrview (next cell) to draw the CST ROIs manually.")

# %% [Step 5e-CST-edit] Open mrview for manual CST ROI drawing — run on BROADMANN126
# Opens mrview with:
#   - t2w_highres_flip.mif as the base image
#   - Existing manual ROIs as overlays for reference (Graft, Thalamus)
#   - Empty CST ROI loaded in the ROI editor for drawing
#
# Instructions:
#   1. Tool → ROI editor
#   2. The CST ROI is pre-loaded — select it in the ROI editor panel
#   3. Draw the CST region (left or right) on coronal/axial slices
#   4. SAVE (Ctrl+S or File → Save) before closing
#   5. Repeat for the other hemisphere
print(f"{'=' * 70}")
print(f" MANUAL CST ROI DRAWING — mrview on BROADMANN126")
print(f"{'=' * 70}")
print()
print("INSTRUCTIONS:")
print("  1. In mrview: Tool → ROI editor")
print("  2. The CST ROI is pre-loaded in the editor")
print("  3. Draw the CST region on coronal/axial slices")
print("  4. SAVE before closing mrview")
print("  5. Repeat for each hemisphere & subject")
print()

# Collect all AP/APTR subjects needing manual drawing
all_need_draw = []
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    # Skip if FOD-space ROIs already exist
    if (cst_dir / "roi_CST_left_fod.mif").exists() and (cst_dir / "roi_CST_right_fod.mif").exists():
        continue

    src_left, src_right = resolve_cst_source_rois(s['id'], s['base'] / 'manual_roi', CST_LEGACY_DIR)
    for hemi, roi_file in [('left', src_left), ('right', src_right)]:
        if roi_file is not None and roi_file.exists():
            all_need_draw.append((s, hemi, roi_file))

for s, hemi, roi_file in all_need_draw:
    sid = s['id']
    t2w = s['t2w']

    # Build overlay list: show existing manual ROIs for reference
    overlays = ""
    manual_roi_dir = s['base'] / 'manual_roi'
    for roi_name in ['roi_Graft_left_flip.mif', 'roi_Graft_right_flip.mif',
                     'roi_Thalamus_left_flip.mif', 'roi_Thalamus_right_flip.mif']:
        roi_path = manual_roi_dir / roi_name
        if roi_path.exists():
            overlays += f" -overlay.load {roi_path} -overlay.opacity 0.3"

    print(f"\n{'=' * 60}")
    print(f"  EDITING: {sid} ({s['group']}) — {hemi} CST")
    print(f"  ROI file: {roi_file}")
    print(f"{'=' * 60}")
    print("  Instructions:")
    print("    1. In mrview: Tool → ROI editor")
    print("    2. The CST ROI is pre-loaded — select it in the editor panel")
    print("    3. Draw the CST region on coronal/axial slices")
    print("    4. Save (Ctrl+S or File → Save) when done")
    print("    5. Close mrview and press Enter to continue")

    launch_gui(f"mrview {t2w}{overlays} -roi.load {roi_file} -roi.opacity 0.5")

    input(f"  Press Enter when done editing {sid} {hemi} CST...")

    if roi_file.exists() and roi_file.stat().st_size > 1000:
        print(f"  ✅ Saved: {roi_file.name} ({roi_file.stat().st_size} bytes)")
    else:
        print(f"  ⚠️  ROI file may still be empty — check before regridding")

if not all_need_draw:
    print("  All AP/APTR subjects already have CST ROIs — nothing to draw!")
else:
    print(f"\n  ✅ Finished launching mrview for {len(all_need_draw)} CST ROIs.")
    print("  Run the next cell to regrid them to FOD space.")

# %% [Step 5e-CST-regrid] Regrid CST ROIs to FOD space and copy to passthrough_cst/
# After manual drawing, regrid source ROIs from t2w space to FOD space using
# mrgrid regrid with the subject's wm_fod_norm_group.mif as template.
os.environ['MRTRIX_TMPDIR'] = '/tmp'

print(f"{'=' * 70}")
print(f" REGRIDDING CST ROIs TO FOD SPACE")
print(f"{'=' * 70}")

regrid_plan = build_cst_regrid_plan(subjects, legacy_dir=CST_LEGACY_DIR)
all_to_regrid = regrid_plan['to_regrid']
for s in regrid_plan['already_complete']:
    print(f"  SKIP: {s['id']} — FOD-space CST ROIs already exist")
for s in regrid_plan['missing_sources']:
    print(f"  ❌ {s['id']}: Source ROIs not found (checked manual_roi/ and legacy dir)")

regrid_ok = 0
for s, src_left, src_right in all_to_regrid:
    sid = s['id']
    fod_template = s['dwi_dir'] / "wm_fod_norm_group.mif"
    cst_dir = s['wb_dir'] / "passthrough_cst"
    cst_dir.mkdir(parents=True, exist_ok=True)

    if not fod_template.exists():
        print(f"  ❌ {sid}: No FOD template at {fod_template}")
        continue

    print(f"  {sid}:")
    ok = True
    for hemi, src in [('left', src_left), ('right', src_right)]:
        dst = cst_dir / f"roi_CST_{hemi}_fod.mif"
        cmd = (f"mrgrid {src} regrid {dst} "
               f"-template {fod_template} -interp nearest -force -quiet")
        print(f"    {hemi}: {src.parent.name}/{src.name} → {dst.name} ... ", end="", flush=True)
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and dst.exists():
            print("✅")
        else:
            print(f"❌ {result.stderr[:200]}")
            ok = False
    if ok:
        regrid_ok += 1

print(f"\n  Regridded {regrid_ok}/{len(all_to_regrid)} subjects successfully")

# %% [Step 5e-CST-QC] QC: Verify all AP/APTR subjects have CST ROIs in FOD space
print(f"{'=' * 70}")
print(f" CST ROI FINAL STATUS — ALL AP/APTR SUBJECTS")
print(f"{'=' * 70}")
qc_summary = build_cst_final_qc_summary(subjects, legacy_dir=CST_LEGACY_DIR)
all_ok = qc_summary['all_ok']
for row in qc_summary['rows']:
    s = row['subject']
    print(f"  {row['icon']} [{s['group']:>4}] {s['id']:<30} "
          f"src={'YES' if row['has_src'] else 'NO':>3}  "
          f"fod_L={'YES' if row['has_fod_left'] else 'NO':>3}  "
          f"fod_R={'YES' if row['has_fod_right'] else 'NO':>3}")
if all_ok:
    print("\n  ✅ All AP/APTR subjects have CST ROIs — ready for Step 5f")
else:
    print("\n  ⚠️  Some subjects still missing CST ROIs — complete drawing before Step 5f")

# %% [Step 5e-CST-extract] Extract CST Passthrough Tracts (AP/APTR only) — PARALLEL
# ────────────────────────────────────────────────────────────────────────
# For each AP/APTR subject, extract streamlines that pass through BOTH:
#   1. The CST ROI  (passthrough_cst/roi_CST_{hemi}_fod.mif)
#   2. The WHOLE graft ROI (manual_roi/fod_space/roi_Graft_{hemi}_flip_fod.mif)
#
# Uses tckedit with dual -include from the whole-brain tractogram.
# Also extracts using shrunk graft ROI if available.
# Output → passthrough_cst/graft_{hemi}_to_CST_{hemi}.tck
# ────────────────────────────────────────────────────────────────────────
print(f"{'=' * 70}")
print(f" EXTRACTING CST PASSTHROUGH TRACTS (AP/APTR only)")
print(f"{'=' * 70}")

HEMIS = ['left', 'right']

# First, delete old segmented CST tracts (from previous dorsal-segment approach)
for s in subjects:
    if s['group'] == 'CTRL':
        continue
    cst_dir = s['wb_dir'] / "passthrough_cst"
    if not cst_dir.exists():
        continue
    old_tcks = [t for t in cst_dir.glob('seg9_graft_*.tck')]
    if old_tcks:
        for t in old_tcks:
            t.unlink()
        print(f"  🗑  {s['id']}: deleted {len(old_tcks)} old segmented CST tracts")

cst_plan = build_cst_passthrough_extraction_tasks(
    subjects,
    modules=MODULES,
    hemis=tuple(HEMIS),
    timeout=None,
)
for msg in cst_plan['skips']:
    print(f"  {msg}")
cst_tasks = cst_plan['tasks']

if cst_tasks:
    cst_results = run_parallel(cst_tasks, label="CST_passthrough")
else:
    print("  Nothing to extract.")

# Show results
print(f"\n  CST tract counts per subject:")
for row in collect_cst_raw_tract_counts(subjects):
    s = row['subject']
    print(f"  {s['id']:<30}: {row['n_normal']} normal, {row['n_shrink']} shrink")

# %% [Step 5e-CST-tckgen] Targeted CST tractography: seed CST → include Graft (AP/APTR) — PARALLEL
# ────────────────────────────────────────────────────────────────────────
# Runs tckgen seeded from the CST ROI, with the whole graft as an -include ROI.
# Uses the same iFOD2 parameters as the whole-brain tractogram (step=0.02,
# angle=60, cutoff=0.1, minlength=2, maxlength=15) but with -select 10000
# to get a reasonable number of targeted streamlines.
#
# Output → passthrough_cst/tckgen_CST_{hemi}_to_graft_{hemi}.tck
# ────────────────────────────────────────────────────────────────────────
print(f"{'=' * 70}")
print(f" TARGETED CST TRACTOGRAPHY: seed CST → include Graft (AP/APTR)")
print(f"{'=' * 70}")

HEMIS = ['left', 'right']
CST_TCKGEN_SELECT = 10000  # number of streamlines to select

tckgen_plan = build_cst_tckgen_tasks(
    subjects,
    modules=MODULES,
    hemis=tuple(HEMIS),
    select=CST_TCKGEN_SELECT,
    timeout=None,
)
for msg in tckgen_plan['skips']:
    print(f"  {msg}")
tckgen_tasks = tckgen_plan['tasks']

if tckgen_tasks:
    tckgen_results = run_parallel(tckgen_tasks, label="CST_tckgen")
else:
    print("  Nothing to generate.")

# Show results
print(f"\n  tckgen CST tract counts per subject:")
for row in collect_cst_tckgen_counts(subjects, hemis=tuple(HEMIS)):
    s = row['subject']
    if row['exists']:
        print(f"  {s['id']:<30} {row['hemi']:<6}: {row['count']} streamlines")
    else:
        print(f"  {s['id']:<30} {row['hemi']:<6}: NOT FOUND")

# %% [TEMP] CST tckgen — exclude Striatum/Cerebellum/Paleocortex/Thalamus (Apallial17NmF only)
# ────────────────────────────────────────────────────────────────────────
# Takes the tckgen CST→Graft tracts and removes any streamline passing
# through Striatum, Cerebellum, Paleocortex or Thalamus (both hemispheres).
# Output saved to /tmp/cst_viz/ for quick visualization — not stored officially.
# ────────────────────────────────────────────────────────────────────────
import subprocess

sid = "sub-Apallial17NmF"
ses = "ses-240816091633"
base = DERIV_BASE / sid / ses
roi_dir = base / "manual_roi" / "fod_space"
cst_dir = base / "dwi" / "wholebrain" / "passthrough_cst"

# Exclusion ROIs (all 8: 4 regions × 2 hemispheres)
roi_summary = collect_temp_exclusion_rois(roi_dir)
for roi in roi_summary['found']:
    print(f"  ✅ Exclude: {roi.name}")
for name in roi_summary['missing']:
    print(f"  ⚠️  Missing: {name}")

out_dir = Path("/tmp/cst_viz")
out_dir.mkdir(exist_ok=True)

temp_plan = build_temp_cst_exclusion_tasks(
    cst_dir=cst_dir,
    out_dir=out_dir,
    output_prefix="Apallial17_CST",
    modules=MODULES,
    exclude_rois=roi_summary['found'],
)
for src_name in temp_plan['missing_sources']:
    print(f"  ❌ Source not found: {src_name}")

for task in temp_plan['tasks']:
    print(f"\n  Processing {task['hemi']}...")
    r = subprocess.run(task['cmd'], shell=True, capture_output=True, text=True, timeout=120)
    if r.returncode == 0:
        # Count remaining streamlines
        info = subprocess.run(f"tckinfo {task['dst']} 2>/dev/null | awk '/^ *count:/{{print $2}}'",
                              shell=True, capture_output=True, text=True)
        cnt = info.stdout.strip() or '?'
        print(f"  ✅ {task['dst'].name}: {cnt} streamlines remaining")
    else:
        print(f"  ❌ Error: {r.stderr[:200]}")

print(f"\n  Output in: {out_dir}/")

# %% [Step 5e-CST-clean] Clean CST Passthrough Tracts with MATLAB AFQ — PARALLEL
# Same MATLAB AFQ cleaning as Step 5g, but for CST tracts in passthrough_cst/
print(f"{'=' * 70}")
print(f" CLEANING CST PASSTHROUGH TRACTS WITH MATLAB AFQ")
print(f"{'=' * 70}")

cst_clean_tasks = []
cst_clean_plan = build_cst_cleaning_tasks(
    subjects,
    modules=MODULES,
    module_load_matlab=MODULE_LOAD_MATLAB,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=None,
)
for msg in cst_clean_plan['skips']:
    print(f"  {msg}")
cst_clean_tasks = cst_clean_plan['tasks']

if cst_clean_tasks:
    cst_clean_results = run_parallel(cst_clean_tasks, label="CST_cleaning")
else:
    print("  Nothing to clean (all already done or no raw CST tracts).")

# Show cleaned counts
for row in collect_cst_cleaned_counts(subjects):
    s = row['subject']
    print(f"  {s['id']:<30}: {row['n_cleaned']} cleaned CST tracts")

# %% [markdown] Step 5f: Passthrough Tract Extraction
# ---
# ## Step 5f: Passthrough Tract Extraction
#
# Extract streamlines that pass through specific ROI pairs from the whole-brain tractogram.
# Uses `tckedit` with `-include` flags.
#
# Parameters from `config1.yaml → passthrough_extraction`.
# For AP/APTR: also extracts using **shrunk** (eroded) graft segments.

# %% [Step 5f] Show Passthrough Extraction Parameters
print("Passthrough extraction parameters (from config1.yaml):")
print(f"  Targets CTRL:    {pt.get('targets', {}).get('CTRL', [])}")
print(f"  Targets AP/APTR: {pt.get('targets', {}).get('AP_APTR', [])}")
print(f"  Segments CTRL:   {pt.get('segments', {}).get('CTRL', [])}")
print(f"  Segments AP:     {pt.get('segments', {}).get('AP_APTR', [])}")
print(f"  Hemispheres:     {pt.get('hemispheres', [])}")
print(f"  Cerebellum incl thalamus: {pt.get('cerebellum_include_thalamus', False)}")
print(f"  Shrink enabled:  {sk.get('enabled', False)}")
print(f"  Erosion voxels:  {sk.get('erosion_voxels', 'N/A')}")
print(f"  SEG_MODE:        {SEG_MODE}")
print(f"  To change: edit config1.yaml → passthrough_extraction / shrink_analysis")

# %% [Step 5f] Extract Passthrough Tracts — PARALLEL (both hemispheres + shrink)
# IMPORTANT: --create-shrink-segments processes ALL subjects globally, so it must
# run ONCE sequentially before the parallel extraction to avoid race conditions.
print(f"{'=' * 70}")
print(f" EXTRACTING PASSTHROUGH TRACTS (PARALLEL) — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

# --- Step 1: Create shrunk segments ONCE (sequential, affects all AP/APTR) ---
shrink_cmd = build_shrink_segments_command(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
)
if shrink_cmd is not None:
    print("  Creating shrunk graft segments (sequential, all AP/APTR subjects)...")
    shrink_result = subprocess.run(shrink_cmd, shell=True, capture_output=True, text=True)
    if shrink_result.returncode == 0:
        print("  ✅ Shrunk segments created successfully")
    else:
        print(f"  ❌ Shrunk segments FAILED (rc={shrink_result.returncode})")
        print(f"     {shrink_result.stderr[-300:]}")

# --- Step 2: Extract passthrough tracts in parallel (left + right + shrink) ---
pt_tasks = build_passthrough_extraction_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=None,
)

pt_results = run_parallel(pt_tasks, label="Passthrough extraction")

# %% [Step 5f-QC] QC: List tract counts per subject
print(f"{'=' * 70}")
print(f" PASSTHROUGH TRACT COUNTS")
print(f"{'=' * 70}")
for row in collect_passthrough_tract_counts(subjects):
    s = row['subject']
    if row['has_dir']:
        print(f"  {s['id']:<30} ({s['group']:<4}): {row['n_raw']} raw tracts")
    else:
        print(f"  {s['id']:<30} ({s['group']:<4}): NO passthrough dir")

# %% [Step 5f-QC2] QC: View sample tract for first subject
s = subjects[0]
sample_tract, roi_overlay = select_sample_tract(s, tract_pattern="*thalamus_left*.tck")
if sample_tract and roi_overlay:
    launch_gui(
        f"mrview {s['t2w']} "
        f"-overlay.load {roi_overlay} -overlay.opacity 0.4 "
        f"-tractography.load {sample_tract}"
    )
    print(f"Showing: {sample_tract.name} for {s['id']}")
elif sample_tract and not roi_overlay:
    launch_gui(f"mrview {s['t2w']} -tractography.load {sample_tract}")
    print(f"Showing: {sample_tract.name} (ROI overlay not found)")
else:
    print("No thalamus tract found for visualization")

# %% [markdown] Step 5g: MATLAB AFQ Cleaning
# ---
# ## Step 5g: MATLAB AFQ Cleaning
#
# Remove outlier streamlines using MATLAB `AFQ_removeFiberOutliers`.
# All parameters from `config1.yaml → matlab_cleaning`.

# %% [Step 5g] Show MATLAB AFQ Cleaning Parameters
print("MATLAB AFQ Cleaning Parameters (from config1.yaml → matlab_cleaning):")
print(f"  max_dist:        {mat['max_dist']}  ← Gaussian distance from core (SD)")
print(f"  max_len:         {mat['max_len']}  ← Length deviation (SD)")
print(f"  num_nodes:       {mat['num_nodes']}  ← Resampling nodes per fiber")
print(f"  min_streamlines: {mat['min_streamlines']}  ← Skip cleaning if fewer")
print(f"  max_iter:        {mat.get('max_iter', 'N/A')}  ← Max cleaning iterations")
print(f"  core_method:     {mat.get('core_method', 'N/A')}")
print(f"  matlab_bin:      {mat.get('matlab_bin', 'N/A')}")
print(f"  afq_path:        {mat.get('afq_path', 'N/A')}")
print(f"  To change: edit config1.yaml → matlab_cleaning")

# %% [Step 5g] Run MATLAB AFQ Cleaning — PARALLEL
# NOTE: MATLAB uses ~2 cores per instance => 12 subjects × 2 = 24 cores (fine for lmc02)
print(f"{'=' * 70}")
print(f" CLEANING PASSTHROUGH TRACTS WITH MATLAB AFQ (PARALLEL)")
print(f"{'=' * 70}")

clean_tasks = build_passthrough_cleaning_tasks(
    subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout_ctrl=3600,
    timeout_ap_aptr=7200,
)

clean_results = run_parallel(clean_tasks, label="MATLAB cleaning")

# %% [Step 5g-QC] QC: Compare raw vs cleaned for first subject
s = subjects[0]
raw, cleaned = select_raw_cleaned_comparison_pair(s)
if raw and cleaned:
    print(f"Comparing: {raw.name}")
    run(f"{MODULES} && tckinfo {raw} -count")
    print(f"Cleaned: {cleaned.name}")
    run(f"{MODULES} && tckinfo {cleaned} -count")
    launch_gui(f"mrview {s['t2w']} -tractography.load {cleaned}")
    print("Cleaned tract should be a tight, clean bundle.")
else:
    print("No cleaned tracts found yet.")

# %% [Step 5fg-CTRL-rerun] Re-run passthrough + cleaning for CTRL only (after ROI L/R fix)
# ────────────────────────────────────────────────────────────────────────
# After fixing isocortex & olfactory L/R labels for CTRL subjects, we must:
#   1. Delete old passthrough tracts (raw + cleaned + SF)
#   2. Re-extract passthrough tracts  (Step 5f logic, CTRL only)
#   3. Re-run MATLAB AFQ cleaning     (Step 5g logic, CTRL only)
#
# AP and APTR are NOT affected — their ROIs were verified correct.
# ────────────────────────────────────────────────────────────────────────
import shutil

ctrl_subjects = [s for s in subjects if s['group'] == 'CTRL']
print(f"{'=' * 70}")
print(f" CTRL-ONLY RERUN: PASSTHROUGH + CLEANING ({len(ctrl_subjects)} subjects)")
print(f"{'=' * 70}")

# ── Part 1: Delete old CTRL passthrough tracts ──
print(f"\n{'─' * 50}")
print(f" Part 1: Deleting old CTRL passthrough tracts")
print(f"{'─' * 50}")
for s in ctrl_subjects:
    pt_d = s['pt_dir']  # .../wholebrain/passthrough_tracts/
    if pt_d.exists():
        old_tcks = list(pt_d.glob('*.tck'))
        if old_tcks:
            for f in old_tcks:
                f.unlink()
            print(f"  🗑  {s['id']}: deleted {len(old_tcks)} old .tck files")
        else:
            print(f"  ⚪ {s['id']}: no .tck files to delete")
    else:
        print(f"  ⚪ {s['id']}: passthrough_tracts/ does not exist")

# ── Part 2: Re-extract passthrough tracts (CTRL only, left + right) ──
print(f"\n{'─' * 50}")
print(f" Part 2: Extracting passthrough tracts (CTRL × left + right)")
print(f"{'─' * 50}")

pt_ctrl_tasks = build_ctrl_passthrough_extraction_tasks(
    ctrl_subjects,
    modules=MODULES,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=None,
)

pt_ctrl_results = run_parallel(pt_ctrl_tasks, label="CTRL_passthrough")

# Show tract counts after extraction
for s in ctrl_subjects:
    pt_d = s['pt_dir']
    if pt_d.exists():
        n = len(list(pt_d.glob('*.tck')))
        print(f"  {s['id']}: {n} tracts extracted")

# ── Part 3: Re-run MATLAB AFQ cleaning (CTRL only) ──
print(f"\n{'─' * 50}")
print(f" Part 3: MATLAB AFQ cleaning (CTRL only)")
print(f"{'─' * 50}")

clean_ctrl_tasks = build_ctrl_cleaning_tasks(
    ctrl_subjects,
    modules=MODULES,
    module_load_matlab=MODULE_LOAD_MATLAB,
    repo_dir=REPO_DIR,
    derivatives_dir=DERIV_BASE,
    timeout=None,
)

clean_ctrl_results = run_parallel(clean_ctrl_tasks, label="CTRL_cleaning")

# Show cleaned counts
for row in collect_ctrl_passthrough_counts(ctrl_subjects):
    s = row['subject']
    if row['has_dir']:
        print(f"  {s['id']}: {row['n_raw']} raw, {row['n_cleaned']} cleaned")
    else:
        print(f"  {s['id']}: NO passthrough dir")

print(f"\n{'=' * 70}")
print(f" CTRL RERUN COMPLETE")
print(f"{'=' * 70}")

# %% [markdown] Step 5h: Generate Metrics Report
# ---
# ## Step 5h: Generate Tract Metrics Report (v2)
#
# Compute tract metrics and generate CSV + HTML report:
# - **Streamline count** (from cleaned passthrough tracts)
# - **SIFT2 weighted** (sum of wholebrain SIFT2 weights for matched raw streamlines)
# - **Density** (streamline count / target ROI volume in mm³)
#
# Source ROIs: CTRL → isocortex {ant, cen, pos}; AP/APTR → seg9 graft dorsal {ant, cen, pos} (shrunk)
#
# Target ROIs: L/R × {Striatum, Thalamus, Paleocortex, Cerebellum}
#
# Output: CSV (main + 3 companion) + HTML (tables + bar charts with error bars + individual points)

# %% [Step 5h] Generate Tract Metrics Report v2 (CSV + HTML)
import datetime as _dt
report_version = f"v{_dt.datetime.now().strftime('%Y%m%d')}_newmask"

print(f"{'=' * 70}")
print(f" GENERATING TRACT METRICS REPORT v2")
print(f" Version: {report_version}")
print(f"{'=' * 70}")

report_cmd = (
    f"{MODULES} && cd {REPO_DIR} && "
    f"python functions/generate_tract_metrics_report_v2.py "
    f"--version {report_version} "
    f"--output-dir {REPO_DIR / 'reports'} "
    f"--derivatives-dir {DERIV_BASE}"
)
print(f"Command:\n  {report_cmd[:300]}...")
result = run(report_cmd, timeout=7200)

reports_dir = REPO_DIR / 'reports'
check_file(reports_dir / f'tract_metrics_{report_version}.csv', 'Main CSV')
check_file(reports_dir / f'tract_metrics_count_{report_version}.csv', 'Count CSV')
check_file(reports_dir / f'tract_metrics_sift2_{report_version}.csv', 'SIFT2 CSV')
check_file(reports_dir / f'tract_metrics_density_{report_version}.csv', 'Density CSV')
check_file(reports_dir / f'tract_metrics_report_{report_version}.html', 'HTML report')

# %% [Step 5h] Preview CSV Report
import pandas as pd
csv_path = REPO_DIR / "reports" / f"tract_metrics_{report_version}.csv"
if csv_path.exists():
    df = pd.read_csv(csv_path)
    print(f"Report: {len(df)} rows\n")
    print("Columns:", list(df.columns))
    print(f"\nFirst 20 rows:")
    print(df.head(20).to_string(index=False))
else:
    print("Report CSV not found.")

# %% [markdown] Step 5i: Fixel-Based Analysis (FBA)
# ---
# ## Step 5i: Fixel-Based Analysis (FBA)
#
# Decompose FODs into discrete fixels and compute microstructural metrics.
#
# **FBA pipeline:**
# 1. **Generate fixels** — `fod2fixel` on `wm_fod_norm_group.mif` → `fixels_group/`
# 2. **Voxel AFD map** — `fixel2voxel` collapses fixel AFD to voxel space → `afd_voxel_sum.mif`
# 3. **ROI fixel metrics** — directional FD within each ROI (AP/DV/LR decomposition)
# 4. **Tract AFD sampling** — `tcksample` samples AFD along cleaned passthrough tracts
#
# **Key metrics:**
# - `fd_total` — Total Fiber Density in ROI
# - `fd_ap`, `fd_dv`, `fd_lm` — Directional FD (Anterior-Posterior, Dorsal-Ventral, Lateral-Medial)
# - `mean_afd` — Mean AFD along a tract path
#
# **Axis convention (after flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 5i-a] Generate Fixels (fod2fixel) — PARALLEL
print(f"{'=' * 70}")
print(f" GENERATING FIXELS (fod2fixel) — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

fixel_tasks = []
for s in subjects:
    DWI = s['dwi_dir']
    FIXEL_DIR = s['fixel_dir']
    fod_file = DWI / 'wm_fod_norm_group.mif'
    mask_file = DWI / 'mask_anat.mif'
    if not fod_file.exists():
        print(f"  ❌ {s['id']}: no group FOD, skipping")
        continue
    # fod2fixel requires an empty output dir — remove stale fixels first
    fixel_tasks.append((
        s['id'],
        f"{MODULES} && rm -rf {FIXEL_DIR} && "
        f"fod2fixel {fod_file} {FIXEL_DIR} "
        f"-mask {mask_file} -afd afd.mif -peak_amp peak.mif -force -quiet",
        600
    ))

fixel_results = run_parallel(fixel_tasks, label="Fixel generation")

# %% [Step 5i-b] Create Voxel AFD Maps (fixel2voxel) — PARALLEL
print(f"{'=' * 70}")
print(f" CREATING VOXEL AFD MAPS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

afd_tasks = []
for s in subjects:
    FIXEL_DIR = s['fixel_dir']
    AFD_VOXEL = s['afd_voxel']
    if (FIXEL_DIR / 'afd.mif').exists():
        afd_tasks.append((
            s['id'],
            f"{MODULES} && fixel2voxel {FIXEL_DIR}/afd.mif sum {AFD_VOXEL} -force",
            120
        ))
    else:
        print(f"  ❌ {s['id']}: no fixels, skipping")

afd_voxel_results = run_parallel(afd_tasks, label="Voxel AFD maps")

# %% [Step 5i-c] ROI-Based Fixel Metrics (Directional FD) — PARALLEL
print(f"{'=' * 70}")
print(f" COMPUTING ROI FIXEL METRICS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

roi_fixel_tasks = []
for s in subjects:
    if s['group'] == 'CTRL':
        roi_names = ['isocortex_left', 'isocortex_right',
                     'thalamus_left', 'thalamus_right',
                     'striatum_left', 'striatum_right',
                     'cerebellum_left', 'cerebellum_right']
    else:
        roi_names = ['graft_left', 'graft_right',
                     'thalamus_left', 'thalamus_right',
                     'striatum_left', 'striatum_right',
                     'palecortex_left', 'palecortex_right',
                     'cerebellum_left', 'cerebellum_right']
    fixel_csv = REPO_DIR / 'reports' / f'fixel_metrics_{s["id"]}.csv'
    # Chain all ROIs for this subject into one compound command
    cmds = [f"{MODULES} && cd {REPO_DIR} && python functions/analyze_roi_fixels.py "
            f"--subject {s['id']} --roi {roi} --output {fixel_csv} --force --derivatives-dir {DERIV_BASE}"
            for roi in roi_names]
    roi_fixel_tasks.append((s['id'], ' && '.join(cmds), 1200))

roi_fixel_results = run_parallel(roi_fixel_tasks, label="ROI fixel metrics")

# %% [Step 5i-d] Tract AFD Sampling (tcksample) — PARALLEL
import numpy as np
print(f"{'=' * 70}")
print(f" SAMPLING AFD ALONG CLEANED TRACTS — PARALLEL — {len(subjects)} SUBJECTS")
print(f"{'=' * 70}")

tcksample_tasks = []
for s in subjects:
    AFD_VOXEL = s['afd_voxel']
    if not AFD_VOXEL.exists():
        print(f"  ❌ {s['id']}: no afd_voxel_sum.mif, skipping")
        continue
    # Collect all cleaned tracts for this subject
    cleaned = sorted(glob.glob(str(s['pt_dir'] / '*_cleaned.tck'))) if s['pt_dir'].exists() else []
    if s['group'] in ['AP', 'APTR']:
        cleaned += sorted(glob.glob(str(s['pt_shrink_dir'] / '*_cleaned.tck'))) if s['pt_shrink_dir'].exists() else []
    if not cleaned:
        print(f"  ⚠️  {s['id']}: no cleaned tracts")
        continue
    # Build compound tcksample command for all tracts of this subject
    cmds = []
    for tp in cleaned:
        tn = Path(tp).stem.replace('_cleaned', '')
        afd_out = Path(tp).parent / f'{tn}_afd.txt'
        cmds.append(f"{MODULES} && tcksample {tp} {AFD_VOXEL} {afd_out} -stat_tck mean -force")
    tcksample_tasks.append((s['id'], ' && '.join(cmds), 600))

tcksample_results = run_parallel(tcksample_tasks, label="Tract AFD sampling")

# %% [Step 5i-e] Slicewise Fixel Analysis + 3-Group Comparison
print(f"{'=' * 70}")
print(f" SLICEWISE FIXEL ANALYSIS — 3-GROUP COMPARISON")
print(f"{'=' * 70}")

# Run the 3-group comparison (compares CTRL vs AP vs APTR per slice)
for hemi in ['left', 'right']:
    print(f"\n  --- {hemi.upper()} hemisphere ---")
    run(f"{MODULES} && cd {REPO_DIR} && python functions/analyze_roi_fixels.py "
        f"--compare --force --derivatives-dir {DERIV_BASE}",
        timeout=600)

print("\n--- Slicewise fixel reports ---")
for hemi in ['left', 'right']:
    for shrink_tag in ['', '_shrink']:
        fig_path = REPO_DIR / 'reports' / f'fixel_3group_comparison_{hemi}{shrink_tag}.png'
        if fig_path.exists():
            check_file(fig_path)

# %% [Step 5i-km] 3-Group Fixel K-means Direction Comparison
import sys, yaml, importlib
import numpy as np
from pathlib import Path
sys.path.insert(0, str(REPO_DIR))

# Ensure derivatives override is set for direct Python imports
os.environ['MBRAIN_DERIVATIVES_DIR'] = str(DERIV_BASE)

import functions.fixel_roi_analysis as _fri
importlib.reload(_fri)
from functions.fixel_roi_analysis import analyze_3group_kmeans

print('Step 5i-km: 3-Group Fixel K-means Direction Comparison\n')

with open(REPO_DIR / 'config1.yaml') as f:
    cfg_km = yaml.safe_load(f)

# Override derivatives paths in the loaded config
cfg_km.setdefault('paths', {})['derivatives_flip1'] = str(DERIV_BASE)
cfg_km['paths']['derivatives_dir'] = str(DERIV_BASE)

report_dir = REPO_DIR / 'reports'
combos = [(h, True) for h in ['left', 'right']]

# Pass 1: compute data to find unified y_max
print('Pass 1: computing data for all combinations...')
all_results_km = {}
global_max = 0.0
for hemi, use_shrink in combos:
    tag = f'{hemi}_shrink'
    r = analyze_3group_kmeans(
        cfg_km, hemisphere=hemi, use_shrink=True,
        output_dir=report_dir, n_clusters=2, conda_env='mbrain',
        force=False, show=False, verbose=False, plot=False)
    all_results_km[tag] = (hemi, True, r)
    if r:
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, data in r['all_group_data'].get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get('scaled_directions', []):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val

y_max_unified = max(0.05, np.ceil(global_max * 20) / 20)
print(f'  Unified y_max = {y_max_unified:.4f} (global max = {global_max:.4f})\n')

# Pass 2: generate plots with forced y_max
for tag, (hemi, use_shrink, r) in all_results_km.items():
    print(f'\n{"="*60}')
    print(f'  {tag.upper()}')
    print(f'{"="*60}')
    if not r:
        print(f'  ⚠ No data for {tag}')
        continue
    from functions.fixel_roi_analysis import (
        plot_3group_kmeans_comparison, plot_3group_component_lines
    )
    suffix = '_shrink' if use_shrink else ''
    tns = r['target_n_slices']

    fig_path = report_dir / f'fixel_3group_comparison_{hemi}{suffix}.png'
    plot_3group_kmeans_comparison(
        r['all_group_data'], hemi, tns,
        use_shrink=use_shrink, save_path=fig_path, show=True,
        y_max=y_max_unified)

    lines_path = report_dir / f'fixel_3group_components_{hemi}{suffix}.png'
    plot_3group_component_lines(
        r['all_group_data'], hemi, tns,
        use_shrink=use_shrink, save_path=lines_path, show=True,
        y_max=y_max_unified)

    print(f'  ✅ Radar: {fig_path}')
    print(f'  ✅ Lines: {lines_path}')

# %% [Step 5i-fig] Generate Fixel 3×3 Composite Figure v08 (newmask)
# Generates the 3×3 composite figure:
#   Row 1: FOD screenshots (CTRL | AP | APTR)
#   Row 2: Fixel screenshots
#   Row 3: Per-mouse radar (4 mice × DV-dominant direction, first-thalamus slice)
# Uses mrtrix_newmask derivatives. Output: reports/fixel_4x3_comparison_v08.{png,svg}

print(f"{'=' * 70}")
print(f" GENERATING FIXEL 3×3 COMPOSITE FIGURE v08 [newmask]")
print(f"{'=' * 70}")

fig_cmd = (
    f"{MODULES} && cd {REPO_DIR} && "
    f"python tmp/make_4x3_figure_v08.py"
)
print(f"Command:\n  {fig_cmd}\n")
result = run(fig_cmd, timeout=3600)

check_file(REPO_DIR / 'reports' / 'fixel_4x3_comparison_v08.png', 'Fixel 4×3 PNG')
check_file(REPO_DIR / 'reports' / 'fixel_4x3_comparison_v08.svg', 'Fixel 4×3 SVG')

# %% [markdown] QC 5i: Fixel Analysis
# ### QC Checkpoint: Fixel Analysis
#
# **Verify:**
# - Fixels generated without errors
# - AFD voxel map looks reasonable (higher in WM, low in CSF)
# - ROI metrics show expected directional patterns
# - Tract AFD values consistent across subjects in same group

# %% [QC 5i] View AFD Voxel Map — check all subjects
print('QC: Checking AFD voxel maps for all subjects\n')
found_afd = 0
missing_afd = 0
for s in subjects:
    if s['afd_voxel'].exists():
        print(f"  ✅ {s['id']}: afd_voxel_sum.mif exists")
        found_afd += 1
    else:
        print(f"  ❌ {s['id']}: afd_voxel_sum.mif NOT FOUND")
        missing_afd += 1

print(f"\n  Summary: {found_afd} found, {missing_afd} missing")

# Open first available for visual check
for s in subjects:
    if s['afd_voxel'].exists():
        print(f"\nOpening {s['id']} for visual check...")
        launch_gui(f"mrview {s['t2w']} -overlay.load {s['afd_voxel']} -overlay.opacity 0.5")
        print("AFD voxel map: should be bright in WM regions, dark in CSF/outside brain.")
        break

# %% [markdown] Step 5j: Fixel Square ROI Analysis
# ---
# ## Step 5j: Fixel Square ROI Direction Analysis
#
# Analyze the **main fiber directions** within a user-drawn square ROI.
#
# **Pipeline:**
# 1. Run `fod2fixel` masked to the square ROI → local fixel directory
# 2. Extract fixel directions (unit vectors) and AFD values
# 3. Cluster fixels into 2 populations (K-means on direction vectors)
# 4. Report **Main Direction 1** and **Main Direction 2** with mean AFD
#
# **Prerequisites:**
# - `square_roi_fixel.mif` drawn in mrview ROI editor and saved in `dwi/`
# - `wm_fod_norm_group.mif` (from Step 3b / re-run in Step 5b)
#
# **Axis convention (flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 5j] Fixel Square ROI — Generate fixels for each subject
from functions.fixel_roi_analysis import generate_fixel_peaks

print(f"{'=' * 70}")
print(f" FIXEL SQUARE ROI ANALYSIS")
print(f"{'=' * 70}")

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    DWI = s['dwi_dir']
    SQUARE_ROI_T2W = DWI / 'square_roi_fixel.mif'
    FOD_PATH = DWI / 'wm_fod_norm_group.mif'
    
    if not SQUARE_ROI_T2W.exists():
        print(f"  ⏭️ {label}: no square_roi_fixel.mif, skipping")
        continue
    if not FOD_PATH.exists():
        print(f"  ❌ {label}: no wm_fod_norm_group.mif, skipping")
        continue
    
    peaks_nii = generate_fixel_peaks(
        roi_path=SQUARE_ROI_T2W,
        fod_path=FOD_PATH,
        output_dir=DWI,
        conda_env='mbrain',
        force=True,
    )
    print(f"  ✅ {label}: fixel extraction complete")

# %% [Step 5j-2] Fixel Square ROI — K-means clustering for each subject
import numpy as np
from functions.fixel_roi_analysis import extract_peak_vectors, cluster_fixel_directions, plot_fixel_directions

print(f"{'=' * 70}")
print(f" FIXEL SQUARE ROI — K-MEANS CLUSTERING")
print(f"{'=' * 70}")

for i, s in enumerate(subjects):
    label = f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})"
    DWI = s['dwi_dir']
    PEAKS_NII = DWI / 'square_roi_peaks.nii'
    
    if not PEAKS_NII.exists():
        print(f"  ⏭️ {label}: no peaks file, skipping")
        continue
    
    all_vectors, all_peak_indices, n_peaks = extract_peak_vectors(PEAKS_NII)
    if len(all_vectors) == 0:
        print(f"  ⚠️ {label}: no valid fixel vectors")
        continue
    
    print(f"\n  {label}: {len(all_vectors)} fixel vectors")
    cluster_result = cluster_fixel_directions(
        all_vectors, all_peak_indices, n_peaks,
        n_clusters=2, verbose=True,
    )
    
    if cluster_result is not None:
        fig_path = REPO_DIR / 'reports' / f'fixel_square_roi_{s["id"]}.png'
        fig = plot_fixel_directions(
            cluster_result,
            subject_id=s['id'],
            save_path=fig_path,
            show=True,
        )
        print(f"  ✅ {s['id']}: saved {fig_path.name}")

# %% [markdown] Step 6: QC — Verify New Masks in Context
# ---
# ## Step 6: Final QC — Verify New Masks + Downstream Outputs
#
# Quick visual checks after re-running the pipeline.

# %% [Step 6a] QC: View new mask on registered DWI for each subject
for i, s in enumerate(subjects):
    b0 = s['dwi_dir'] / "b0_final.mif"
    mask = s['dwi_dir'] / "mask_anat.mif"
    if b0.exists() and mask.exists():
        print(f"[{i+1}/{len(subjects)}] {s['id']} ({s['group']})")
        launch_gui(
            f"mrview {s['t2w']} "
            f"-overlay.load {mask} -overlay.opacity 0.3 -overlay.colour 0,1,0"
        )
        answer = input("  PASS (Enter) or FAIL (f): ").strip().lower()
        if answer == 'f':
            print(f"  ⚠️ {s['id']} FAILED final QC — may need re-editing")
    else:
        print(f"  ⏭️ {s['id']}: missing b0_final or mask_anat")

# %% [Step 6b] QC: Summary table — mask sizes
print(f"{'=' * 80}")
print(f" {'Subject':<30} {'Group':<6} {'Old mask MB':>12} {'New mask MB':>12} {'Ratio':>8}")
print(f"{'=' * 80}")

for s in subjects:
    old = s['dwi_dir'] / f"mask_old_{BACKUP_SUFFIX}.mif"
    new = s['dwi_dir'] / "mask.mif"
    
    old_sz = old.stat().st_size / (1024*1024) if old.exists() else 0
    new_sz = new.stat().st_size / (1024*1024) if new.exists() else 0
    ratio = new_sz / old_sz if old_sz > 0 else float('inf')
    
    flag = "" if 0.5 < ratio < 2.0 else " ⚠️"
    print(f" {s['id']:<30} {s['group']:<6} {old_sz:>10.3f}  {new_sz:>10.3f}  {ratio:>6.2f}x{flag}")

print(f"{'=' * 80}")

# %% [markdown] Step 7: Final Visualization & Summary
# ---
# ## Step 7: Final Visualization & Summary

# %% [Step 7a] Pipeline Summary Table — ALL subjects
print(f"\n{'=' * 80}")
print(f" PIPELINE COMPLETE: {len(subjects)} SUBJECTS")
print(f"{'=' * 80}\n")

for s in subjects:
    pt_d = s['pt_dir']
    raw_tracts = glob.glob(str(pt_d / '*.tck')) if pt_d.exists() else []
    raw_tracts = [t for t in raw_tracts if '_cleaned' not in t]
    cleaned_tracts = glob.glob(str(pt_d / '*_cleaned.tck')) if pt_d.exists() else []
    
    wb_ok = "✅" if s['wb_tck'].exists() else "❌"
    sift2_ok = "✅" if s['sift2_weights'].exists() else "❌"
    
    print(f"  {s['id']:<30} ({s['group']:<4})")
    print(f"    WB tractogram: {wb_ok}  |  SIFT2: {sift2_ok}")
    print(f"    Raw tracts: {len(raw_tracts):>3}  |  Cleaned: {len(cleaned_tracts):>3}")
    
    if s['group'] in ['AP', 'APTR']:
        shrink_d = s['pt_shrink_dir']
        shrink_raw = glob.glob(str(shrink_d / '*.tck')) if shrink_d.exists() else []
        shrink_raw = [t for t in shrink_raw if '_cleaned' not in t]
        shrink_cleaned = glob.glob(str(shrink_d / '*_cleaned.tck')) if shrink_d.exists() else []
        print(f"    Shrunk raw: {len(shrink_raw):>3}  |  Shrunk cleaned: {len(shrink_cleaned):>3}")

print(f"\nReport: reports/tract_metrics_{report_version}.csv")
print(f"HTML:   reports/tract_metrics_report_{report_version}.html")
print(f"\n{'=' * 80}")

# %% [Step 7b] Final Visualization — Cleaned Tracts for first subject per group
for group_name in ['CTRL', 'AP', 'APTR']:
    group_subjects = [s for s in subjects if s['group'] == group_name]
    if not group_subjects:
        continue
    s = group_subjects[0]
    pt_d = s['pt_dir']
    
    if s['group'] == 'CTRL':
        seg_roi = s['base'] / 'atlas' / 'rois' / 'fod_space' / 'roi_isocortex_left.nii.gz'
        target_roi = s['base'] / 'atlas' / 'rois' / 'fod_space' / 'roi_thalamus_left.nii.gz'
    else:
        seg_roi = s['base'] / 'manual_roi' / 'fod_space' / 'roi_Graft_left_flip_fod.mif'
        target_roi = s['base'] / 'manual_roi' / 'fod_space' / 'roi_Thalamus_left_flip_fod.mif'
    
    cleaned_left = sorted(glob.glob(str(pt_d / '*left*_cleaned.tck')))
    cmd_viz = (f"mrview {s['t2w']} "
               f"-overlay.load {seg_roi} -overlay.opacity 0.3 "
               f"-overlay.load {target_roi} -overlay.opacity 0.3")
    for t in cleaned_left[:6]:
        cmd_viz += f" -tractography.load {t}"
    
    print(f"Opening final visualization for {s['id']} ({group_name})...")
    launch_gui(cmd_viz)
    print(f"  Showing left-hemisphere cleaned tracts for {s['id']}")

# %% [markdown] Bar-plots: APTR SIFT2 Density
# ---
# ## Bar-plots: APTR SIFT2 Density (Left Hemisphere & Central Graft)
#
# Bar plots showing SIFT2 density (connectivity / seed_voxels × 1000) for all
# APTR mice. Individual data points overlaid on group means.

# %% [Bar-plots] APTR Density — Left hemisphere & Central summary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Load metrics CSV (cleaned tracts)
csv_path = REPO_DIR / "reports" / f"tract_metrics_{report_version}.csv"
if not csv_path.exists():
    print(f"❌ Report CSV not found: {csv_path}")
    print("  Run Step 5h first to generate the report.")
else:
    df_plot = pd.read_csv(csv_path)
    
    # Filter APTR, left hemisphere only
    aptr_left = df_plot[(df_plot['group'] == 'APTR') & (df_plot['hemisphere'] == 'left')].copy()
    targets_order = ['olfactory_cortex', 'striatum', 'thalamus', 'cerebellum']
    target_labels = {'thalamus': 'Thalamus', 'striatum': 'Striatum',
                     'olfactory_cortex': 'Paleocortex', 'cerebellum': 'Cerebellum'}
    segments_order = ['anterior', 'central', 'posterior']
    seg_colors = {'anterior': '#E24A33', 'central': '#348ABD', 'posterior': '#988ED5'}
    
    plot_subjects = aptr_left['subject'].unique()
    marker_list = ['o', 's', 'D', '^']
    
    import matplotlib.gridspec as gridspec
    from matplotlib.image import imread as mpl_imread
    
    fig = plt.figure(figsize=(20, 14))
    fig.suptitle('APTR — SIFT2 Density (Left Hemisphere, Cleaned Tracts)',
                 fontsize=14, fontweight='bold')
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3,
                           height_ratios=[1, 0.9])
    
    # Row 0: 3 segment bar subplots
    ax_ant = fig.add_subplot(gs[0, 0])
    draw_single_seg_panel(
        ax_ant,
        aptr_left,
        'anterior',
        seg_colors['anterior'],
        'Anterior Graft → Destinations',
        targets_order=targets_order,
        target_labels=target_labels,
        plot_subjects=plot_subjects,
        marker_list=marker_list,
    )
    ax_ant.set_ylabel('SIFT2 Density\n(conn / seed_vox × 1000)', fontsize=10)
    
    ax_cen = fig.add_subplot(gs[0, 1])
    draw_single_seg_panel(
        ax_cen,
        aptr_left,
        'central',
        seg_colors['central'],
        'Central Graft → Destinations',
        targets_order=targets_order,
        target_labels=target_labels,
        plot_subjects=plot_subjects,
        marker_list=marker_list,
    )
    ax_cen.set_ylabel('SIFT2 Density', fontsize=10)
    
    ax_pos = fig.add_subplot(gs[0, 2])
    draw_single_seg_panel(
        ax_pos,
        aptr_left,
        'posterior',
        seg_colors['posterior'],
        'Posterior Graft → Destinations',
        targets_order=targets_order,
        target_labels=target_labels,
        plot_subjects=plot_subjects,
        marker_list=marker_list,
    )
    ax_pos.set_ylabel('SIFT2 Density', fontsize=10)
    
    # Row 1: mrview screenshots by destination
    destinations = ['paleocortex', 'striatum', 'thalamus']
    dest_labels = {'paleocortex': 'Paleocortex', 'striatum': 'Striatum', 'thalamus': 'Thalamus'}
    
    for col_idx, dest in enumerate(destinations):
        ax = fig.add_subplot(gs[1, col_idx])
        img_path = REPO_DIR / 'reports' / f'mrview_most_{dest}.png'
        title = f'{dest_labels[dest]}'
        if img_path.exists():
            img = mpl_imread(str(img_path))
            ax.imshow(img)
            ax.set_title(title, fontsize=10, fontweight='bold')
        else:
            ax.text(0.5, 0.5, f'Screenshot not yet saved:\n{img_path.name}',
                    ha='center', va='center', fontsize=11, color='gray',
                    transform=ax.transAxes)
            ax.set_title(title, fontsize=10, fontweight='bold', color='gray')
        ax.set_xticks([])
        ax.set_yticks([])
    
    save_path = REPO_DIR / "reports" / "barplot_aptr_density_left_newmask.png"
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()

# %% [markdown] QC Summary & Final Screenshot
# ---
# ## QC Summary
#
# | QC  | Check | Result |
# |-----|-------|--------|
# | 1c  | Brain mask covers entire brain (atlas-based), no holes | ☐ PASS / ☐ FAIL |
# | 2   | DWI-T2w registration — b0 brain edges match T2w in all 3 planes | ☐ PASS / ☐ FAIL |
# | 4   | FOD quality (directional glyphs in WM, small in ventricles) | ☐ PASS / ☐ FAIL |
# | 5   | ROIs in FOD space (aligned on FOD grid) | ☐ PASS / ☐ FAIL |
# | 6   | Segmentation (no overlaps, tiles seed ROI, SEG_MODE=9 for AP/APTR) | ☐ PASS / ☐ FAIL |
# | 7   | Whole-brain tractogram (fills brain volume) | ☐ PASS / ☐ FAIL |
# | 8   | Passthrough tracts (connect correct regions) | ☐ PASS / ☐ FAIL |
# | 9   | Cleaned tracts (tight bundles, outliers removed) | ☐ PASS / ☐ FAIL |
# | 10  | Tract metrics report generated (CSV + HTML) | ☐ PASS / ☐ FAIL |
# | 11  | AFD voxel map reasonable; ROI fixel metrics consistent | ☐ PASS / ☐ FAIL |
# | 11f | Fixel square ROI directions cluster correctly | ☐ PASS / ☐ FAIL |
#
# **Processed by:** _Your name_
# **Date:** _Date_
# **Notes:** _Any issues or observations_
