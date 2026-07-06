# %% [markdown] MBRAIN Pipeline: Processing a New Mouse Subject
# # MBRAIN Pipeline: Processing a New Mouse Subject
# This notebook walks through **every step** of the MBRAIN tractography pipeline.
#
# **How to use this notebook:**
# 1. Fill in the subject info in Section 0
# 2. Run each cell sequentially
# 3. At each QC checkpoint, inspect the output
#
# **Pipeline Overview:**
# ```
# 0. Setup & Config
# 1. T2w + DWI Preprocessing (flip1 transformation, denoise, register)
# 2. ROI Preparation (CTRL: Atlas | AP/APTR: Manual ROI flip)
# 3. FOD Estimation (response functions, multi-tissue CSD, normalization)
# 4. ROI Resampling to FOD space
# 5. Segmentation (CTRL: 5 A-P | AP/APTR: 9 segments)
# 6. Whole-Brain Tractography (1M streamlines)
# 7. SIFT2 Weights (per-streamline connectivity weights)
# 8. Passthrough Tract Extraction (tckedit with -include)
# 9. MATLAB AFQ Cleaning
# 10. Metrics Report Generation
# 11. Final Visualization & Summary
# ```
#
# ---

# %% [markdown] Step 0: Setup & Configuration
# ## Step 0: Setup & Configuration
#
# Load the config and set up all paths. **Edit the cell below** with your subject info.

# %% [Step 0] Subject Info — EDIT THIS CELL
# ============================================================
# EDIT THIS CELL: Fill in your subject information
# ============================================================

SUBJECT_ID = "sub-APtestmouseNmM"       # Must match config1.yaml groups → id
SESSION_ID = "ses-240911111853"           # Must match config1.yaml groups → session
GROUP      = "AP"                  # "CTRL", "AP", or "APTR"

# For AP/APTR only: segmentation mode
SEG_MODE = 9  # 6 (3 A-P x 2 D-V) or 9 (3 A-P x 3 D-V) — recommended: 9

print(f"{'=' * 60}")
print(f"  SUBJECT_ID : {SUBJECT_ID}")
print(f"  SESSION_ID : {SESSION_ID}")
print(f"  GROUP      : {GROUP}")
print(f"  SEG_MODE   : {SEG_MODE}")
print(f"{'=' * 60}")

# %% [Step 0] Load Config & Derive Paths
import os, sys, subprocess, yaml
from pathlib import Path
from IPython.display import display, HTML, Image, Markdown

# Add project root to path
REPO_DIR = Path("/Users/glerma/soft/mbrain")
sys.path.insert(0, str(REPO_DIR))
sys.path.insert(0, str(REPO_DIR / "functions"))

# Load config
from functions.config_loader import load_config
cfg = load_config()

# Clean subject and session IDs for BIDS paths
from utils.file_utils import clean_subject_id, extract_session_info
cleaned_subject = clean_subject_id(SUBJECT_ID.replace("sub-", ""))
_, cleaned_session = extract_session_info(SESSION_ID)

# Derived paths
DERIV = Path(cfg['paths']['derivatives_dir'])
BIDS  = Path(cfg['paths']['bids_dir'])
ATLAS = Path(cfg['paths']['atlas_dir'])
BASE  = DERIV / SUBJECT_ID / SESSION_ID
DWI   = BASE / "dwi"
BIDS_SUBJ_DIR = BIDS / f"sub-{cleaned_subject}" / f"ses-{cleaned_session}"
BIDS_DWI  = BIDS_SUBJ_DIR / "dwi"
BIDS_ANAT = BIDS_SUBJ_DIR / "anat"

# Tractography params from config
tck = cfg['tractography']
mat = cfg['matlab_cleaning']

print(f"Subject:    {SUBJECT_ID} → sub-{cleaned_subject}")
print(f"Session:    {SESSION_ID} → ses-{cleaned_session}")
print(f"Group:      {GROUP}")
print(f"Base path:  {BASE}")
print(f"BIDS DWI:   {BIDS_DWI}")
print(f"BIDS Anat:  {BIDS_ANAT}")
print(f"Config loaded from: {REPO_DIR / 'config1.yaml'}")

# %% [markdown] Parameter Review (all from config1.yaml)
# ### Parameter Review (all from config1.yaml)
#
# All parameters below are read from `config1.yaml`. To change any value,
# edit that file (not this notebook) and re-run the setup cell.
#
# **To edit:** `open config1.yaml` in VS Code, find the section, save, then re-run Step 0.

# %% [Step 0] Parameter Review — All Config Values
# ============================================================
# PARAMETER REVIEW — all values read from config1.yaml
# If you want to change any parameter, edit config1.yaml and re-run this cell.
# ============================================================
wb  = cfg.get('wholebrain_tractography', {})
pt  = cfg.get('passthrough_extraction', {})
mat = cfg['matlab_cleaning']
sk  = cfg.get('shrink_analysis', {})

print(f"{'=' * 70}")
print(f" PIPELINE PARAMETERS (from config1.yaml)")
print(f" To change: edit config1.yaml → save → re-run Step 0")
print(f"{'=' * 70}")

print(f"\n--- Tractography (config1.yaml → tractography:) ---")
print(f"  algorithm:   {tck['algorithm']}")
print(f"  select:      {tck['select']}  (streamlines per ROI-to-ROI pair)")
print(f"  seeds:       {tck.get('seeds', 'N/A')}  (max seeds before giving up)")
print(f"  step_size:   {tck['step_size']} mm")
print(f"  angle:       {tck['angle']}°  (max angle between steps)")
print(f"  cutoff:      {tck['cutoff']}  (FOD amplitude stopping criterion)")
print(f"  min_length:  {tck['min_length']} mm")
print(f"  max_length:  {tck['max_length']} mm")

print(f"\n--- Whole-Brain Tractography (config1.yaml → wholebrain_tractography:) ---")
print(f"  num_streamlines: {wb.get('num_streamlines', 'N/A'):,}")
print(f"  (uses same step/angle/cutoff/length as above)")

print(f"\n--- Passthrough Extraction (config1.yaml → passthrough_extraction:) ---")
print(f"  targets CTRL:    {pt.get('targets', {}).get('CTRL', [])}")
print(f"  targets AP/APTR: {pt.get('targets', {}).get('AP_APTR', [])}")
print(f"  hemispheres:     {pt.get('hemispheres', [])}")
print(f"  segments CTRL:   {pt.get('segments', {}).get('CTRL', [])}")
print(f"  segments AP:     {pt.get('segments', {}).get('AP_APTR', [])}")
print(f"  cerebellum incl thalamus: {pt.get('cerebellum_include_thalamus', False)}")
print(f"  seg_mode:        {pt.get('seg_mode', 'N/A')}")

print(f"\n--- MATLAB AFQ Cleaning (config1.yaml → matlab_cleaning:) ---")
print(f"  max_dist:        {mat['max_dist']}  (Gaussian distance from core, in SD)")
print(f"  max_len:         {mat['max_len']}  (length deviation, in SD)")
print(f"  num_nodes:       {mat['num_nodes']}  (resampling nodes per fiber)")
print(f"  min_streamlines: {mat['min_streamlines']}  (skip cleaning if fewer)")
print(f"  max_iter:        {mat.get('max_iter', 'N/A')}  (cleaning iterations)")
print(f"  core_method:     {mat.get('core_method', 'N/A')}")
print(f"  matlab_bin:      {mat.get('matlab_bin', 'N/A')}")

print(f"\n--- Shrink Analysis (config1.yaml → shrink_analysis:) ---")
print(f"  enabled:         {sk.get('enabled', False)}")
print(f"  erosion_voxels:  {sk.get('erosion_voxels', 'N/A')}")
print(f"  primary for AP:  {sk.get('primary_for_ap_aptr', False)}")

print(f"\n--- Segmentation ---")
print(f"  CTRL:    5 A-P segments (slicewise isocortex)")
print(f"  AP/APTR: {SEG_MODE} segments (set above as SEG_MODE)")

print(f"\n{'=' * 70}")
print("Review the values above. If anything needs changing, edit config1.yaml.")
print(f"{'=' * 70}")

# %% [Step 0] Helper Functions (run, check_file, mrinfo)
# Helper: run a shell command and print output
def run(cmd, check=True, timeout=None):
    """Run a shell command, print it, and return the result."""
    print(f"\n>>> {cmd}\n")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print("STDERR:", result.stderr.strip())
    if check and result.returncode != 0:
        print(f"⚠️ Command failed with exit code {result.returncode}")
    return result

# Helper: launch a GUI app (mrview) without blocking the notebook
def launch_gui(cmd):
    """Launch a GUI command (mrview, etc.) in the background without blocking."""
    print(f"\n>>> {cmd}\n")
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("(launched in background)")


# Helper: check if file exists and show info
def check_file(path, label=""):
    p = Path(path)
    if p.exists():
        size = p.stat().st_size / (1024*1024)
        print(f"✅ {label or p.name}: {size:.1f} MB")
        return True
    else:
        print(f"❌ MISSING: {label or p.name} ({p})")
        return False

# Helper: run mrinfo
def mrinfo(path):
    run(f"mrinfo {path}")

print("Helpers loaded ✅")

# %% [markdown] Step 0b: Flywheel to BIDS Conversion
# ## Step 0b: Flywheel to BIDS Conversion
#
# If your subject's raw data is still in Flywheel format (`local/flywheel/fMRI/SUBJECTS/`),
# run this step to convert it to BIDS format before proceeding.
#
# **Flywheel layout:**
# ```
# local/flywheel/fMRI/SUBJECTS/{SubjectName}/SESSIONS/{SessionName}/ACQUISITIONS/
#     ├── DTI_EPI_seg_30dir_sat/FILES/   → .nii.gz, .bval, .mvec, .json
#     ├── 4_T2w_lowres/FILES/            → .nii.gz, .json
#     ├── 5_T2w_highres/FILES/           → .nii.gz, .json
#     ├── T2star_FID_EPI.../FILES/       → .nii.gz, .json
#     └── 1_Localizer/FILES/             → .nii.gz, .json
# ```
#
# **BIDS output:**
# ```
# data/bids/sub-{ID}/ses-{SessionID}/
#     ├── anat/   → _acq-lowres_T2w.nii.gz, _acq-highres_T2w.nii.gz, _T2starw.nii.gz, ...
#     └── dwi/    → _dwi.nii.gz, _dwi.bval, _dwi.bvec, _dwi.json
# ```
#
# **Skip this cell if your subject is already in `data/bids/`.**

# %% [Step 0b] Flywheel → BIDS Conversion
# ============================================================
# Flywheel → BIDS Conversion (uses utils/bids_converter.py)
# Set DRY_RUN=True to preview, False to actually copy files
# Skip this cell if subject is already in data/bids/
# ============================================================
from utils.bids_converter import FlywheelToBIDSConverter

DRY_RUN = False  # ← Change to False to actually copy files

FW_ROOT = REPO_DIR / "local" / "flywheel" / "fMRI"
raw_subject_name = SUBJECT_ID.replace("sub-", "")
FW_SUBJECT = FW_ROOT / "SUBJECTS" / raw_subject_name

if not FW_SUBJECT.exists():
    print(f"⚠️  Flywheel source not found: {FW_SUBJECT}")
    print(f"   Available: {[d.name for d in (FW_ROOT / 'SUBJECTS').iterdir() if d.is_dir() and not d.name.startswith('.')]}")
    print(f"   If already in BIDS, skip this cell.")
else:
    converter = FlywheelToBIDSConverter(flywheel_root=FW_ROOT, bids_root=BIDS)
    converter._convert_subject(FW_SUBJECT, dry_run=DRY_RUN)

    if DRY_RUN:
        print("\n🔍 DRY RUN — no files copied. Set DRY_RUN = False to convert.")
    else:
        print(f"\n✅ Done! Files copied: {converter.stats['files_copied']}")
        if converter.stats['errors']:
            print(f"⚠️ Errors: {converter.stats['errors']}")

# %% [Step 0b] QC: Verify BIDS Conversion Output
# 🔍 QC: Verify BIDS conversion output
print("--- BIDS output verification ---")
run(f"find {BIDS_SUBJ_DIR} -type f -not -name '.*' | sort", check=False)

print("\n--- Critical files ---")
bids_prefix = f"sub-{cleaned_subject}_ses-{cleaned_session}"
check_file(BIDS_DWI / f"{bids_prefix}_dwi.nii.gz", "DWI NIfTI")
check_file(BIDS_DWI / f"{bids_prefix}_dwi.bval", "DWI bval")
check_file(BIDS_DWI / f"{bids_prefix}_dwi.bvec", "DWI bvec")
check_file(BIDS_ANAT / f"{bids_prefix}_acq-highres_T2w.nii.gz", "T2w highres")
check_file(BIDS_ANAT / f"{bids_prefix}_acq-lowres_T2w.nii.gz", "T2w lowres")

# %% [Step 0b] Create Output Dirs & Verify BIDS Source
# Create output directories
os.makedirs(DWI, exist_ok=True)
print(f"Output directory ready: {DWI}")

# Verify BIDS source files exist
print("\n--- Checking BIDS source files ---")
print("\nAnatomy:")
run(f"ls -la {BIDS_ANAT}/", check=False)
print("\nDWI:")
run(f"ls -la {BIDS_DWI}/", check=False)

# %% [markdown] Step 1: BIDS to MRtrix Preprocessing (T2w + DWI)
# ---
# ## Step 1: BIDS to MRtrix Preprocessing (T2w + DWI)
#
# This step applies **flip1 orientation correction** and **full DWI preprocessing**.
# All operations are performed by `batch_process_flip1.py --steps 2,4`.
#
# ### Step 2: T2w Flip1 Transformation
# | Sub-step | Operation | Output File |
# |----------|-----------|-------------|
# | 2a | Axis permutation `(2,1,0)` + strides `1,2,3` | `tmp_t2w_perm.mif` (temp) |
# | 2b | Flip axis 1 (Y) to match Allen atlas orientation | `t2w_highres_flip.mif` |
#
# ### Step 4: Full DWI Preprocessing (14 sub-steps)
# | Sub-step | Operation | Command | Output File |
# |----------|-----------|---------|-------------|
# | 4a | **NIfTI→MIF conversion** with gradients | `mrconvert` + `-fslgrad` | `dwi.mif` |
# | 4b | **Axis permutation** `(2,1,0,3)` for 4D data | `mrconvert -axes` | `dwi_perm.mif` (temp) |
# | 4c | **Flip axis 1** (Y) — same transform as T2w | `mrtransform -flip 1` | `dwi_flip.mif` |
# | 4d | **Gradient check** — finds optimal gradient scheme | `dwigradcheck` + `mrconvert -grad` | `dwi_flip_gradcheck.mif` |
# | 4e | **MP-PCA denoising** — removes thermal noise (5×5×5 kernel, Exp2 estimator) | `dwidenoise` | `dwi_flip_gradcheck_denoise.mif` |
# | 4f | **Brain mask** from denoised DWI | `dwi2mask` | `mask.mif` |
# | 4g | **N4 bias field correction** via ANTs (B-spline `[200,3]`, convergence `[10,1e-6]`, sigma=4) | `dwibiascorrect ants` | `dwi_flip_gradcheck_denoise_bias.mif` |
# | 4h | **Extract mean low-b** (b<50) for registration reference | `dwiextract -bzero` + `mrmath mean` | `b0_bias.mif` |
# | 4i | **Prepare for ANTs** — convert T2w+b0 to NIfTI, mask b0 | `mrconvert` + `mrcalc` | `t2w_highres_flip.nii.gz`, `b0_bias_brain.nii.gz` |
# | 4j | **ANTs rigid registration** — align DWI (b0) to T2w anatomical space | `antsRegistrationSyN.sh -t r` | `dwi2anat_0GenericAffine.mat` |
# | 4k | **Convert transform** — ANTs `.mat` → MRtrix `.txt` | `ConvertTransformFile` + `transformconvert` | `dwi2anat_mrtrix.txt` |
# | 4l | **Apply transform** to full 4D DWI | `mrtransform -linear` | `dwi_flip_gradcheck_denoise_bias_anat.mif` |
# | 4m | **Extract final b0** reference from registered DWI | `dwiextract -bzero` + `mrmath mean` | `b0_final.mif` |
# | 4n | **Final brain mask** in anatomical (T2w) space | `dwi2mask` | `mask_anat.mif` |
#
# **Why this order matters:**
# - Denoising (4e) must happen before bias correction (4g) to preserve noise statistics
# - Mask (4f) is needed for bias correction
# - Registration (4j) needs the masked b0 and T2w in NIfTI format
# - Final mask (4n) is computed from registered data so it's in the same space as all downstream analysis

# %% [Step 1] Run Preprocessing (T2w Flip + DWI Pipeline)
# ============================================================
# BIDS to MRtrix Preprocessing (T2w + DWI)
# ============================================================

print("Running BIDS → MRtrix preprocessing...")
print("This includes T2w flipping and full DWI preprocessing pipeline.")
print("Steps: T2w flip1 transformation + DWI convert/flip/denoise/bias-correct/register")
print()

# Run the preprocessing using batch_process_flip1.py
# This script handles all the flip1 pipeline steps automatically
#
# OPTIONS for batch_process_flip1.py:
#   --subject SUBJECT_ID    : Process single subject (required if not using --group)
#   --group GROUP           : Process all subjects in group (AP, APTR, or CTRL)
#   --steps STEPS           : Comma-separated steps to run (default: 1,2,3,4,5,6,7,8)
#                           : 1=create folders, 2=flip T2w, 3=flip ROIs, 4=preprocess DWI,
#                           : 5=generate FOD, 6=resample ROIs, 7=create segments, 8=tractography
#   --dry-run               : Preview what would be done without executing
#   --force                 : Force reprocessing even if files exist
#
# EXAMPLES:
#   python functions/batch_process_flip1.py --subject sub-Apallial20R1F
#   python functions/batch_process_flip1.py --group AP --dry-run
#   python functions/batch_process_flip1.py --group CTRL --steps 2,4
#   python functions/batch_process_flip1.py --subject sub-129-921-R1F --force

result = run(
    f"cd {REPO_DIR} && micromamba run -n mbrain python functions/batch_process_flip1.py "
    f"--subject {SUBJECT_ID} --steps 2,4 --force",
    timeout=1800
)
success = result.returncode == 0

# Check results
if success:
    print("✅ Preprocessing completed successfully!")
else:
    print("❌ Preprocessing failed!")
    print(f"Return code: {result.returncode}")
    if result.stderr:
        print(f"STDERR: {result.stderr}")

# Verify outputs (same for both options)
print("\n--- Verifying preprocessing outputs ---")
check_file(DWI / "t2w_highres_flip.mif", "T2w high-res (flipped)")
check_file(DWI / "dwi_flip_gradcheck_denoise_bias_anat.mif", "Final preprocessed DWI")
check_file(DWI / "mask_anat.mif", "Brain mask (anat space)")
check_file(DWI / "b0_final.mif", "Final b0 reference")

# %% [Step 1] Verify ALL Intermediate Preprocessing Files
# ============================================================
# Detailed verification of every intermediate file from Step 1.
# This checks the full chain so you can spot where a failure occurred.
# ============================================================
import datetime as _dt

print("=" * 75)
print(" PREPROCESSING CHAIN — ALL INTERMEDIATE FILES")
print("=" * 75)

# Define the complete chain with expected pipeline order
_chain = [
    # (filename, substep, description)
    ("t2w_highres_flip.mif",                      "2",  "T2w flip1 (axis permute + flip Y)"),
    ("dwi.mif",                                    "4a", "NIfTI→MIF conversion with gradients"),
    ("dwi_flip.mif",                               "4c", "DWI axis permute + flip Y"),
    ("dwi_flip_gradcheck.mif",                     "4d", "Gradient direction check/correction"),
    ("noise.mif",                                  "4e", "Noise map from MP-PCA denoising"),
    ("dwi_flip_gradcheck_denoise.mif",             "4e", "MP-PCA denoised DWI"),
    ("mask.mif",                                   "4f", "Brain mask (pre-registration)"),
    ("dwi_flip_gradcheck_denoise_bias.mif",        "4g", "N4 bias field corrected DWI"),
    ("b0_bias.mif",                                "4h", "Mean low-b (pre-registration)"),
    ("t2w_highres_flip.nii.gz",                    "4i", "T2w NIfTI for ANTs"),
    ("b0_bias_brain.mif",                          "4i", "Masked b0 (MIF)"),
    ("b0_bias_brain.nii.gz",                       "4i", "Masked b0 NIfTI for ANTs"),
    ("dwi2anat_Warped.nii.gz",                     "4j", "ANTs rigid: b0 warped to T2w space"),
    ("dwi2anat_InverseWarped.nii.gz",              "4j", "ANTs rigid: T2w warped to DWI space"),
    ("dwi2anat_affine.txt",                        "4k", "ANTs affine (text)"),
    ("dwi2anat_mrtrix.txt",                        "4k", "MRtrix transform (converted)"),
    ("dwi_flip_gradcheck_denoise_bias_anat.mif",   "4l", "FINAL: DWI in anatomical (T2w) space"),
    ("b0_final.mif",                               "4m", "Mean low-b (post-registration)"),
    ("mask_anat.mif",                              "4n", "Brain mask in anatomical space"),
]

_ok = 0
_missing = 0
_prev_time = None
_order_ok = True
# noise.mif is written by dwidenoise together with the denoised output,
# so sub-minute timestamp differences between these two are expected.
_timestamp_skip_pairs = {("noise.mif", "dwi_flip_gradcheck_denoise.mif")}
_prev_fname = None
for _fname, _step, _desc in _chain:
    _p = DWI / _fname
    if _p.exists():
        _sz_bytes = _p.stat().st_size
        # Show KB for small files, MB for large
        if _sz_bytes < 1024 * 1024:
            _sz_str = f"{_sz_bytes / 1024:>7.1f} KB"
        else:
            _sz_str = f"{_sz_bytes / (1024 * 1024):>7.1f} MB"
        _mt = _dt.datetime.fromtimestamp(_p.stat().st_mtime)
        _ts = _mt.strftime("%H:%M:%S")
        # Check temporal ordering (skip known simultaneous outputs)
        _flag = ""
        if _prev_time is not None and _mt < _prev_time:
            _pair = (_prev_fname, _fname)
            if _pair not in _timestamp_skip_pairs:
                _flag = " ⚠️ OLDER than previous!"
                _order_ok = False
        _prev_time = _mt
        _prev_fname = _fname
        print(f"  ✅ [{_step:>3}] {_fname:<50} {_sz_str}  ({_ts}){_flag}")
        _ok += 1
    else:
        print(f"  ❌ [{_step:>3}] {_fname:<50} MISSING — {_desc}")
        _missing += 1

print(f"\n  Files: {_ok}/{len(_chain)} present, {_missing} missing")
if _order_ok and _missing == 0:
    print("  ✅ All files present and in correct temporal order.")
else:
    if _missing:
        print(f"  ⚠️ {_missing} file(s) missing — check step output above for errors.")
    if not _order_ok:
        print("  ⚠️ Some files are out of temporal order — may indicate stale outputs.")

# Quick sanity checks on key files
print(f"\n--- Quick sanity checks ---")
# 1. DWI should be larger than b0 (4D vs 3D)
_dwi_final = DWI / "dwi_flip_gradcheck_denoise_bias_anat.mif"
_b0_final  = DWI / "b0_final.mif"
if _dwi_final.exists() and _b0_final.exists():
    _ratio = _dwi_final.stat().st_size / _b0_final.stat().st_size
    print(f"  DWI/b0 size ratio: {_ratio:.0f}x {'✅' if _ratio > 5 else '⚠️ unexpectedly low'}")

# 2. Mask should be much smaller than DWI (binary)
_mask = DWI / "mask_anat.mif"
if _mask.exists() and _dwi_final.exists():
    _pct = (_mask.stat().st_size / _dwi_final.stat().st_size) * 100
    print(f"  Mask/DWI size ratio: {_pct:.2f}% {'✅' if _pct < 5 else '⚠️ mask seems large'}")

# 3. Check that DWI dimensions look reasonable
_mrinfo_result = subprocess.run(
    f"mrinfo {_dwi_final} -size -spacing -datatype",
    shell=True, capture_output=True, text=True
)
if _mrinfo_result.returncode == 0:
    _lines = _mrinfo_result.stdout.strip().split('\n')
    print(f"  Final DWI dimensions: {_lines[0] if _lines else '?'}")
    if len(_lines) > 1:
        print(f"  Voxel spacing:        {_lines[1]}")
    if len(_lines) > 2:
        print(f"  Data type:            {_lines[2]}")

print(f"\n{'=' * 75}")

# %% [markdown] QC 1: T2w Orientation
# ### QC Checkpoint 1: T2w Orientation
#
# Run the mrview command below, then inspect in mrview.
#
# **Verify:**
# - Dorsal (top of brain) is UP
# - Ventral (bottom) is DOWN
# - Anterior (nose) is to the LEFT in sagittal view (config1.yaml: AXIS 0 HIGH index = anterior)
# - Left hemisphere appears on anatomical LEFT

# %% [QC 1] View T2w Orientation in mrview
# Open subject T2w and Allen atlas template side by side for orientation comparison
ALLEN_TEMPLATE = Path(cfg['paths']['allen_template'])
print("Opening subject T2w and Allen atlas template in separate mrview windows.")
print("Compare orientations — they are NOT registered, just check that:")
print("  - Dorsal/Ventral, Anterior/Posterior, Left/Right match between the two")
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif")
launch_gui(f"mrview {ALLEN_TEMPLATE}")
print("\n👆 Two mrview windows should be open:")
print("   1) Subject T2w (flipped)")
print("   2) Allen CCF average template")
print("   Verify BOTH have the same orientation:")
print("     - Dorsal (top of brain) is UP")
print("     - Ventral (bottom) is DOWN")
print("     - Anterior (nose/olfactory) is to the LEFT in sagittal view")
print("     - Left hemisphere appears on anatomical LEFT")

# %% [markdown] QC 1 Result
# **QC 1 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] QC 1b: Denoising Quality (MP-PCA Noise Map)
# ### QC Checkpoint 1b: Denoising Quality
#
# The MP-PCA denoiser produces a **noise map** (`noise.mif`) that estimates
# voxel-wise noise levels. A good noise map should:
# - Be relatively **uniform** inside the brain (homogeneous noise)
# - Show **higher values outside** the brain (air/background)
# - Have **no structure** that mirrors anatomy (that would mean signal leaked into the noise estimate)

# %% [QC 1b] View Noise Map + Compare Raw vs Denoised DWI
# View noise map — should be spatially smooth, no anatomical structure
launch_gui(f"mrview {DWI}/noise.mif")
print("👆 Noise map should be smooth and uniform inside the brain.")
print("   If you see anatomical structure in the noise, denoising may have failed.\n")

# Compare raw vs denoised DWI (overlay denoised on raw to spot artifacts)
print("Opening raw DWI (before denoise) and denoised DWI for comparison:")
launch_gui(f"mrview {DWI}/dwi_flip_gradcheck.mif")
launch_gui(f"mrview {DWI}/dwi_flip_gradcheck_denoise.mif")
print("👆 Compare: denoised should look smoother but preserve structures.")
print("   No ringing artifacts or unnatural blurring should be present.")

# %% [markdown] QC 1b Result
# **QC 1b Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] QC 1c: Brain Mask Quality
# ### QC Checkpoint 1c: Brain Mask Quality
#
# The brain mask (`mask.mif`) is used for bias correction (step 4g).
# A good mask should:
# - **Cover the entire brain** without gaps
# - **Exclude skull/background**
# - Be reasonably **smooth** at the boundary

# %% [QC 1c] View Brain Mask on DWI (pre-registration)
# View mask overlaid on denoised DWI (pre-registration space)
launch_gui(f"mrview {DWI}/dwi_flip_gradcheck_denoise.mif "
           f"-overlay.load {DWI}/mask.mif -overlay.opacity 0.4")
print("👆 Mask (yellow/red) should cover the whole brain, no holes.")
print("   If mask is poor, bias correction and registration will suffer.\n")

# Also check the FINAL mask in anatomical space
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif "
           f"-overlay.load {DWI}/mask_anat.mif -overlay.opacity 0.4")
print("👆 mask_anat.mif on T2w — should cover brain cleanly in registered space.")

# %% [markdown] QC 1c Result
# **QC 1c Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] QC 2: DWI-T2w Registration
# ### QC Checkpoint 2: DWI-T2w Registration
#
# The ANTs rigid registration (step 4j) aligns the b0 to the T2w.
#
# **Verify:** The b0 image should be perfectly aligned with the T2w:
# - Brain boundaries should match exactly
# - Ventricles, cortical surface, and subcortical structures should overlap
# - Check in **all three planes** (axial, coronal, sagittal)

# %% [QC 2] View DWI-T2w Registration (b0 overlay)
# View b0 overlaid on T2w to check registration
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif "
           f"-overlay.load {DWI}/b0_final.mif -overlay.opacity 0.5")
print("👆 Verify b0 aligns with T2w. Brain edges should match in all 3 planes.\n")

# Also view the ANTs warped result directly (should match T2w grid)
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif "
           f"-overlay.load {DWI}/dwi2anat_Warped.nii.gz -overlay.opacity 0.5")
print("👆 ANTs warped b0 (dwi2anat_Warped.nii.gz) overlaid on T2w.")
print("   This is the direct ANTs output — should match b0_final.mif alignment.")

# %% [markdown] QC 2 Result
# **QC 2 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 2: ROI Preparation
# ## Step 2: ROI Preparation
#
# This step differs by group:
# - **CTRL**: Register Allen CCF atlas to subject T2w, extract ROIs
# - **AP/APTR**: Flip manually-drawn ROIs, correct hemisphere labels
#
# **Run only the section that matches your group.**

# %% [markdown] Step 2A: CTRL — Atlas Coregistration + ROI Extraction
# ### Step 2A: CTRL — Allen Atlas Coregistration + ROI Extraction
#
# Skip this if your subject is AP or APTR.

# %% [Step 2A] CTRL: Atlas Coregistration + ROI Extraction
# ============================================================
# CTRL ONLY: Atlas coregistration and ROI extraction
if GROUP == 'CTRL':
    print('Running atlas ROI extraction for CTRL subject...')
    print('This takes ~5 minutes for ANTs registration.\n')
    result_1 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/atlas_roi_extraction.py --subject {SUBJECT_ID} --session {SESSION_ID}', timeout=600)
    atlas_dir = BASE / 'atlas'
    rois_dir = atlas_dir / 'rois'
    print('\n--- Atlas ROI outputs ---')
    check_file(atlas_dir / 'annotation_warped.nii.gz', 'Warped annotation')
    for roi in ['isocortex', 'striatum', 'thalamus', 'cerebellum', 'olfactory_cortex']:
        check_file(rois_dir / f'roi_{roi}_left.nii.gz', f'{roi} L')
        check_file(rois_dir / f'roi_{roi}_right.nii.gz', f'{roi} R')
else:
    print(f'Skipping Step 2A (GROUP={GROUP}, not CTRL)')

# %% [markdown] Step 2B: AP/APTR — Manual ROI Flip + Hemisphere Fix
# ### Step 2B: AP/APTR — Manual ROI Flip + Hemisphere Correction
#
# Skip this if your subject is CTRL.
#
# **Source:** Original 3D Slicer `.seg.nrrd` files from `paths.manual_roi_slicer_dir`
# (Google Drive shared folder). Each AP/APTR subject has a `slicer_file` field in config1.yaml.
#
# The script (`batch_process_flip1.py --steps 3`) will:
# 1. Read the `.seg.nrrd` and extract each segment as a NIfTI mask
# 2. Swap Left/Right labels (original Slicer labels were inverted)
# 3. Apply flip1 transformation (axes permutation + flip axis 1)
# 4. Write `roi_{Name}_{hemi}_flip.mif` into `manual_roi/`
#
# **Output:** `{subject}/{session}/manual_roi/roi_{Name}_{hemi}_flip.mif`
#
# **Note:** If `slicer_file` is `null` for this subject, ROIs must be drawn
# in 3D Slicer first and the `.seg.nrrd` placed in the shared folder.

# %% [Step 2B] AP/APTR: Manual ROI Flip + Hemisphere Fix
# ============================================================
# AP/APTR ONLY: Manual ROI flip and hemisphere correction
# Uses batch_process_flip1.py step 3 (handles flip + L/R swap)
# ============================================================
if GROUP in ["AP", "APTR"]:
    print('Running ROI flip via batch_process_flip1.py --steps 3...\n')
    result_2b = run(
        f'cd {REPO_DIR} && micromamba run -n mbrain python functions/batch_process_flip1.py '
        f'--subject {SUBJECT_ID} --steps 3',
        timeout=600
    )
    # Verify outputs
    MANUAL_ROI = BASE / 'manual_roi'
    roi_names = ['Graft', 'Striatum', 'Thalamus', 'Palecortex', 'cerebellum']
    print('\n--- Verifying flipped ROIs ---')
    for name in roi_names:
        for hemi in ['left', 'right']:
            check_file(MANUAL_ROI / f'roi_{name}_{hemi}_flip.mif', f'{name} {hemi} (flipped)')
else:
    print(f'Skipping Step 2B (GROUP={GROUP}, not AP/APTR)')

# %% [markdown] QC 3: ROI Alignment on T2w
# ### QC Checkpoint 3: ROI Alignment on T2w
#
# **Verify:** ROIs should match anatomical structures on T2w.

# %% [QC 3] View ROIs on T2w (Alignment Check)
# QC: View ROIs on T2w
if GROUP == 'CTRL':
    rois_dir_1 = BASE / 'atlas' / 'rois'
    cmd = (f'mrview {DWI}/t2w_highres_flip.mif '
           f'-overlay.load {rois_dir_1}/roi_isocortex_left.nii.gz -overlay.opacity 0.4 '
           f'-overlay.load {rois_dir_1}/roi_thalamus_left.nii.gz -overlay.opacity 0.4 '
           f'-overlay.load {rois_dir_1}/roi_striatum_left.nii.gz -overlay.opacity 0.4')
else:
    mroi = BASE / 'manual_roi'
    cmd = (f'mrview {DWI}/t2w_highres_flip.mif '
           f'-overlay.load {mroi}/roi_Graft_left_flip.mif -overlay.opacity 0.4 '
           f'-overlay.load {mroi}/roi_Thalamus_left_flip.mif -overlay.opacity 0.4 '
           f'-overlay.load {mroi}/roi_Striatum_left_flip.mif -overlay.opacity 0.4')
launch_gui(cmd)
print('\n👆 Check that ROIs sit on correct anatomical structures.')
print('   Left-labeled ROIs should appear on anatomical LEFT.')

# %% [markdown] QC 3 Result
# **QC 3 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 3: FOD Estimation
# ## Step 3: FOD Estimation
#
# Generate Fiber Orientation Distributions for tractography.
#
# **Two-stage process** (matching the fixel-based analysis pipeline):
# 1. **Step 3a – Individual response + FOD** (dhollander → msmt_csd → mtnormalise)
#    Produces per-subject response functions and `wm_fod_norm.mif`.
# 2. **Step 3b – Group-response FOD** (uses pre-computed `group_responses/`):
#    Re-estimates FOD with the group-averaged response → `wm_fod_norm_group.mif`.
#    This is the FOD used for ALL downstream analyses (tractography, fixels, etc.).
#
# The group response files already exist in `derivatives/mrtrix_preproc_flip1/group_responses/`
# (computed with `compute_group_response.py` over all subjects).

# %% [Step 3a] Individual Response + FOD (dhollander + msmt_csd + mtnormalise)
print('Step 3a: Individual response estimation + FOD (this may take several minutes)...\n')
result_2 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/dwi_fod.py '
               f'--subject {SUBJECT_ID} --force', timeout=1800)
print('\n--- Individual FOD outputs ---')
for f in ['response_wm.txt', 'response_gm.txt', 'response_csf.txt',
          'wm_fod.mif', 'gm_fod.mif', 'csf_fod.mif',
          'wm_fod_norm.mif', 'csf_fod_norm.mif']:
    check_file(DWI / f)

# %% [Step 3b] Group-Response FOD (msmt_csd with group average + mtnormalise)
# Re-estimate FOD using the existing group-averaged response functions.
# This ensures all subjects' FODs are comparable for fixel-based analysis.
GROUP_RESP = DERIV / 'group_responses'
print(f'Step 3b: Re-estimating FOD with group response...\n')
print(f'  Group responses: {GROUP_RESP}')
for rf in ['group_response_wm.txt', 'group_response_gm.txt', 'group_response_csf.txt']:
    check_file(GROUP_RESP / rf)

# msmt_csd with group response
result_3b_csd = run(
    f'micromamba run -n mbrain dwi2fod msmt_csd '
    f'{DWI}/dwi_flip_gradcheck_denoise_bias_anat.mif '
    f'{GROUP_RESP}/group_response_wm.txt {DWI}/wm_fod_group.mif '
    f'{GROUP_RESP}/group_response_gm.txt {DWI}/gm_fod_group.mif '
    f'{GROUP_RESP}/group_response_csf.txt {DWI}/csf_fod_group.mif '
    f'-mask {DWI}/mask_anat.mif -force -quiet',
    timeout=1800)

# mtnormalise (WM + CSF only)
result_3b_norm = run(
    f'micromamba run -n mbrain mtnormalise '
    f'{DWI}/wm_fod_group.mif {DWI}/wm_fod_norm_group.mif '
    f'{DWI}/csf_fod_group.mif {DWI}/csf_fod_norm_group.mif '
    f'-mask {DWI}/mask_anat.mif -force -quiet',
    timeout=600)

# If mtnormalise fails (common for mouse brain), copy unnormalized
if not (DWI / 'wm_fod_norm_group.mif').exists():
    print('  ⚠ mtnormalise failed — using unnormalized group FOD')
    run(f'cp {DWI}/wm_fod_group.mif {DWI}/wm_fod_norm_group.mif')

print('\n--- Group FOD outputs ---')
for f in ['wm_fod_group.mif', 'gm_fod_group.mif', 'csf_fod_group.mif',
          'wm_fod_norm_group.mif', 'csf_fod_norm_group.mif']:
    check_file(DWI / f)

# %% [markdown] QC 4: FOD Quality
# ### QC Checkpoint 4: FOD Quality
#
# **Verify:** FOD glyphs should show directionally consistent orientations in WM.

# %% [QC 4] View FOD Glyphs on T2w
# View GROUP FOD glyphs on T2w (wm_fod_norm_group.mif — used for all analyses)
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif -odf.load_sh {DWI}/wm_fod_norm_group.mif")
print("\n👆 Check FOD glyphs (GROUP response). They should be:"
      "\n   - Elongated in WM tracts (not spherical)"
      "\n   - Consistent direction along known tracts (corpus callosum, etc.)"
      "\n   - Small/absent in ventricles and outside brain"
      "\n   - NOTE: Using wm_fod_norm_group.mif (group-averaged response)")

# %% [markdown] QC 4 Result
# **QC 4 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 4: Resample ROIs to FOD Space
# ## Step 4: Resample ROIs to FOD Space
#
# ROIs are in T2w resolution but tractography needs them in FOD resolution.
# Uses `mrgrid regrid` with nearest-neighbor interpolation.

# %% [Step 4] Resample ROIs to FOD Space (mrgrid regrid)
# Resample all ROIs to FOD space
print('Resampling ROIs to FOD space...\n')
result_3 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/resample_rois_to_fod.py '
               f'--subject {SUBJECT_ID} --force', timeout=300)
print('\n--- Resampled ROI outputs ---')
if GROUP == 'CTRL':
    fod_rois = BASE / 'atlas' / 'rois' / 'fod_space'
    for roi_1 in ['isocortex', 'striatum', 'thalamus', 'cerebellum', 'olfactory_cortex']:
        check_file(fod_rois / f'roi_{roi_1}_left.nii.gz', f'{roi_1} L (FOD)')
        check_file(fod_rois / f'roi_{roi_1}_right.nii.gz', f'{roi_1} R (FOD)')
else:
    fod_rois = BASE / 'manual_roi' / 'fod_space'
    for roi_1 in ['Graft', 'Striatum', 'Thalamus', 'Palecortex', 'cerebellum']:
        check_file(fod_rois / f'roi_{roi_1}_left_flip_fod.mif', f'{roi_1} L (FOD)')
        check_file(fod_rois / f'roi_{roi_1}_right_flip_fod.mif', f'{roi_1} R (FOD)')

# %% [markdown] QC 5: ROIs in FOD Space
# ### QC Checkpoint 5: ROIs in FOD Space
#
# **Verify:** Resampled ROIs should still align with anatomy on the FOD grid.

# %% [QC 5] View Resampled ROIs on FOD Grid
# View FOD-space ROIs overlaid on FOD
if GROUP == 'CTRL':
    fod_rois_1 = BASE / 'atlas' / 'rois' / 'fod_space'
    cmd_1 = (f'mrview {DWI}/wm_fod_norm_group.mif '
             f'-overlay.load {fod_rois_1}/roi_isocortex_left.nii.gz -overlay.opacity 0.5 '
             f'-overlay.load {fod_rois_1}/roi_thalamus_left.nii.gz -overlay.opacity 0.5')
else:
    fod_rois_1 = BASE / 'manual_roi' / 'fod_space'
    cmd_1 = (f'mrview {DWI}/wm_fod_norm_group.mif '
             f'-overlay.load {fod_rois_1}/roi_Graft_left_flip_fod.mif -overlay.opacity 0.5 '
             f'-overlay.load {fod_rois_1}/roi_Thalamus_left_flip_fod.mif -overlay.opacity 0.5')
launch_gui(cmd_1)
print('\n👆 Check that ROIs align on the FOD grid.')

# %% [markdown] QC 5 Result
# **QC 5 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 5: Segmentation
# ## Step 5: Segmentation
#
# Divide the seed ROI (isocortex/graft) into spatial segments for detailed tractography.
# Segment definitions come from `config1.yaml → passthrough_extraction → segments`.
#
# - **CTRL**: Isocortex → 5 anterior-posterior segments (no D-V split)
# - **AP/APTR**: Graft → SEG_MODE segments (set in Step 0; default 9 = 3 A-P x 3 D-V)
#
# To change segment count for AP/APTR, edit `SEG_MODE` in Step 0.
# To change which segments are used for passthrough, edit `config1.yaml → passthrough_extraction → segments`.

# %% [Step 5] Show Segmentation Parameters
# Show segmentation parameters
print(f"Segmentation config (from config1.yaml):")
print(f"  Group:      {GROUP}")
if GROUP == 'CTRL':
    print(f"  Segments:   {pt.get('segments', {}).get('CTRL', [])}")
    print(f"  Seed ROI:   isocortex (from Allen atlas)")
else:
    print(f"  SEG_MODE:   {SEG_MODE} (set in Step 0 above)")
    print(f"  Segments:   {pt.get('segments', {}).get('AP_APTR', [])}")
    print(f"  Seed ROI:   graft (manually drawn)")
print(f"  To change:  edit SEG_MODE (Step 0) or config1.yaml → passthrough_extraction → segments")

# %% [Step 5] Run Segmentation (Isocortex or Graft)
# ============================================================
# Segmentation: CTRL or AP/APTR
# ============================================================
if GROUP == 'CTRL':
    print('Segmenting isocortex into 5 A-P segments...\n')
    result_4 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/segment_isocortex_slicewise.py '
                   f'--subject {SUBJECT_ID}', timeout=300)
    seg_dir = DWI / 'isocortex_segments_v2'
    fod_dir = seg_dir / 'fod_space'
    positions = ['anterior', 'anterior_central', 'central', 'posterior_central', 'posterior']
    print('\n--- Isocortex segments ---')
    for pos in positions:
        for hemi_1 in ['left', 'right']:
            check_file(seg_dir / f'isocortex_{hemi_1}_{pos}.nii.gz', f'{pos} {hemi_1}')
elif GROUP in ['AP', 'APTR']:
    # AP/APTR: Graft segmentation (9 segments = 3 A-P x 3 D-V)
    print(f'Segmenting graft into {SEG_MODE} segments...\n')
    result_4 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/segment_graft_slicewise.py '
                   f'--subject {SUBJECT_ID} --session {SESSION_ID} --segments {SEG_MODE}', timeout=300)
    seg_t2w = DWI / 'graft_segments_t2w'
    seg_fod = DWI / 'graft_segments_fod'
    if SEG_MODE == 9:
        ap_labels = ['anterior', 'central', 'posterior']
        dv_labels = ['dorsal', 'medial', 'ventral']
    else:
        ap_labels = ['anterior', 'central', 'posterior']
        dv_labels = ['dorsal', 'ventral']
    print(f'\n--- Graft segments (seg{SEG_MODE}) ---')
    for ap in ap_labels:
        for dv in dv_labels:
            for hemi_1 in ['left', 'right']:
                check_file(seg_t2w / f'seg{SEG_MODE}_graft_{hemi_1}_{ap}_{dv}.nii.gz',
                           f'seg{SEG_MODE} {hemi_1} {ap}_{dv}')

# %% [markdown] QC 6: Segmentation
# ### QC Checkpoint 6: Segmentation
#
# **Verify:** Segments should tile the seed ROI without gaps or overlaps.

# %% [QC 6a] Check Segment Overlaps (Validation)
# Check for overlaps using the validation utility
result_5 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python utils/check_roi_overlap.py '
               f'--subject {SUBJECT_ID} --session {SESSION_ID} --group {GROUP}', timeout=120)
print('\n👆 Verify: total unique voxels == sum of individual segments (no overlap).')

# %% [QC 6b] View Segments on T2w (Visual Check)
# Visual QC: overlay all segments on T2w
if GROUP == 'CTRL':
    seg_dir_1 = DWI / 'isocortex_segments_v2'
    cmd_2 = f'mrview {DWI}/t2w_highres_flip.mif'
    for pos_1 in ['anterior', 'central', 'posterior']:
        cmd_2 += f' -overlay.load {seg_dir_1}/isocortex_left_{pos_1}.nii.gz -overlay.opacity 0.4'
else:
    seg_dir_1 = DWI / 'graft_segments_t2w'
    cmd_2 = f'mrview {DWI}/t2w_highres_flip.mif'
    for seg in ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']:
        cmd_2 += f' -overlay.load {seg_dir_1}/seg{SEG_MODE}_graft_left_{seg}.nii.gz -overlay.opacity 0.4'
launch_gui(cmd_2)
print('\n👆 Segments should tile the seed ROI from anterior to posterior.')

# %% [markdown] QC 6 Result
# **QC 6 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 6: Whole-Brain Tractography
# ## Step 6: Whole-Brain Tractography
#
# Generate a dense whole-brain tractogram seeded from the brain mask.
# This is the basis for passthrough extraction and SIFT2 filtering.
#
# Parameters from `config1.yaml → wholebrain_tractography` + `tractography`.
# To change streamline count, edit `config1.yaml → wholebrain_tractography → num_streamlines`.
# To change step/angle/cutoff, edit `config1.yaml → tractography`.
#
# **This step takes 30-90 minutes.**

# %% [Step 6] Show Whole-Brain Tractography Parameters
# Show whole-brain tractography parameters
print("Whole-brain tractography parameters (from config1.yaml):")
print(f"  num_streamlines: {wb.get('num_streamlines', 'N/A'):,}")
print(f"  algorithm:       {tck['algorithm']}")
print(f"  step_size:       {tck['step_size']} mm")
print(f"  angle:           {tck['angle']}°")
print(f"  cutoff:          {tck['cutoff']}")
print(f"  min_length:      {tck['min_length']} mm")
print(f"  max_length:      {tck['max_length']} mm")
print(f"  FOD used:        wm_fod_norm_group.mif")
print(f"  Mask:            mask_anat.mif")
print(f"  To change:       edit config1.yaml → wholebrain_tractography / tractography")

# %% [Step 6] Run Whole-Brain Tractography (1M streamlines)
# Generate whole-brain streamlines via generate_wholebrain_tractogram.py
# This script reads all params from config1.yaml and uses wm_fod_norm_group.mif
WB_DIR = DWI / 'wholebrain'
WB_TCK = WB_DIR / 'wholebrain_1000k.tck'

print('Running whole-brain tractography via generate_wholebrain_tractogram.py...')
print('This takes 30-90 minutes.\n')
result_6 = run(
    f'cd {REPO_DIR} && micromamba run -n mbrain python functions/generate_wholebrain_tractogram.py '
    f'--subject {SUBJECT_ID}',
    timeout=7200
)
check_file(WB_TCK, 'Whole-brain tractogram')

# %% [markdown] QC 7: Whole-Brain Tractogram
# ### QC Checkpoint 7: Whole-Brain Tractogram
#
# **Verify:** Tractogram should fill the brain volume with reasonable fiber directions.

# %% [QC 7] View Whole-Brain Tractogram in mrview
# Quick QC: tract info + visualization
run(f"tckinfo {WB_TCK}")
print("\nOpening mrview...")
launch_gui(f"mrview {DWI}/t2w_highres_flip.mif -tractography.load {WB_TCK}")
print("\n👆 Fibers should fill the brain. Check for obvious artifacts.")

# %% [markdown] QC 7 Result
# **QC 7 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 7: SIFT2 Weights
# ## Step 7: SIFT2 Weights for Whole-Brain Tractogram
#
# Compute SIFT2 cross-section multipliers so each streamline gets a weight
# proportional to how well the tractogram matches the underlying FOD.

# %% [Step 7] Compute SIFT2 Weights
# SIFT2 via run_passthrough_pipeline.py --step sift2
# This script uses wm_fod_norm_group.mif and produces sift2_weights.txt
SIFT2_WEIGHTS = WB_DIR / 'sift2_weights.txt'
SIFT2_MU = WB_DIR / 'sift2_mu.txt'

print('Running SIFT2 via run_passthrough_pipeline.py...\n')
result_7 = run(
    f'cd {REPO_DIR} && micromamba run -n mbrain python functions/run_passthrough_pipeline.py '
    f'--subject {SUBJECT_ID} --step sift2',
    timeout=1800
)
check_file(SIFT2_WEIGHTS, 'SIFT2 weights')
if SIFT2_WEIGHTS.exists():
    run(f'wc -l {SIFT2_WEIGHTS}')

# %% [markdown] Step 8: Passthrough Tract Extraction
# ---
# ## Step 8: Passthrough Tract Extraction
#
# Extract streamlines that pass through specific ROI pairs from the whole-brain tractogram.
# Uses `tckedit` with `-include` flags.
#
# Parameters from `config1.yaml → passthrough_extraction`.
# To change targets/segments, edit that section.
#
# For AP/APTR: also extracts using **shrunk** (eroded) graft segments
# (controlled by `config1.yaml → shrink_analysis`).

# %% [Step 8] Show Passthrough Extraction Parameters
# Show passthrough extraction parameters
print("Passthrough extraction parameters (from config1.yaml):")
if GROUP == 'CTRL':
    print(f"  Targets:    {pt.get('targets', {}).get('CTRL', [])}")
    print(f"  Segments:   {pt.get('segments', {}).get('CTRL', [])}")
else:
    print(f"  Targets:    {pt.get('targets', {}).get('AP_APTR', [])}")
    print(f"  Segments:   {pt.get('segments', {}).get('AP_APTR', [])}")
print(f"  Hemispheres:  {pt.get('hemispheres', [])}")
print(f"  Cerebellum includes thalamus waypoint: {pt.get('cerebellum_include_thalamus', False)}")
if GROUP in ['AP', 'APTR']:
    print(f"  Shrink analysis enabled: {sk.get('enabled', False)}")
    print(f"  Erosion voxels: {sk.get('erosion_voxels', 'N/A')}")
print(f"  To change: edit config1.yaml → passthrough_extraction / shrink_analysis")

# %% [Step 8] Extract Passthrough Tracts — LEFT Hemisphere
# Extract passthrough tracts — LEFT hemisphere
print('Extracting passthrough tracts (LEFT hemisphere)...\n')
result_8 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/extract_passthrough_tracts.py '
               f'--subject {SUBJECT_ID} --hemi left', timeout=1800)
pt_dir = WB_DIR / 'passthrough_tracts'
run(f'ls -la {pt_dir}/ | head -20', check=False)

# %% [Step 8] Extract Passthrough Tracts — RIGHT Hemisphere
# Extract passthrough tracts — RIGHT hemisphere
print('Extracting passthrough tracts (RIGHT hemisphere)...\n')
result_9 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/extract_passthrough_tracts.py '
               f'--subject {SUBJECT_ID} --hemi right', timeout=1800)
run(f'ls -la {pt_dir}/ | head -20', check=False)

# %% [Step 8] AP/APTR: Extract with Shrunk (Eroded) ROIs
# ============================================================
# AP/APTR ONLY: Also extract with shrunk ROIs
# ============================================================
if GROUP in ['AP', 'APTR']:
    # First create shrunk segments if they don't exist
    print('Creating shrunk (eroded) graft segments...\n')
    run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/extract_passthrough_tracts.py '
        f'--subject {SUBJECT_ID} --create-shrink-segments', timeout=300)
    # Extract with shrunk ROIs — both hemispheres
    for hemi_2 in ['left', 'right']:
        print(f'\nExtracting shrunk passthrough tracts ({hemi_2.upper()})...')
        run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/extract_passthrough_tracts.py '
            f'--subject {SUBJECT_ID} --hemi {hemi_2} --shrink', timeout=1800)
    shrink_dir = WB_DIR / 'passthrough_tracts_shrink'
    run(f'ls -la {shrink_dir}/ | head -20', check=False)
else:
    print("Skipping shrunk extraction (CTRL subjects don't use shrunk ROIs)")

# %% [markdown] QC 8: Passthrough Tracts
# ### QC Checkpoint 8: Passthrough Tracts
#
# **Verify:** Check a few tracts to make sure they connect the right regions.

# %% [QC 8a] List All Extracted Tracts + Streamline Counts
# List all extracted tracts with streamline counts
import glob

pt_dir_1 = WB_DIR / 'passthrough_tracts'
if pt_dir_1.exists():
    tracts = sorted(glob.glob(str(pt_dir_1 / '*.tck')))
    print(f'Found {len(tracts)} passthrough tract files:\n')
    for t in tracts:
        r = subprocess.run(f'tckinfo {t} -count', shell=True, capture_output=True, text=True)
        count = r.stdout.strip().split()[-1] if r.returncode == 0 else '?'
        print(f'  {Path(t).name}: {count} streamlines')
else:
    print('No passthrough tracts found yet.')

# %% [QC 8b] View Sample Tract (Thalamus Left) in mrview
# Visual QC: view a sample tract (thalamus left)
if GROUP == 'CTRL':
    sample_tract = pt_dir_1 / 'isocortex_left_anterior_to_thalamus_left.tck'
    roi_overlay = BASE / 'atlas' / 'rois' / 'fod_space' / 'roi_thalamus_left.nii.gz'
else:
    # For AP/APTR use first thalamus tract found
    sample_tract = list(sorted(pt_dir_1.glob('*thalamus_left*.tck')))
    sample_tract = sample_tract[0] if sample_tract else None
    roi_overlay = BASE / 'manual_roi' / 'fod_space' / 'roi_Thalamus_left_flip_fod.mif'

if sample_tract and Path(sample_tract).exists():
    cmd_5 = (f'mrview {DWI}/t2w_highres_flip.mif '
             f'-overlay.load {roi_overlay} -overlay.opacity 0.4 '
             f'-tractography.load {sample_tract}')
    launch_gui(cmd_5)
    print(f'\n👆 Showing: {Path(sample_tract).name}')
    print('   Fibers should connect seed region to thalamus.')
else:
    print('No sample tract found for visualization.')

# %% [markdown] QC 8 Result
# **QC 8 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 9: MATLAB AFQ Cleaning
# ## Step 9: MATLAB AFQ Cleaning
#
# Remove outlier streamlines using MATLAB `AFQ_removeFiberOutliers`.
# All parameters are read from `config1.yaml → matlab_cleaning`.
# To change thresholds, edit that section and re-run Step 0.

# %% [Step 9] Show MATLAB AFQ Cleaning Parameters
# Show MATLAB AFQ cleaning parameters (from config1.yaml → matlab_cleaning)
print("MATLAB AFQ Cleaning Parameters (from config1.yaml → matlab_cleaning):")
print(f"  max_dist:        {mat['max_dist']}  ← Gaussian distance from core (SD). Higher=more permissive")
print(f"  max_len:         {mat['max_len']}  ← Length deviation (SD). Higher=more permissive")
print(f"  num_nodes:       {mat['num_nodes']}  ← Resampling nodes per fiber")
print(f"  min_streamlines: {mat['min_streamlines']}  ← Skip cleaning if tract has fewer")
print(f"  max_iter:        {mat.get('max_iter', 'N/A')}  ← Max cleaning iterations")
print(f"  core_method:     {mat.get('core_method', 'N/A')}  ← How tract core is computed")
print(f"  matlab_bin:      {mat.get('matlab_bin', 'N/A')}")
print(f"  afq_path:        {mat.get('afq_path', 'N/A')}")
print(f"\n  To change: edit config1.yaml → matlab_cleaning")

# %% [Step 9] Run MATLAB AFQ Cleaning (Passthrough Tracts)
# Clean passthrough tracts
print('Cleaning passthrough tracts with MATLAB AFQ...\n')
print('(This runs MATLAB in batch mode — may take 5-20 minutes)\n')
result_10 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/clean_tracts.py '
                f'--subject {SUBJECT_ID} --passthrough', timeout=3600)
n_cleaned = len(glob.glob(str(pt_dir_1 / '*_cleaned.tck')))
print(f'\nCleaned files created: {n_cleaned}')

# %% [Step 9] AP/APTR: Clean Shrunk Passthrough Tracts
# ============================================================
# AP/APTR ONLY: Also clean shrunk tracts
# ============================================================
if GROUP in ['AP', 'APTR']:
    print('Cleaning shrunk passthrough tracts...\n')
    result_11 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/clean_tracts.py '
                    f'--subject {SUBJECT_ID} --passthrough --shrink', timeout=3600)
    shrink_dir_2 = WB_DIR / 'passthrough_tracts_shrink'
    n_cleaned_shrink = len(glob.glob(str(shrink_dir_2 / '*_cleaned.tck')))
    print(f'\nCleaned shrunk files: {n_cleaned_shrink}')
else:
    print('Skipping shrunk cleaning (CTRL)')

# %% [markdown] QC 9: Cleaned vs Raw Tracts
# ### QC Checkpoint 9: Cleaned vs Raw Tracts
#
# **Verify:** Cleaned tracts should retain the main bundle but remove outliers.

# %% [QC 9] Compare Raw vs Cleaned Tract (Before/After)
# Compare raw vs cleaned for one tract
cleaned_files = sorted(glob.glob(str(pt_dir_1 / '*_cleaned.tck')))
if cleaned_files:
    cleaned = Path(cleaned_files[0])
    raw = cleaned.parent / cleaned.name.replace('_cleaned.tck', '.tck')
    print(f'Comparing: {raw.name}')
    run(f'tckinfo {raw} -count')
    print(f'Cleaned: {cleaned.name}')
    run(f'tckinfo {cleaned} -count')
    launch_gui(f'mrview {DWI}/t2w_highres_flip.mif -tractography.load {cleaned}')
    print('\n👆 Cleaned tract should be a tight, clean bundle.')
else:
    print('No cleaned tracts found yet.')

# %% [markdown] QC 9 Result
# **QC 9 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 10: Generate Metrics Report
# ## Step 10: Generate Metrics Report
#
# Compute tract metrics and generate CSV + HTML report:
# - Streamline count (from cleaned tracts)
# - SIFT2 connectivity (from raw tracts matched to whole-brain weights)
# - Mean AFD (Apparent Fiber Density) along tracts
#
# **Note:** AP/APTR automatically uses shrunk tracts for the report.

# %% [Step 10] Generate Tract Metrics Report (CSV + HTML)
# Generate metrics report
import datetime
version = f"v{datetime.datetime.now().strftime('%Y%m%d')}_{SUBJECT_ID}"
print(f'Generating metrics report (version: {version})...\n')
result_12 = run(f'cd {REPO_DIR} && micromamba run -n mbrain python functions/generate_tract_metrics_report.py '
                f'--cleaned --subject {SUBJECT_ID} --version {version}', timeout=600)
reports_dir = REPO_DIR / 'reports'
check_file(reports_dir / f'tract_metrics_{version}.csv', 'CSV report')
check_file(reports_dir / f'tract_metrics_report_{version}.html', 'HTML report')

# %% [Step 10] Preview CSV Report
# Preview the CSV report
import pandas as pd

csv_path = REPO_DIR / "reports" / f"tract_metrics_{version}.csv"
if csv_path.exists():
    df = pd.read_csv(csv_path)
    print(f"Report: {len(df)} rows\n")
    print("Columns:", list(df.columns))
    print(f"\nSummary for {SUBJECT_ID}:")
    display(df[df['subject'] == SUBJECT_ID].to_string(index=False))
else:
    print("Report CSV not found.")

# %% [markdown] Step 11: Fixel-Based Analysis (FBA)
# ---
# ## Step 11: Fixel-Based Analysis (FBA)
#
# Decompose FODs into discrete fixels and compute microstructural metrics.
#
# **FBA pipeline:**
# 1. **Generate fixels** — `fod2fixel` on `wm_fod_norm_group.mif` → `fixels_group/`
#    (produces `afd.mif`, `peak.mif`, `directions.mif`, `index.mif`)
# 2. **Voxel AFD map** — `fixel2voxel` collapses fixel AFD to voxel space → `afd_voxel_sum.mif`
# 3. **ROI fixel metrics** — directional FD within each ROI (AP/DV/LR decomposition)
# 4. **Tract AFD sampling** — `tcksample` samples AFD along cleaned passthrough tracts
#
# **Key metrics:**
# - `fd_total` — Total Fiber Density in ROI
# - `fd_ap`, `fd_dv`, `fd_lm` — Directional FD (Anterior-Posterior, Dorsal-Ventral, Lateral-Medial)
# - `fd_ap_frac`, `fd_dv_frac`, `fd_lm_frac` — Directional fractions
# - `mean_afd` — Mean AFD along a tract path
#
# **Axis convention (after flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 11a] Generate Fixels (fod2fixel)
# Generate fixel decomposition from group-calibrated FOD
FIXEL_DIR = DWI / 'fixels_group'
print(f'Step 11a: Generating fixels from group FOD...\n')
print(f'  Input:  {DWI}/wm_fod_norm_group.mif')
print(f'  Mask:   {DWI}/mask_anat.mif')
print(f'  Output: {FIXEL_DIR}/')

fod_file = DWI / 'wm_fod_norm_group.mif'
mask_file = DWI / 'mask_anat.mif'

if not fod_file.exists():
    print(f'\n  ❌ ERROR: Group FOD not found! Run Step 3b first.')
elif FIXEL_DIR.exists() and (FIXEL_DIR / 'afd.mif').exists():
    print(f'\n  ✅ Fixels already exist. Use --force in manual run to regenerate.')
    for ff in ['index.mif', 'directions.mif', 'afd.mif', 'peak.mif']:
        check_file(FIXEL_DIR / ff)
else:
    result_11a = run(
        f'micromamba run -n mbrain fod2fixel '
        f'{DWI}/wm_fod_norm_group.mif {FIXEL_DIR} '
        f'-mask {DWI}/mask_anat.mif '
        f'-afd afd.mif -peak peak.mif -force -quiet',
        timeout=600)
    print('\n--- Fixel outputs ---')
    for ff in ['index.mif', 'directions.mif', 'afd.mif', 'peak.mif']:
        check_file(FIXEL_DIR / ff)

# %% [Step 11b] Create Voxel AFD Map (fixel2voxel)
# Collapse fixel AFD values to a voxel-level sum for tract sampling
AFD_VOXEL = DWI / 'afd_voxel_sum.mif'
print(f'Step 11b: Creating voxel AFD map...\n')

if AFD_VOXEL.exists():
    print(f'  ✅ afd_voxel_sum.mif already exists.')
    check_file(AFD_VOXEL)
elif not (FIXEL_DIR / 'afd.mif').exists():
    print(f'  ❌ ERROR: Fixels not found! Run Step 11a first.')
else:
    result_11b = run(
        f'micromamba run -n mbrain fixel2voxel '
        f'{FIXEL_DIR}/afd.mif sum {AFD_VOXEL} -force',
        timeout=120)
    check_file(AFD_VOXEL)

# %% [Step 11c] ROI-Based Fixel Metrics (Directional FD)
# Compute directional fiber density within each ROI using the fixel decomposition.
# This runs analyze_roi_fixels.py on each ROI, producing CSV with FD_total, FD_ap, FD_dv, FD_lm.
print(f'Step 11c: Computing ROI-based fixel metrics...\n')

if GROUP == 'CTRL':
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

fixel_csv = REPO_DIR / 'reports' / f'fixel_metrics_{SUBJECT_ID}.csv'
for roi_name in roi_names:
    print(f'  Processing ROI: {roi_name}')
    result_11c = run(
        f'cd {REPO_DIR} && micromamba run -n mbrain python functions/analyze_roi_fixels.py '
        f'--subject {SUBJECT_ID} --roi {roi_name} '
        f'--output {fixel_csv} --force',
        timeout=300)

if fixel_csv.exists():
    print(f'\n--- ROI Fixel Metrics ---')
    import pandas as pd
    df_fixel = pd.read_csv(fixel_csv)
    print(f'  {len(df_fixel)} ROI entries')
    display(df_fixel.to_string(index=False))
else:
    print(f'\n  ⚠ Fixel CSV not generated (check errors above)')

# %% [Step 11d] Tract AFD Sampling (tcksample along cleaned tracts)
# Sample the voxel AFD map along each cleaned passthrough tract.
# This gives mean_afd per tract — a microstructural metric independent of streamline count.
import numpy as np
print(f'Step 11d: Sampling AFD along cleaned tracts...\n')

AFD_VOXEL = DWI / 'afd_voxel_sum.mif'
if not AFD_VOXEL.exists():
    print(f'  ❌ ERROR: afd_voxel_sum.mif not found! Run Step 11b first.')
else:
    # Find all cleaned tract files
    pt_dir_afd = WB_DIR / 'passthrough_tracts'
    cleaned_tracts_afd = sorted(glob.glob(str(pt_dir_afd / '*_cleaned.tck'))) if pt_dir_afd.exists() else []
    
    # Also check shrunk tracts for AP/APTR
    if GROUP in ['AP', 'APTR']:
        shrink_dir_afd = WB_DIR / 'passthrough_tracts_shrink'
        cleaned_tracts_afd += sorted(glob.glob(str(shrink_dir_afd / '*_cleaned.tck'))) if shrink_dir_afd.exists() else []
    
    afd_results = []
    for tract_path in cleaned_tracts_afd:
        tract_name = Path(tract_path).stem.replace('_cleaned', '')
        afd_out = Path(tract_path).parent / f'{tract_name}_afd.txt'
        
        result_afd = run(
            f'micromamba run -n mbrain tcksample '
            f'{tract_path} {AFD_VOXEL} {afd_out} '
            f'-stat_tck mean -force',
            timeout=120)
        
        if afd_out.exists():
            # Read mean AFD values
            with open(afd_out) as f:
                lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
            if lines:
                values = [float(v) for v in lines if v.replace('.','').replace('-','').replace('e','').isdigit() or '.' in v]
                if values:
                    mean_afd = np.mean(values)
                    afd_results.append({'tract': tract_name, 'mean_afd': mean_afd, 'n_streamlines': len(values)})
                    print(f'  ✅ {tract_name}: mean_afd={mean_afd:.6f} (n={len(values)})')
    
    if afd_results:
        import pandas as pd
        df_afd = pd.DataFrame(afd_results)
        print(f'\n--- Tract AFD Summary ---')
        display(df_afd.to_string(index=False))
    else:
        print('  ⚠ No AFD results (are cleaned tracts available?)')

# %% [Step 11e] Slicewise Fixel Analysis (Directional FD by coronal slice)
# Compute fixel metrics per coronal slice within the seed ROI.
# This produces the "radar plot" data used in the group comparison reports.
print(f'Step 11e: Slicewise fixel analysis...\n')

# Run the 3-group comparison analysis (includes this subject's data)
# This generates radar plots comparing CTRL vs AP vs APTR per slice
for hemi in ['left', 'right']:
    print(f'\n  --- {hemi.upper()} hemisphere ---')
    result_11e = run(
        f'cd {REPO_DIR} && micromamba run -n mbrain python functions/analyze_roi_fixels.py '
        f'--compare --force',
        timeout=600)

print('\n--- Slicewise fixel reports ---')
for hemi in ['left', 'right']:
    for shrink_tag in ['', '_shrink']:
        fig_path = REPO_DIR / 'reports' / f'fixel_3group_comparison_{hemi}{shrink_tag}.png'
        if fig_path.exists():
            check_file(fig_path)

# %% [Step 11e-km] 3-Group Fixel K-means Direction Comparison
# Same layout as Step 11e, but uses K-means clustering + scatter matrix
# eigenvector method to identify the 2 main fiber directions per slice.
# Solid = Dir 1 (largest cluster), Dashed = Dir 2.
import sys, yaml, importlib
import numpy as np
from pathlib import Path
sys.path.insert(0, str(REPO_DIR))
import functions.fixel_roi_analysis as _fri
importlib.reload(_fri)
from functions.fixel_roi_analysis import analyze_3group_kmeans

print('Step 11e-km: 3-Group Fixel K-means Direction Comparison\n')

# Load config
with open(REPO_DIR / 'config1.yaml') as f:
    cfg = yaml.safe_load(f)

report_dir = REPO_DIR / 'reports'
combos = [(h, True) for h in ['left', 'right']]  # always use shrink ROIs

# --- Pass 1: compute data only (no plots) to find unified y_max ---
print('Pass 1: computing data for all combinations...')
all_results = {}
global_max = 0.0
for hemi, use_shrink in combos:
    tag = f'{hemi}_shrink'
    r = analyze_3group_kmeans(
        cfg, hemisphere=hemi, use_shrink=True,
        output_dir=report_dir, n_clusters=2, conda_env='mbrain',
        force=False, show=False, verbose=False, plot=False)
    all_results[tag] = (hemi, True, r)
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

# --- Pass 2: generate plots with forced y_max ---
for tag, (hemi, use_shrink, r) in all_results.items():
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

# %% [markdown] QC 11: Fixel Analysis
# ### QC Checkpoint 11: Fixel Analysis
#
# **Verify:**
# - Fixels generated without errors (Step 11a)
# - AFD voxel map looks reasonable (higher in WM, low in CSF)
# - ROI metrics show expected directional patterns:
#   - Graft/isocortex: mixed directions (cortical fanning)
#   - Thalamus: strong AP component (internal capsule)
#   - Cerebellum: characteristic DV pattern
# - Tract AFD values consistent with other subjects in same group

# %% [QC 11] View AFD Voxel Map on T2w
# Check afd_voxel_sum.mif across all 12 real subjects (excluding test mice)
print('QC 11: Checking AFD voxel maps for all subjects\n')

skip_ids = {'sub-APtestmouseNmM'}
found = 0
missing = 0
for group_name in ['CTRL', 'AP', 'APTR']:
    print(f'\n--- {group_name} ---')
    for entry in cfg['groups'][group_name]:
        sid, ses = entry['id'], entry['session']
        if sid in skip_ids:
            continue
        afd_path = DERIV / sid / ses / 'dwi' / 'afd_voxel_sum.mif'
        if afd_path.exists():
            print(f'  ✅ {sid}: afd_voxel_sum.mif exists')
            found += 1
        else:
            print(f'  ❌ {sid}: afd_voxel_sum.mif NOT FOUND — run Steps 11a+11b for this subject')
            missing += 1

print(f'\n--- Summary: {found} found, {missing} missing ---')
if found > 0:
    # Open the first available one for visual check
    for group_name in ['CTRL', 'AP', 'APTR']:
        for entry in cfg['groups'][group_name]:
            sid, ses = entry['id'], entry['session']
            if sid in skip_ids:
                continue
            afd_path = DERIV / sid / ses / 'dwi' / 'afd_voxel_sum.mif'
            t2w_path = DERIV / sid / ses / 'dwi' / 't2w_highres_flip.mif'
            if afd_path.exists():
                print(f'\nOpening {sid} for visual check...')
                launch_gui(f"mrview {t2w_path} -overlay.load {afd_path} -overlay.opacity 0.5")
                print('👆 AFD voxel map. Should be bright in WM regions, dark in CSF/outside brain.')
                break
        else:
            continue
        break

# %% [markdown] QC 11 Result
# **QC 11 Result:** ☐ PASS / ☐ FAIL
#
# ---

# %% [markdown] Step 11f: Fixel Square ROI Analysis
# ---
# ## Step 11f: Fixel Square ROI Direction Analysis
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
# - `wm_fod_norm_group.mif` (from Step 3b)
#
# **Antipodal symmetry:** Fiber directions are sign-ambiguous (v ≡ −v).
# Before clustering, all vectors are flipped to the same hemisphere (z ≥ 0).
#
# **Axis convention (flip1):** axis 0 = AP, axis 1 = DV, axis 2 = LR

# %% [Step 11f] Fixel Square ROI — Generate fixels & extract directions
from functions.fixel_roi_analysis import generate_fixel_peaks

SQUARE_ROI_T2W = DWI / 'square_roi_fixel.mif'
FOD_PATH = DWI / 'wm_fod_norm_group.mif'

print('Step 11f: Fixel Square ROI Direction Analysis\n')

if not SQUARE_ROI_T2W.exists():
    print('  ❌ ERROR: square_roi_fixel.mif not found!')
    print('  Draw a square ROI in mrview (Tool → ROI editor) and save as:')
    print(f'  {SQUARE_ROI_T2W}')
elif not FOD_PATH.exists():
    print('  ❌ ERROR: wm_fod_norm_group.mif not found! Run Step 3b first.')
else:
    peaks_nii = generate_fixel_peaks(
        roi_path=SQUARE_ROI_T2W,
        fod_path=FOD_PATH,
        output_dir=DWI,
        conda_env='mbrain',
        force=True,
    )
    print('\n  ✅ Fixel extraction complete.')

# %% [Step 11f-2] Fixel Square ROI — K-means clustering of directions
import numpy as np
from functions.fixel_roi_analysis import extract_peak_vectors, cluster_fixel_directions

PEAKS_NII = DWI / 'square_roi_peaks.nii'
print('Step 11f-2: Clustering fixel directions within square ROI\n')

if not PEAKS_NII.exists():
    print('  ❌ Peaks file not found — run the cell above first.')
else:
    all_vectors, all_peak_indices, n_peaks = extract_peak_vectors(PEAKS_NII)
    if len(all_vectors) == 0:
        print('  ⚠ No valid fixel vectors found in the square ROI.')
    else:
        print(f'  Total fixel vectors: {len(all_vectors)}')
        print(f'    Peak 1: {np.sum(all_peak_indices == 1)}')
        if n_peaks > 1:
            print(f'    Peak 2: {np.sum(all_peak_indices == 2)}')

        cluster_result = cluster_fixel_directions(
            all_vectors, all_peak_indices, n_peaks,
            n_clusters=2, verbose=True,
        )

# %% [Step 11f-3] Fixel Square ROI — 3D Direction Visualization
from functions.fixel_roi_analysis import plot_fixel_directions

if 'cluster_result' in dir() and cluster_result is not None:
    fig_path = REPO_DIR / 'reports' / f'fixel_square_roi_{SUBJECT_ID}.png'
    fig = plot_fixel_directions(
        cluster_result,
        subject_id=SUBJECT_ID,
        save_path=fig_path,
        show=True,
    )
else:
    print('  ❌ Run Step 11f-2 first.')

# %% [markdown] Step 12: Final Visualization & Summary
# ## Step 12: Final Visualization & Summary
#
# View the complete results for this subject.

# %% [Step 12] Pipeline Summary Table
# Summary table: all tracts for this subject
print(f"\n{'=' * 80}")
print(f' PIPELINE COMPLETE: {SUBJECT_ID} ({GROUP})')
print(f"{'=' * 80}\n")

pt_dir_2 = WB_DIR / 'passthrough_tracts'
raw_tracts = glob.glob(str(pt_dir_2 / '*.tck')) if pt_dir_2.exists() else []
raw_tracts = [t for t in raw_tracts if '_cleaned' not in t]
cleaned_tracts = glob.glob(str(pt_dir_2 / '*_cleaned.tck')) if pt_dir_2.exists() else []

print(f"Whole-brain tractogram:  {'✅' if WB_TCK.exists() else '❌'}")
print(f"SIFT2 weights:          {'✅' if SIFT2_WEIGHTS.exists() else '❌'}")
print(f'Raw passthrough tracts:  {len(raw_tracts)}')
print(f'Cleaned tracts:          {len(cleaned_tracts)}')

if GROUP in ['AP', 'APTR']:
    shrink_dir_3 = WB_DIR / 'passthrough_tracts_shrink'
    shrink_raw = glob.glob(str(shrink_dir_3 / '*.tck')) if shrink_dir_3.exists() else []
    shrink_raw = [t for t in shrink_raw if '_cleaned' not in t]
    shrink_cleaned = glob.glob(str(shrink_dir_3 / '*_cleaned.tck')) if shrink_dir_3.exists() else []
    print(f'Shrunk raw tracts:       {len(shrink_raw)}')
    print(f'Shrunk cleaned tracts:   {len(shrink_cleaned)}')

print(f'\nReport: reports/tract_metrics_{version}.csv')
print(f'HTML:   reports/tract_metrics_report_{version}.html')
print(f"\n{'=' * 80}")

# %% [Step 12] Final Visualization — All Cleaned Tracts in mrview
# Final comprehensive visualization
# Open mrview with T2w + best tracts
if GROUP == 'CTRL':
    seg_roi = BASE / 'atlas' / 'rois' / 'fod_space' / 'roi_isocortex_left.nii.gz'
    target_roi = BASE / 'atlas' / 'rois' / 'fod_space' / 'roi_thalamus_left.nii.gz'
else:
    seg_roi = BASE / 'manual_roi' / 'fod_space' / 'roi_Graft_left_flip_fod.mif'
    target_roi = BASE / 'manual_roi' / 'fod_space' / 'roi_Thalamus_left_flip_fod.mif'

cleaned_1 = sorted(glob.glob(str(pt_dir_2 / '*left*_cleaned.tck')))
cmd_6 = (f'mrview {DWI}/t2w_highres_flip.mif '
         f'-overlay.load {seg_roi} -overlay.opacity 0.3 '
         f'-overlay.load {target_roi} -overlay.opacity 0.3')
# Load up to 6 tracts
for t_1 in cleaned_1[:6]:
    cmd_6 += f' -tractography.load {t_1}'
launch_gui(cmd_6)
print('\n👆 Final visualization with all left-hemisphere cleaned tracts.')

# %% [markdown] Bar-plots: APTR SIFT2 Density
# ## Bar-plots: APTR SIFT2 Density (Left Hemisphere & Central Graft)
#
# Bar plots showing SIFT2 density (connectivity / seed_voxels × 1000) for all
# APTR mice. Individual data points overlaid on group means.

# %% [Bar-plots] APTR Density — Left hemisphere & Central summary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── Load metrics CSV (cleaned tracts) ────────────────────────────────────
csv_path = Path("reports/tract_metrics_v05_cleaned.csv")
df = pd.read_csv(csv_path)

# Filter APTR, left hemisphere only
aptr_left = df[(df['group'] == 'APTR') & (df['hemisphere'] == 'left')].copy()
targets_order = ['olfactory_cortex', 'striatum', 'thalamus', 'cerebellum']
target_labels = {'thalamus': 'Thalamus', 'striatum': 'Striatum',
                 'olfactory_cortex': 'Paleocortex', 'cerebellum': 'Cerebellum'}
segments_order = ['anterior', 'central', 'posterior']
seg_colors = {'anterior': '#E24A33', 'central': '#348ABD', 'posterior': '#988ED5'}
panel_color = '#555555'

subjects = aptr_left['subject'].unique()
n_subj = len(subjects)
marker_list = ['o', 's', 'D', '^']  # one per mouse
n_tgt = len(targets_order)

# ── Helper: draw one bar subplot ─────────────────────────────────────────
def _bar_panel(ax, data_df, seg_list, title, show_legend=False):
    """Draw grouped bars (one bar per segment) for 4 targets, with individual points."""
    n_seg = len(seg_list)
    bar_width = 0.22
    group_width = n_seg * bar_width + 0.15
    x_tgt = np.arange(n_tgt) * group_width

    for s_idx, seg in enumerate(seg_list):
        seg_data = data_df[data_df['segment'] == seg]
        means, sems = [], []
        for tgt in targets_order:
            vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
            means.append(np.mean(vals) if len(vals) > 0 else 0)
            sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)

        x_pos = x_tgt + s_idx * bar_width
        ax.bar(x_pos, means, bar_width * 0.85, yerr=sems, capsize=4,
               color=seg_colors[seg], alpha=0.55, edgecolor='black', linewidth=0.5,
               label=seg.capitalize(), zorder=2)

        for tgt_idx, tgt in enumerate(targets_order):
            tgt_rows = seg_data[seg_data['target'] == tgt]
            for _, row in tgt_rows.iterrows():
                subj_idx = list(subjects).index(row['subject'])
                jitter = (np.random.RandomState(subj_idx).rand() - 0.5) * bar_width * 0.5
                ax.scatter(x_pos[tgt_idx] + jitter, row['sift2_density'],
                           marker=marker_list[subj_idx % len(marker_list)],
                           color=seg_colors[seg], edgecolor='black', s=40, zorder=3,
                           linewidth=0.6)

    ax.set_xticks(x_tgt + bar_width * (n_seg - 1) / 2)
    ax.set_xticklabels([target_labels[t] for t in targets_order], fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    if show_legend:
        ax.legend(title='Segment', fontsize=9)


def _single_seg_panel(ax, data_df, seg, color, title):
    """Draw single-segment bars for 4 targets, with individual points + subject legend."""
    seg_data = data_df[data_df['segment'] == seg]
    means, sems, all_vals = [], [], []
    for tgt in targets_order:
        vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
        means.append(np.mean(vals) if len(vals) > 0 else 0)
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        all_vals.append(vals)

    x_c = np.arange(n_tgt)
    ax.bar(x_c, means, 0.55, yerr=sems, capsize=4,
           color=color, alpha=0.55, edgecolor='black', linewidth=0.5, zorder=2)

    for tgt_idx in range(n_tgt):
        tgt_subjs = seg_data[seg_data['target'] == targets_order[tgt_idx]]['subject'].values
        for val_idx, val in enumerate(all_vals[tgt_idx]):
            subj_idx = list(subjects).index(tgt_subjs[val_idx])
            short = tgt_subjs[val_idx].replace('sub-', '')
            jitter = (np.random.RandomState(subj_idx).rand() - 0.5) * 0.3
            label = short if tgt_idx == 0 else None
            ax.scatter(x_c[tgt_idx] + jitter, val,
                       marker=marker_list[subj_idx % len(marker_list)],
                       color=color, edgecolor='black', s=50, zorder=3,
                       linewidth=0.6, label=label)

    ax.set_xticks(x_c)
    ax.set_xticklabels([target_labels[t] for t in targets_order], fontsize=9,
                       rotation=30, ha='right')
    ax.set_title(title, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(title='Subject', fontsize=6, loc='upper right')


# ── Build figure: 3 segment bar plots + 2 mrview screenshots ─────────────
# Layout: top row = Anterior | Central | Posterior
#         bottom row = mrview least connected | mrview most connected
import matplotlib.gridspec as gridspec
from matplotlib.image import imread as mpl_imread

fig = plt.figure(figsize=(20, 14))
fig.suptitle('APTR — SIFT2 Density (Left Hemisphere, Cleaned Tracts)',
             fontsize=14, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3,
                       height_ratios=[1, 0.9])

# Row 0: 3 segment bar subplots
ax_ant = fig.add_subplot(gs[0, 0])
_single_seg_panel(ax_ant, aptr_left, 'anterior', seg_colors['anterior'],
                  'Anterior Graft → Destinations')
ax_ant.set_ylabel('SIFT2 Density\n(conn / seed_vox × 1000)', fontsize=10)

ax_cen = fig.add_subplot(gs[0, 1])
_single_seg_panel(ax_cen, aptr_left, 'central', seg_colors['central'],
                  'Central Graft → Destinations')
ax_cen.set_ylabel('SIFT2 Density', fontsize=10)

ax_pos = fig.add_subplot(gs[0, 2])
_single_seg_panel(ax_pos, aptr_left, 'posterior', seg_colors['posterior'],
                  'Posterior Graft → Destinations')
ax_pos.set_ylabel('SIFT2 Density', fontsize=10)

# Row 1: mrview screenshots by destination for sub-Apallial20R1F
# Save screenshots from mrview as:
#   reports/mrview_most_paleocortex.png  reports/mrview_most_striatum.png  reports/mrview_most_thalamus.png
destinations = ['paleocortex', 'striatum', 'thalamus']
dest_labels  = {'paleocortex': 'Paleocortex', 'striatum': 'Striatum', 'thalamus': 'Thalamus'}

for col_idx, dest in enumerate(destinations):
    ax = fig.add_subplot(gs[1, col_idx])
    img_path = Path(f'reports/mrview_most_{dest}.png')
    title = f'{dest_labels[dest]}\nsub-Apallial20R1F'
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

save_path = Path("reports/barplot_aptr_density_left.png")
fig.savefig(save_path, dpi=150, bbox_inches='tight')
print(f"Saved: {save_path}")
plt.show()

# %% [markdown] QC Summary & Final Screenshot
# **Final Screenshot:**
#
# _Paste screenshot here_
#
# ---
#
# ## QC Summary
#
# | QC  | Check | Result |
# |-----|-------|--------|
# | 1   | T2w orientation matches Allen atlas (dorsal up, anterior left in sagittal) | ☐ PASS / ☐ FAIL |
# | 1b  | Noise map uniform, no anatomical structure; denoised DWI smooth without ringing | ☐ PASS / ☐ FAIL |
# | 1c  | Brain mask covers entire brain, no holes; mask_anat clean on T2w | ☐ PASS / ☐ FAIL |
# | 2   | DWI-T2w registration — b0 brain edges match T2w in all 3 planes | ☐ PASS / ☐ FAIL |
# | 3   | ROIs on correct anatomy (L on L, R on R) | ☐ PASS / ☐ FAIL |
# | 4   | FOD quality (directional glyphs in WM, small in ventricles) | ☐ PASS / ☐ FAIL |
# | 5   | ROIs in FOD space (aligned on FOD grid) | ☐ PASS / ☐ FAIL |
# | 6   | Segmentation (no overlaps, tiles seed ROI) | ☐ PASS / ☐ FAIL |
# | 7   | Whole-brain tractogram (fills brain volume) | ☐ PASS / ☐ FAIL |
# | 8   | Passthrough tracts (connect correct regions) | ☐ PASS / ☐ FAIL |
# | 9   | Cleaned tracts (tight bundles, outliers removed) | ☐ PASS / ☐ FAIL |
# | 10  | Tract metrics report generated (CSV + HTML) | ☐ PASS / ☐ FAIL |
# | 11  | AFD voxel map reasonable (bright WM, dark CSF); ROI fixel metrics consistent | ☐ PASS / ☐ FAIL |
#
# **Processed by:** _Your name_
# **Date:** _Date_
# **Notes:** _Any issues or observations_
