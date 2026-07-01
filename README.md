# mbrain - Mouse Brain fMRI and DTI Analysis

Analysis pipeline for mouse brain imaging data using MRtrix3 and PyBIDS.

## New-Mask Pipeline: Detailed Guide and Running Schema

This repository currently contains a large operational notebook/script at
`notebooks/pipeline_groups_new_mask.py` that performs mask regeneration,
downstream reruns, CST processing, tract QC, and reporting.

To make execution safer and easier to resume, the workflow is now documented
as staged notebooks with an explicit QC gate at the end of each stage.

### Preconditions (must pass before Stage 1)

1. You are running on lmc02 with required modules available:
    - freesurfer/7.4.1-Rocky8
    - mrtrix3/3.0.3
    - ants/latest
    - matlab/R2024A
2. `config1.yaml` paths are correct for the current machine.
3. `DERIV_SRC` and destination `DERIV_BASE` are accessible and writable.
4. Subject roster in config is final for the run.

### Stage-by-stage running schema

#### Stage 1: Setup and Atlas Mask Generation

- Notebook: `notebooks/newmask_split/01_setup_and_masks.py`
- Scope:
   - Environment setup, config load, subject manifest
   - Copy derivatives tree
   - Group atlas preparation
   - Atlas-to-subject mask generation
- Expected outputs:
   - subject manifest in memory (`subjects`)
   - `mask_new.mif` and `mask_new.nii.gz` per subject
- End-of-stage QC gate:
   - all subjects have `mask_new.mif`
   - no failed registrations remain unresolved

#### Stage 2: Mask QC, Manual Edits, and DWI-Space Install

- Notebook: `notebooks/newmask_split/02_mask_qc_and_install.py`
- Scope:
   - visual QC (old vs new)
   - manual correction loop
   - mask regrid to DWI and install
   - checkpoint readiness report
- Expected outputs:
   - `mask_new_edited.mif` where needed
   - installed `mask.mif` and `mask_anat.mif` on DWI grid
- End-of-stage QC gate:
   - DWI-grid match for masks
   - checkpoint indicates Stage 3 can proceed

#### Stage 3: FOD and ROI Preparation

- Notebook: `notebooks/newmask_split/03_fod_and_roi_prep.py`
- Scope:
   - FOD rerun (5b)
   - group response and group-normalized FODs (5b-group)
   - ROI resampling + segmentation (5c)
   - ROI space/laterality QC
- Expected outputs:
   - `wm_fod_norm_group.mif` per subject
   - updated FOD-space ROIs
- End-of-stage QC gate:
   - all required group-FOD outputs exist
   - ROI QC issues reviewed/accepted

#### Stage 4: Whole-Brain Tractography and SIFT2

- Notebook: `notebooks/newmask_split/04_tractography_and_sift2.py`
- Scope:
   - whole-brain tractography rerun
   - 10k downsample and TRK conversion
   - SIFT2 recomputation
- Expected outputs:
   - wholebrain tractograms + downsample + trk
   - sift2 weights
- End-of-stage QC gate:
   - expected tractogram files are present
   - SIFT2 weights are present and up-to-date

#### Stage 5: CST Pipeline (AP/APTR)

- Notebook: `notebooks/newmask_split/05_cst_pipeline.py`
- Scope:
   - CST ROI status/create/edit/regrid
   - CST passthrough extraction
   - targeted CST tckgen
   - CST cleaning
- Expected outputs:
   - FOD-space CST ROIs per AP/APTR subject
   - CST raw/tckgen/cleaned tracts
- End-of-stage QC gate:
   - CST ROIs complete in FOD space
   - CST tract counts look consistent

#### Stage 6: Passthrough and Cleaning (Global + CTRL rerun)

- Notebook: `notebooks/newmask_split/06_passthrough_and_cleaning.py`
- Scope:
   - passthrough extraction (5f)
   - AFQ cleaning (5g)
   - CTRL-only rerun branch (5fg)
- Expected outputs:
   - raw and cleaned passthrough tracts
   - refreshed CTRL tract sets if rerun used
- End-of-stage QC gate:
   - per-subject raw/cleaned counts are non-empty where expected
   - spot checks raw vs cleaned streamline counts

#### Stage 7: Reports, Fixel Analyses, and Final QC

- Notebook: `notebooks/newmask_split/07_reports_and_fixel.py`
- Scope:
   - tract metrics report generation (5h)
   - fixel workflows (5i, 5j)
   - final visual QC and summary tables
- Expected outputs:
   - CSV/HTML metrics reports
   - fixel figures and per-subject outputs
   - final summary tables
- End-of-stage QC gate:
   - metrics and fixel artifacts present
   - final QC table completed

### Handoff contract between stages

At the top of each stage notebook:

1. Validate required input artifacts from previous stage.
2. Abort early if any critical dependency is missing.
3. Write a stage checkpoint JSON to:
    - `DERIV_BASE/logs/checkpoints/NN_stage_name.json`

Suggested checkpoint JSON fields:

- `stage`
- `timestamp`
- `operator`
- `status` (pass/fail)
- `required_inputs`
- `produced_outputs`
- `qc_summary`
- `notes`

### Current split files

The split scaffolds are located at:

- `notebooks/newmask_split/01_setup_and_masks.py`
- `notebooks/newmask_split/02_mask_qc_and_install.py`
- `notebooks/newmask_split/03_fod_and_roi_prep.py`
- `notebooks/newmask_split/04_tractography_and_sift2.py`
- `notebooks/newmask_split/05_cst_pipeline.py`
- `notebooks/newmask_split/06_passthrough_and_cleaning.py`
- `notebooks/newmask_split/07_reports_and_fixel.py`

## Project Structure

```
mbrain/
├── config1.yaml             # Central configuration (ALL parameters)
├── config_loader.py         # Python config access utility
├── main.py                  # Main analysis script
│
├── # Generic Pipeline Scripts (parameter-driven, no hardcoding)
├── atlas_coregistration.py  # Atlas registration (all subjects/groups)
├── generate_tracts.py       # Tractography generation (all tracts)
├── clean_tracts.py          # MATLAB AFQ fiber cleaning
├── extract_top_tracts.py    # Extract top N by SIFT2 weight
├── split_rois.py            # Split ROIs into hemispheres
├── transform_tracts.py      # Transform tracts to templates
├── preprocess_subjects.py   # DWI preprocessing
├── view_tracts.py           # Visualization with mrview
├── build_group_templates.py # Build group average templates
│
├── utils/                   # Utility modules
│   ├── __init__.py
│   ├── config.py           # Legacy config (use config1.yaml)
│   ├── file_utils.py       # File handling utilities
│   ├── bids_converter.py   # Flywheel to BIDS converter
│   └── dwi_preprocessor.py # MRtrix3 preprocessing pipeline
├── local/                   # Local data (not in git)
│   └── atlas/allen_ccf/    # Allen CCF atlas files
└── data/                    # Processed data
    ├── bids/               # BIDS-formatted dataset
    └── derivatives/        # Analysis outputs
        └── mrtrix_preproc/ # Preprocessed DWI data
```

## Configuration

All parameters are centralized in `config1.yaml`:

```yaml
# Paths, groups, tractography params, cleaning params, ROIs, tracts
# Edit this file to change ANY parameter - no code changes needed
```

Key sections:
- **paths**: All directory paths
- **groups**: Subject lists (CTRL, AP, APTR) with sessions
- **tractography**: tckgen parameters (algorithm, select, step, angle, etc.)
- **matlab_cleaning**: AFQ_removeFiberOutliers parameters
- **rois**: Allen CCF structure IDs for each ROI
- **tracts**: Tract pair definitions (seed → target)

## Generic Pipeline Usage

All scripts use the same parameter pattern:
- `--group CTRL/AP/APTR` - Process specific group
- `--subject sub-xxx` - Process specific subject
- `--tract Isocortex_Striatum` - Process specific tract
- `--hemisphere left/right/both` - Hemisphere selection
- `--dry-run` - Preview commands without executing

### Examples

```bash
# Atlas coregistration
python atlas_coregistration.py                    # All subjects
python atlas_coregistration.py --group CTRL       # CTRL only
python atlas_coregistration.py --subject sub-129-921-R1F

# Generate tracts
python generate_tracts.py                         # All tracts, all subjects
python generate_tracts.py --tract Isocortex_Striatum --hemisphere left
python generate_tracts.py --group AP --dry-run

# Clean tracts (MATLAB AFQ)
python clean_tracts.py --subject sub-129-921-R1F
python clean_tracts.py --group CTRL --tract Isocortex_Thalamus

# Extract top N tracts
python extract_top_tracts.py --top-n 100 --cleaned
python extract_top_tracts.py --group CTRL --tract Isocortex_Striatum

# View tracts
python view_tracts.py --subject sub-129-921-R1F
python view_tracts.py --group CTRL --template
python view_tracts.py --list-tracts
```

## Setup

### 1. Create Environment

```bash
micromamba create -n mbrain python=3.12 numpy pandas scipy matplotlib seaborn \
  scikit-learn nibabel nilearn ipywidgets ipykernel jupyter notebook pyarrow \
  pybids pyyaml -c conda-forge

micromamba activate mbrain
```

### 2. Install ANTs (Required for bias field correction)

```bash
# macOS with Homebrew
brew install ants

# Or download from: http://stnava.github.io/ANTs/
```

### 3. Verify MRtrix3 Installation

```bash
mrinfo --version
```


## Usage

### Command-Line Options

```bash
# Run all steps (BIDS conversion + preprocessing)
python main.py --all
# or simply
python main.py

# Only run BIDS conversion
python main.py --bids

# Only run preprocessing
python main.py --preprocess

# Preview BIDS conversion (dry-run)
python main.py --bids --dry-run

# Run specific preprocessing steps
python main.py --preprocess --steps denoise unring
python main.py --preprocess --steps denoise unring mask bias_correct

# Get help
python main.py --help
```

### Pipeline Steps

**Step 1: BIDS Conversion**
- Converts Flywheel data structure to BIDS format
- Organizes subjects, sessions, and acquisitions
- Handles anatomical (T2w, T2*) and diffusion (DTI) data
- Copies JSON sidecars with acquisition metadata
- Creates BIDS metadata files

**Step 2: DWI Preprocessing (MRtrix3)**
1. **Denoising** - Marchenko-Pastur PCA (`dwidenoise`)
2. **Unringing** - Gibbs artifact removal (`mrdegibbs`)
3. **Brain Mask** - Mask estimation (`dwi2mask` or threshold-based)
4. **Bias Correction** - B1 field correction (`dwibiascorrect ants`)

## Data Organization

### Input (Flywheel Structure)
```
local/flywheel/fMRI/SUBJECTS/
  ├── {subject_id}/
  │   └── SESSIONS/
  │       └── {session_name}/
  │           └── ACQUISITIONS/
  │               ├── T2w_highres/
  │               ├── T2w_lowres/
  │               ├── DTI_EPI_seg_30dir_sat/
  │               └── ...
```

### Output (BIDS Structure)
```
data/bids/
  ├── dataset_description.json
  ├── participants.tsv
  └── sub-{subject_id}/
      └── ses-{session_id}/
          ├── anat/
          │   ├── sub-{id}_ses-{ses}_acq-highres_T2w.nii.gz
          │   ├── sub-{id}_ses-{ses}_acq-highres_T2w.json
          │   └── sub-{id}_ses-{ses}_acq-lowres_T2w.nii.gz
          └── dwi/
              ├── sub-{id}_ses-{ses}_dwi.nii.gz
              ├── sub-{id}_ses-{ses}_dwi.json
              ├── sub-{id}_ses-{ses}_dwi.bval
              └── sub-{id}_ses-{ses}_dwi.bvec
```

### Preprocessed Data Structure
```
data/derivatives/mrtrix_preproc/
  └── sub-{subject_id}/
      └── ses-{session_id}/
          └── dwi/
              ├── sub-{id}_ses-{ses}_denoised.mif
              ├── sub-{id}_ses-{ses}_noise.mif
              ├── sub-{id}_ses-{ses}_unringed.mif
              ├── sub-{id}_ses-{ses}_mask.mif
              ├── sub-{id}_ses-{ses}_preprocessed.mif
              └── sub-{id}_ses-{ses}_preprocessed.nii.gz
```

## Troubleshooting

### Common Errors and Solutions

#### 1. ANTs Not Found (`dwibiascorrect` fails)

**Error:**
```
dwibiascorrect: [ERROR] Unable to find ANTS program "N4BiasFieldCorrection"
```

**Solution:**
```bash
# Install ANTs
brew install ants

# Verify installation
which N4BiasFieldCorrection

# If not in PATH, add to ~/.zshrc:
export PATH="/opt/homebrew/bin:$PATH"
source ~/.zshrc
```

#### 2. Brain Mask Fails for Mouse Data

**Error:**
```
dwi2mask: [WARNING] ... failed (expected for ex-vivo/mouse data)
```

**Solution:**
This is expected for mouse/ex-vivo data. The pipeline automatically creates a simple threshold-based mask as fallback. However, you should **manually check and correct masks** if needed:

```bash
# View the mask
mrview data/derivatives/mrtrix_preproc/sub-XXX/ses-YYY/dwi/sub-XXX_ses-YYY_mask.mif

# If mask is poor, create manually:
mrthreshold data/bids/sub-XXX/ses-YYY/dwi/sub-XXX_ses-YYY_dwi.nii.gz \
  custom_mask.mif -abs 20

# Or use FSL's bet (may need tuning for mouse brains)
bet data/bids/sub-XXX/ses-YYY/dwi/sub-XXX_ses-YYY_dwi.nii.gz \
  brain_mask -m -f 0.3
```

#### 3. Missing Gradient Information

**Error:**
```
dwidenoise: [ERROR] no diffusion gradient table found
```

**Solution:**
Ensure `.bval` and `.bvec` files are present alongside the DWI NIfTI:
```bash
# Check files exist
ls data/bids/sub-XXX/ses-YYY/dwi/

# Should see:
# sub-XXX_ses-YYY_dwi.nii.gz
# sub-XXX_ses-YYY_dwi.bval
# sub-XXX_ses-YYY_dwi.bvec
```

If missing, check the BIDS conversion step or original Flywheel data.

#### 4. Preprocessing Fails Mid-Pipeline

**Solution:**
The pipeline saves intermediate outputs. You can continue from where it failed:

```bash
# Run only specific steps that failed
python main.py --preprocess --steps mask bias_correct

# Or manually run MRtrix commands on intermediate files
dwibiascorrect ants \
  data/derivatives/mrtrix_preproc/sub-XXX/ses-YYY/dwi/sub-XXX_ses-YYY_unringed.mif \
  output.mif \
  -mask mask.mif
```

#### 5. Out of Memory Errors

**Solution:**
Mouse brain data can be high-resolution. Reduce memory usage:

```bash
# Limit number of threads
export MRTRIX_NTHREADS=2

# Or in commands
dwidenoise input.mif output.mif -nthreads 2
```

#### 6. File Permission Errors

**Solution:**
```bash
# Ensure write permissions
chmod -R u+w data/derivatives/mrtrix_preproc/

# Or run with force flag (already included in pipeline)
dwidenoise input.mif output.mif -force
```

### Quality Control Checklist

After preprocessing, always check:

1. **Residuals from denoising** - Should show noise, not anatomy:
   ```bash
   mrcalc dwi.mif denoised.mif -subtract residuals.mif
   mrview residuals.mif
   ```

2. **Brain mask coverage** - Should cover entire brain:
   ```bash
   mrview preprocessed.mif -roi.load mask.mif
   ```

3. **Bias field correction** - Compare before/after:
   ```bash
   mrview unringed.mif preprocessed.mif -mode 2
   ```

## Next Steps

After successful preprocessing:

1. **Validate preprocessed data**
   ```bash
   mrview data/derivatives/mrtrix_preproc/sub-XXX/ses-YYY/dwi/sub-XXX_ses-YYY_preprocessed.mif
   ```

2. **Estimate response functions**
   ```bash
   dwi2response tournier preprocessed.mif response.txt
   ```

3. **Compute fiber orientation distributions**
   ```bash
   dwi2fod csd preprocessed.mif response.txt fod.mif -mask mask.mif
   ```

4. **Perform tractography**
   ```bash
   tckgen fod.mif tracks.tck -seed_image mask.mif -select 10000
   ```

## Configuration

Edit `utils/config.py` to customize:
- Data paths
- BIDS dataset metadata
- Modality mappings
- Acquisition naming conventions

## Getting Help

- MRtrix3 Documentation: https://mrtrix.readthedocs.io/
- MRtrix3 Community Forum: https://community.mrtrix.org/
- PyBIDS Documentation: https://bids-standard.github.io/pybids/
