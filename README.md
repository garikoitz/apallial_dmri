# apallial_dmri - Apallial Mouse Brain dMRI Analysis

Analysis pipeline for apallial mouse brain imaging data using MRtrix3.

## Two Ways to Run This Pipeline

You can use this repository in one of two ways, depending on your goal:

- **Option A — Generate the paper figures from installed data (fast).**
  Install the pre-computed derivatives with `python utils/setup_data.py`, then
  run only the figure-generation cells in
  `notebooks/generate_figures_only.py` (or `.ipynb`). This reads the
  already-processed data and produces the images used in the paper — one
  figure per cell, no preprocessing required.

- **Option B — Run the full pipeline from BIDS data (complete reproduction).**
  Install the BIDS + local components, then run the preprocessing and
  downstream stages in `notebooks/pipeline_groups_new_mask.py` starting from
  raw BIDS data. This reproduces every intermediate output (masks, FODs,
  tractograms, cleaning, reports) end to end. This option also requires
  external tools not bundled in the conda env: `antspyx`
  (`pip install antspyx`).

See [SETUP.md](SETUP.md) for the component-by-component install commands and
the "Basic Replication Workflow" section below for which notebook steps to run
for each option.

## Setup

For a complete step-by-step guide, see [SETUP.md](SETUP.md).

### 1. Clone and Create Environment

If you have micromamba, this is the code.  Equivalently, you can use conda instead of micromamba (or any other environment manager). This step takes a few minutes.

```bash
git clone https://github.com/garikoitz/apallial_dmri.git
cd apallial_dmri
micromamba env create -f mbrain_env.yml
micromamba activate mbrain
```

### 2. Download Test Data

```bash
python utils/setup_data.py
```

This downloads and extracts all 4 OSF components to `local/` and `data/`.
See [SETUP.md](SETUP.md) for component-by-component download options.

### 3. Verify External Dependencies

Full pipeline runs require MRtrix3 3.0.4 or newer.

On macOS, use the official MRtrix3 standalone installer:

```bash
sudo bash -c "$(curl -fsSL https://raw.githubusercontent.com/MRtrix3/macos-installer/master/install)"
```

Then verify that MRtrix3 is available:

```bash
mrinfo --version
```

We recommend this installer on macOS because Homebrew and Conda may not provide
a usable current MRtrix3 package, especially on Apple Silicon.

```bash
# ANTs (required for bias field correction)
# macOS: brew install ants
# Or download from: http://stnava.github.io/ANTs/
```


## Pipeline Overview

The complete pipeline is in a single notebook:

- **`notebooks/pipeline_groups_new_mask.py`** — The authoritative, complete, runnable pipeline for all subjects/groups.

The notebook imports reusable helper functions from `functions/pipeline_new_mask/` (11 modules with unit tests). All parameters come from `config1.yaml`.

### Pipeline Stages

The notebook runs through these stages:

1. **Setup & Atlas Mask Generation** — Environment setup, config load, subject manifest, group atlas preparation, atlas-to-subject mask generation
2. **Mask QC, Manual Edits & DWI-Space Install** — Visual QC, manual correction loop, mask regrid to DWI and install, checkpoint readiness
3. **FOD & ROI Preparation** — FOD rerun, group response and group-normalized FODs, ROI resampling + segmentation, ROI space/laterality QC
4. **Whole-Brain Tractography & SIFT2** — Whole-brain tractography rerun, 10k downsample and TRK conversion, SIFT2 recomputation
5. **CST Pipeline (AP/APTR)** — CST ROI status/create/edit/regrid, CST passthrough extraction, targeted CST tckgen, CST cleaning
6. **Passthrough & Cleaning (Global + CTRL rerun)** — Passthrough extraction, AFQ cleaning, CTRL-only rerun branch
7. **Reports, Fixel Analyses & Final QC** — Tract metrics report generation, fixel workflows, final visual QC and summary tables

### Preconditions

1. Required modules available: MRtrix3 (for both options); MATLAB and FreeSurfer only for Option B (full pipeline)
2. `config1.yaml` paths are correct for the current machine
3. Subject roster in config is final for the run
4. Data downloaded via `python utils/setup_data.py`

## Project Structure

```
apallial_dmri/
├── config1.yaml             # Central configuration (ALL parameters)
├── mbrain_env.yml           # Micromamba environment definition
├── README.md                # This file
├── SETUP.md                 # Quick setup guide (new computer onboarding)
│
├── functions/               # Pipeline scripts (parameter-driven, no hardcoding)
│   ├── __init__.py
│   ├── config_loader.py     # Python config access utility
│   ├── atlas_coregistration.py  # Atlas registration (all subjects/groups)
│   ├── atlas_roi_extraction.py  # ROI extraction from warped atlas
│   ├── batch_process_flip1.py   # Batch DWI preprocessing
│   ├── build_group_templates.py # Build group average templates
│   ├── clean_tracts.py          # MATLAB AFQ fiber cleaning
│   ├── compute_group_response.py # Group RF computation
│   ├── dwi_fod.py               # DWI preprocessing + FOD estimation
│   ├── generate_tracts.py       # Tractography generation
│   ├── generate_wholebrain_tractogram.py # Whole-brain tractography
│   ├── resample_rois_to_fod.py  # ROI resampling to FOD space
│   ├── segment_graft_slicewise.py   # Graft segmentation (AP/APTR)
│   ├── segment_isocortex_slicewise.py # Isocortex segmentation (CTRL)
│   ├── view_tracts.py           # Visualization with mrview
│   ├── qc_viewer.py             # QC visualization
│   ├── generate_qc_report.py    # QC report generation
│   ├── generate_connectivity_report.py # Connectivity reporting
│   ├── generate_tract_metrics_report.py # Tract metrics reporting
│   ├── fixel_roi_analysis.py    # Fixel-based ROI analysis
│   ├── fixel_cluster_analysis.py # Fixel cluster analysis
│   ├── fod_direction_analysis.py # FOD direction analysis
│   ├── collect_bilateral_data.py # Bilateral data collection
│   ├── extract_passthrough_tracts.py # Passthrough tract extraction
│   ├── run_passthrough_pipeline.py   # Full passthrough pipeline
│   ├── matlab_clean_standalone.py    # Standalone MATLAB cleaning
│   ├── validate_config.py       # Config validation
│   └── _tests_/                 # Test suite
│
├── notebooks/               # Pipeline notebook + documentation
│   ├── pipeline_groups_new_mask.py  # Main pipeline notebook (all subjects/groups)
│   ├── generate_figures_only.py     # Figures-only notebook (Option A, one figure per cell)
│   ├── generate_figures_only.ipynb  # Same as .py, in native Jupyter format
│   └── methods.md                   # Methods documentation
│
├── utils/                   # Utility modules
│   ├── setup_data.py        # OSF data download & extraction
│   └── ...
│
├── docs/                    # Documentation
│   ├── installation.md
│   ├── manual_preprocessing_commands.md
│   ├── manual_tractography_commands.md
│   ├── OSF_DATA_PACKAGING.md
│   ├── OSF_WORKFLOW_CHECKLIST.md
│   ├── pipeline_all_subjects.md
│   ├── two_stage_coregistration.md
│   └── visualize_graft_tracts.md
│
├── for_checking/            # Archived/legacy scripts (not active)
│
├── local/                   # Downloaded data: atlas + manual ROIs (gitignored)
│   └── atlas/allen_ccf/    # Allen CCF atlas files
│
└── data/                    # Downloaded data: BIDS + derivatives (gitignored)
    ├── bids/               # BIDS-formatted dataset
    └── derivatives/        # Analysis outputs
        └── mrtrix_newmask/ # New-mask pipeline derivatives
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

## Script Usage Reference

All scripts are in `functions/` and use the same parameter pattern:

- `--group CTRL/AP/APTR` - Process specific group
- `--subject sub-xxx` - Process specific subject
- `--tract Isocortex_Striatum` - Process specific tract
- `--hemisphere left/right/both` - Hemisphere selection
- `--dry-run` - Preview commands without executing

### Key Scripts

```bash
# Atlas coregistration
python functions/atlas_coregistration.py                    # All subjects
python functions/atlas_coregistration.py --group CTRL       # CTRL only
python functions/atlas_coregistration.py --subject sub-129-921-R1F

# Generate tracts
python functions/generate_tracts.py                         # All tracts, all subjects
python functions/generate_tracts.py --tract Isocortex_Striatum --hemisphere left
python functions/generate_tracts.py --group AP --dry-run

# Clean tracts (MATLAB AFQ)
python functions/clean_tracts.py --subject sub-129-921-R1F
python functions/clean_tracts.py --group CTRL --tract Isocortex_Thalamus

# View tracts
python functions/view_tracts.py --subject sub-129-921-R1F
python functions/view_tracts.py --group CTRL --template
python functions/view_tracts.py --list-tracts

# Whole-brain tractography
python functions/generate_wholebrain_tractogram.py --subject sub-129-921-R1F

# ROI resampling to FOD space
python functions/resample_rois_to_fod.py --group CTRL

# Segmentation
python functions/segment_isocortex_slicewise.py --group CTRL
python functions/segment_graft_slicewise.py --group AP

# Fixel analysis
python functions/fixel_roi_analysis.py --group CTRL
python functions/fixel_cluster_analysis.py --group CTRL

# Reports
python functions/generate_tract_metrics_report.py --group CTRL
python functions/generate_connectivity_report.py --group CTRL
```

## Usage

### Main Pipeline Notebook

The primary workflow is in `notebooks/pipeline_groups_new_mask.py`. It is a
**jupytext percent-format** Python script — a plain `.py` file that VS Code
renders as an interactive notebook (cells, outputs, kernel) thanks to the
`jupytext` package (already in `mbrain_env.yml`).

#### Opening it as an interactive notebook in VS Code (recommended)

1. Install the **Jupyter** and **Python** extensions in VS Code
   (Extensions marketplace → search "Jupyter" / "Python").
2. Select the `mbrain` micromamba environment as the Python interpreter:
   `Cmd+Shift+P` → `Python: Select Interpreter` → pick `mbrain`.
3. Open `notebooks/pipeline_groups_new_mask.py` in VS Code.
4. VS Code detects the jupytext percent format and shows an
   **"Open in Notebook"** button at the top-right of the editor (or run the
   command `Jupyter: Open as Notebook`). Click it.
5. The file opens as an interactive notebook: each `# %%` block becomes a cell,
   you can run cells with `Shift+Enter`, see inline outputs, and the kernel is
   the `mbrain` env. Edits save back to the same `.py` file (no `.ipynb` needed).

> No `jupyter notebook` server is required — everything runs inside VS Code
> using the `mbrain` kernel.

#### Alternative: JupyterLab / classic server

If you prefer a browser-based server instead of VS Code:

```bash
micromamba run -n mbrain jupyter lab notebooks/pipeline_groups_new_mask.py
```

`jupytext` (in `mbrain_env.yml`) lets JupyterLab render the `.py` file as a
notebook. The earlier `jupyter notebook ...` form also works, but `jupyter lab`
is recommended (notebook 7 / JupyterLab 4.6+).

This notebook runs the complete new-mask pipeline (see stage-by-stage guide above).

### Basic Replication Workflow (Recommended)

If your goal is to replicate the project results with minimal ambiguity, use this sequence.

1. Edit `config1.yaml` first

   - `paths.repo_dir`: set to your local repository path
   - `paths.data_base`: set to your data location
   - `paths.matlab_path`: verify MATLAB binary path for your machine
   - `groups`: confirm the subject/session roster you want to run
2. Install data with `setup_data.py` (preferred over manual OSF download)

```bash
python utils/setup_data.py --check-only
python utils/setup_data.py --components local pre_bids pre_preproc post
```

Component names are exactly: `local`, `pre_bids`, `pre_preproc`, `post`.

3. Open `notebooks/pipeline_groups_new_mask.py` in VS Code as an interactive notebook

   - Use **Open in Notebook** and select the `mbrain` kernel.
4. Run notebook steps according to your starting point

   - If `data/derivatives/mrtrix_newmask/` is already populated:
     - Skip Steps 1-4
     - Start at Step 0, then continue through the downstream pipeline
   - If starting from BIDS + local data only:
     - Run Step 1 (verify/download checks)
     - Run Step 2 (batch preprocessing)
     - Run Step 3 (FOD + ROI resampling)
     - Run Step 4 (segmentation + tractography)
     - Then run Step 0 and continue
5. Validate key outputs before long downstream reruns

   - Confirm subject folders exist under `data/derivatives/mrtrix_newmask/`
   - Confirm masks and FODs are present for all subjects
   - Do visual QC before launching heavy tractography/cleaning branches

### Individual Scripts

All pipeline scripts are in `functions/` and use `config1.yaml` for parameters:

```bash
# Atlas coregistration
python functions/atlas_coregistration.py                    # All subjects
python functions/atlas_coregistration.py --group CTRL       # CTRL only
python functions/atlas_coregistration.py --subject sub-129-921-R1F

# Generate tracts
python functions/generate_tracts.py                         # All tracts, all subjects
python functions/generate_tracts.py --tract Isocortex_Striatum --hemisphere left
python functions/generate_tracts.py --group AP --dry-run

# Clean tracts (MATLAB AFQ)
python functions/clean_tracts.py --subject sub-129-921-R1F
python functions/clean_tracts.py --group CTRL --tract Isocortex_Thalamus

# View tracts
python functions/view_tracts.py --subject sub-129-921-R1F
python functions/view_tracts.py --group CTRL --template
python functions/view_tracts.py --list-tracts
```

## Configuration

Edit `config1.yaml` to customize paths, groups, processing parameters, and ROI/tract definitions.

## Getting Help

- MRtrix3 Documentation: https://mrtrix.readthedocs.io/
- MRtrix3 Community Forum: https://community.mrtrix.org/
- PyBIDS Documentation: https://bids-standard.github.io/pybids/

