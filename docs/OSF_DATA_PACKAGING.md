# MBRAIN Test Data Packaging Guide

This guide explains how to package and distribute the MBRAIN test data via OSF.io.

## Overview

The packaged data enables users to run the MBRAIN pipeline without lengthy preprocessing steps. It includes:
- **local/** — Allen CCF atlas + manual ROI delineations (~500 MB)
- **data/bids/** — Raw BIDS data for all 12 subjects (~1.5 GB)
- **data/derivatives/mrtrix_newmask/** — Pre-computed derivatives used by the notebook (~4 GB)

**Total size:** ~6 GB

## Step 1: Prepare the ZIP File

### 1.1 Gather Data

Create a staging directory and copy all required files:

```bash
# Create staging area
mkdir -p /tmp/mbrain_package
cd /tmp/mbrain_package

# Copy local directory (atlas + ROIs)
cp -r /path/to/apallial_dmri/local ./

# Copy BIDS data (all 12 subjects)
mkdir -p data/bids
cp -r /export/home/glerma/glerma/PROYECTOS/mbrain/data/bids/sub-* data/bids/
cp /export/home/glerma/glerma/PROYECTOS/mbrain/data/bids/dataset_description.json data/bids/
cp /export/home/glerma/glerma/PROYECTOS/mbrain/data/bids/participants.tsv data/bids/

# Copy pre-computed newmask derivatives (all 12 subjects)
mkdir -p data/derivatives/mrtrix_newmask
cp -r /export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/sub-* \
   data/derivatives/mrtrix_newmask/
cp /export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/*.json \
   data/derivatives/mrtrix_newmask/ 2>/dev/null || true
```

### 1.2 Verify Structure

```bash
# Check what we have
du -sh .   # Total size

# Verify critical files
ls -la local/atlas/                              # Allen templates
ls -la local/manual_roi_Slicer/                  # Slicer ROIs
ls -la local/manual_roi_flip/                    # Flipped ROIs

# Count BIDS subjects
ls -d data/bids/sub-* | wc -l                    # Should be 12

# Count derivatives
ls -d data/derivatives/mrtrix_newmask/sub-* | wc -l  # Should be 12
```

Expected file count:
```
local/atlas/: 5+ files (template, annotation, structures.json, rois/, etc.)
local/manual_roi_Slicer/: 8 .seg.nrrd files
local/manual_roi_flip/: 5 .mif files
data/bids/: 12 subject folders
data/derivatives/mrtrix_newmask/: 12 subject folders
```

### 1.3 Create ZIP

```bash
cd /tmp/mbrain_package

# Create ZIP with compression
zip -r -9 mbrain_test_data.zip local/ data/

# Verify ZIP integrity
unzip -t mbrain_test_data.zip | head -20
unzip -t mbrain_test_data.zip | tail -5  # Check final summary

# Check size
ls -lh mbrain_test_data.zip

# Calculate checksum (for documentation)
sha256sum mbrain_test_data.zip > mbrain_test_data.sha256
cat mbrain_test_data.sha256
```

Expected output:
```
adding: local/ (stored 0%)
adding: local/atlas/ (stored 0%)
adding: local/atlas/average_template_25.nii.gz (deflated 65%)
...
     0 extra bytes (100 compression = 0%)
Archive:  mbrain_test_data.zip
    testing: local/ ... OK
    testing: local/atlas/ ... OK
    ...

-rw-r--r-- 1 user staff 5.8G mbrain_test_data.zip  # Size should be ~5-6 GB
```

## Step 2: Upload to OSF.io

### 2.1 Create OSF Project

1. Go to https://osf.io/
2. Sign in (or create account if needed)
3. Click "Create a project"
4. Fill in:
   - **Title:** "MBRAIN Test Data — Mouse Brain Tractography"
   - **Description:** "Pre-processed DWI data and derivatives for the MBRAIN pipeline (12 subjects: 4 CTRL, 4 AP, 4 APTR). Enables reproducible analysis without lengthy preprocessing."
   - **Category:** "Neuroscience"
   - **License:** "CC-By Attribution 4.0 International"

### 2.2 Add Files to Project

1. In the project, click the "Files" tab
2. Click "Upload" or drag `mbrain_test_data.zip` into the window
3. Wait for upload to complete (may take 10-30 minutes depending on connection)
4. Once uploaded, click the file to get viewing options
5. Copy the download URL (should be like `https://osf.io/XXXXX/download`)

### 2.3 Get Permanent Download Link

Right-click the file → "Copy link address" or look for the permanent URL badge on the file page.

**Format:** `https://osf.io/XXXXX/download` where `XXXXX` is the unique OSF file ID.

## Step 3: Update Repository Configuration

### 3.1 Update setup_data.py

Edit `utils/setup_data.py` and replace the URL:

```python
# OLD:
DEFAULT_OSF_URL = "https://osf.io/XXXXX/download"  # ← Replace with actual OSF project ID

# NEW:
DEFAULT_OSF_URL = "https://osf.io/ABC123/download"  # Example
```

### 3.2 Document in README

Add to main `README.md`:

```markdown
## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/garikoitz/apallial_dmri.git
   cd apallial_dmri
   ```

2. The notebook will automatically download test data on first run (6 GB)
   - Alternatively, run: `python utils/setup_data.py`

3. Open the notebook:
   ```bash
   jupyter notebook notebooks/pipeline_groups_new_mask.py
   ```

4. Run Pre-Step -2 through -1c to regenerate from BIDS if desired, OR skip Pre-Steps and go to Step 0 when newmask data is already present
```

## Step 4: User Workflow

When a user clones the repo and runs the notebook:

### 4A: First Time (Auto-Download)

```python
# User runs notebook
# Pre-Step -2 detects missing local/ directory
# Automatically downloads mbrain_test_data.zip from OSF.io (~6 GB, ~20-40 min)
# Extracts to repo root
# User proceeds to PreStep -1a or Step 0
```

### 4B: Already Has Data

```python
# If local/ exists and newmask derivatives are present, skip to Step 0 (analysis)
```

### 4C: Manual Download (if auto fails)

```bash
# Download manually from OSF.io
# Extract to repo root:
unzip mbrain_test_data.zip

# Verify
python utils/setup_data.py --check-only
```

## Troubleshooting

### Download Times Out

```bash
# Increase timeout and retry
MBRAIN_DATA_URL="https://osf.io/ABC123/download" timeout 600 python utils/setup_data.py
```

### ZIP Corrupted After Download

```bash
# Verify checksum
sha256sum -c mbrain_test_data.sha256

# If fails, re-download from OSF
rm -f mbrain_test_data.zip
python utils/setup_data.py --force
```

### Insufficient Disk Space

```bash
# Check available space (need ~15 GB: 6 GB ZIP + 9 GB unzipped)
df -h ~

# If needed, download to external drive and extract there
# Then create a symlink
# Example: ln -s /Volumes/ExternalDrive/mbrain_data ~/mbrain_data
```

## Maintenance

### Updating Data

If preprocessing is re-run and outputs change:

```bash
# Create new ZIP (same as Step 1.3)
rm /tmp/mbrain_package/data/derivatives/mrtrix_newmask/*
cp -r /export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/sub-* \
   /tmp/mbrain_package/data/derivatives/mrtrix_newmask/

# Re-create ZIP and upload to OSF (create a new version)
zip -r -9 /tmp/mbrain_test_data_v2.zip /tmp/mbrain_package/*

# Update DEFAULT_OSF_URL in setup_data.py with new link
```

### Documentation for Paper

Include in paper's Methods section:

> "Test data (BIDS format + pre-computed mrtrix_newmask derivatives for 12 subjects) is available at https://osf.io/ABC123/. Data includes raw DWI, T2-weighted anatomical images, Allen CCF atlas, and manually delineated ROIs. Pre-computed derivatives (FOD estimates, brain masks, whole-brain tractograms) enable reproduction of downstream analysis without lengthy preprocessing (estimated 4-6 hours for preprocessing pipeline on a 40-core workstation)."

### Data Citation

Using the OSF citation generator, your full citation would be:

```bibtex
@misc{glerma_2026_mbrain,
  author = {Glerma, Garikoitz},
  title = {MBRAIN Test Data},
  year = {2026},
  url = {https://osf.io/ABC123/},
  publisher = {Open Science Framework}
}
```

## File Manifest

Quick reference for what's in the ZIP:

```
mbrain_test_data.zip (5.8 GB)
├── local/ (540 MB)
│   ├── atlas/
│   │   ├── allen_ccf/
│   │   │   ├── annotation_25.nii.gz
│   │   │   ├── average_template_25.nii.gz
│   │   │   ├── structures.json
│   │   │   ├── rois/
│   │   │   └── lut.txt
│   ├── manual_roi_Slicer/ (8 .seg.nrrd files)
│   └── manual_roi_flip/ (5 .mif files)
│
└── data/ (5.3 GB)
    ├── bids/
    │   ├── sub-129-921-R1F/
    │   ├── sub-129-921-R1L1F/
    │   ├── sub-129-1014-L1M/
    │   ├── sub-138-1010-L1F/
    │   ├── sub-596-C-OO-NmM/
    │   ├── sub-769-C-O5-R1M/
    │   ├── sub-769-C-OO-NmM/
    │   ├── sub-875-C-OO-NM/
    │   ├── sub-AP72-862-901-NmF/
    │   ├── sub-AP72-862-901-R1M/
    │   ├── sub-APtestmouseNmM/
    │   ├── sub-Apallial17NmF/
    │   ├── sub-Apallial20R1F/
    │   ├── dataset_description.json
    │   └── participants.tsv
    │
   └── derivatives/mrtrix_newmask/
        ├── sub-129-921-R1F/ses-*/dwi/ (t2w_flip, FODs, masks, ROIs, wholebrain.tck)
        ├── sub-129-921-R1L1F/ses-*/dwi/
        ├── ... (12 subjects total)
        └── passthrough_summary.json
```

---

**Questions?** Contact: garikoitz@gmail.com
