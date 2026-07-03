# OSF.io Data Package Workflow — Complete Checklist

This is a step-by-step checklist to prepare and distribute MBRAIN test data.

## 🎯 Summary

**What:** Package 12 subjects (BIDS + pre-computed mrtrix_newmask derivatives) into a single 6GB ZIP
**Where:** Upload to OSF.io for public download
**Why:** Enable reproducible research without lengthy preprocessing
**Who:** Anyone who wants to run the MBRAIN pipeline

---

## ✅ Step-by-Step Workflow

### Phase 1: Prepare Data (On Linux/Mac with data access)

#### [ ] 1.1 Verify you have all data
```bash
# Check that these directories exist and have content:
ls -la /export/home/glerma/glerma/PROYECTOS/mbrain/local/atlas/
ls -la /export/home/glerma/glerma/PROYECTOS/mbrain/local/manual_roi_Slicer/
ls -la /export/home/glerma/glerma/PROYECTOS/mbrain/data/bids/
ls -la /export/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask/
```
Expected: 12 subjects in bids/, 12 subjects in newmask derivatives

#### [ ] 1.2 Verify via script
```bash
cd /bcbl/home/home_g-m/glerma/soft/apallial_dmri
bash utils/create_data_package.sh --verify-only
```
Expected: All files found, $n_bids/$n_deriv match required counts

#### [ ] 1.3 Update script paths (if needed)
```bash
# Edit line 8 of utils/create_data_package.sh:
nano utils/create_data_package.sh
# Change: DATA_BASE="/export/home/glerma/glerma/PROYECTOS/mbrain"
# To match your actual data path
```

#### [ ] 1.4 Create ZIP
```bash
cd /bcbl/home/home_g-m/glerma/soft/apallial_dmri
bash utils/create_data_package.sh
# Output: /tmp/mbrain_test_data.zip (5-6 GB)
# Also creates: /tmp/mbrain_test_data.zip.sha256

# Wait for completion (~10-20 minutes on fast machine)
```

Expected output:
```
✅ ZIP created
  Size: 5.8G
✅ ZIP integrity verified
  ZIP contains 45,000+ files/folders
```

#### [ ] 1.5 Verify ZIP contents
```bash
# List key files to confirm they're inside
unzip -l /tmp/mbrain_test_data.zip | grep -E "(average_template|annotation_25|dataset_description)" | head -3
```

Expected: See Allen atlas files, BIDS metadata

---

### Phase 2: Upload to OSF.io (One-time setup)

#### [ ] 2.1 Create OSF.io Account (if needed)
- Go to https://osf.io/
- Sign up with email / institutional login

#### [ ] 2.2 Create New Project
- Click "Create a project"
- **Title:** `MBRAIN Test Data — Mouse Brain Tractography`
- **Description:** 
  ```
  Pre-processed diffusion MRI and anatomical data for all 12 subjects 
  (4 CTRL + 4 AP + 4 APTR) from the MBRAIN mouse brain tractography pipeline. 
  
  Contents:
  - Raw BIDS format: DWI and T2w anatomical images for all subjects
  - Pre-computed mrtrix_newmask derivatives: preprocessed DWI, FOD estimates, ROIs, whole-brain tractograms
  - Allen CCF atlas and manually delineated ROI delineations
  
  Enables reproducible analysis without lengthy preprocessing (4-6 hours per subject).
  ```
- **Category:** Neuroscience
- **License:** CC-By Attribution 4.0 International
- Click "Create"

#### [ ] 2.3 Upload ZIP to Project
- In the project page, click the "Files" tab
- Click "Upload" or drag `mbrain_test_data.zip`
- Wait for upload (10-30 minutes depending on connection)
- **Watch the progress bar** until 100%

#### [ ] 2.4 Get Download URL
Once upload is complete:
- Right-click the file → "Copy link address"
- OR click the file and look for a "Download" or "View" button
- Copy the URL (should look like `https://osf.io/v1d2e/download`)

#### [ ] 2.5 Document the URL
```
OSF Project ID: https://osf.io/XXXXX/
File URL: https://osf.io/XXXXX/download
File size: 5.8 GB
SHA256: [from .sha256 file]
```

---

### Phase 3: Update Repository (On development machine)

#### [ ] 3.1 Update setup_data.py
```bash
cd /bcbl/home/home_g-m/glerma/soft/apallial_dmri

# Edit utils/setup_data.py
nano utils/setup_data.py

# Find line ~20:
# DEFAULT_OSF_URL = "https://osf.io/XXXXX/download"

# Replace XXXXX with your actual OSF project ID (from Phase 2.5)
# Example:
# DEFAULT_OSF_URL = "https://osf.io/v1d2e/download"

# Save and exit (Ctrl+X, Y, Enter)
```

#### [ ] 3.2 Test auto-download on fresh machine (optional but recommended)
```bash
# On a different machine or in a fresh directory:
git clone https://github.com/garikoitz/apallial_dmri.git mbrain_test
cd mbrain_test

# This should auto-download
python utils/setup_data.py --check-only

# Expected: Downloads ZIP, extracts, validates all files found
```

#### [ ] 3.3 Commit and push changes
```bash
cd /bcbl/home/home_g-m/glerma/soft/apallial_dmri

git add -A
git commit -m "docs: add OSF data download configuration and setup scripts"
git push origin main
```

---

### Phase 4: Document for Paper

#### [ ] 4.1 Add to Methods section of paper
```markdown
**Data Availability:**

Test data (BIDS format + pre-computed mrtrix_newmask derivatives for 12 subjects)
is available at https://osf.io/XXXXX/. The dataset includes raw DWI 
and T2-weighted anatomical images, as well as preprocessed derivatives 
(brain masks, FOD estimates, whole-brain tractograms) that enable 
reproduction of downstream analyses without requiring lengthy preprocessing 
(estimated 4-6 hours for preprocessing on a 40-core workstation).

The complete MBRAIN pipeline is available at 
https://github.com/garikoitz/apallial_dmri. Users can run 
`python utils/setup_data.py` to automatically download test data, then execute 
`notebooks/pipeline_groups_new_mask.py` to reproduce all new-mask analyses.
```

#### [ ] 4.2 Add to Supplementary Materials (if applicable)
Create a supplementary document:
```
Supplementary Text 1: Data and Code Access

BIDS Data Structure:
  data/bids/sub-*/ses-*/anat/  (T2w images)
  data/bids/sub-*/ses-*/dwi/   (DWI + gradients)

Pre-computed Derivatives:
  data/derivatives/mrtrix_newmask/sub-*/ses-*/dwi/
    ├── t2w_highres_flip.mif (registered T2w)
    ├── dwi_flip_gradcheck_denoise_bias_anat.mif (final DWI)
    ├── wm_fod_norm_group.mif (normalized FODs)
    ├── mask_anat.mif (brain mask)
    ├── manual_roi/fod_space/ (ROIs in FOD space)
    └── wholebrain/ (whole-brain tractogram 1M streamlines)

Reference Files:
  local/atlas/ (Allen CCF atlas + structures)
  local/manual_roi_Slicer/ (original ROI delineations)
  local/manual_roi_flip/ (corrected ROI masks)
```

#### [ ] 4.3 Citation Format (for reviewers/readers)
```bibtex
@article{glerma2026mbrain,
  title={MBRAIN: Multi-subject Mouse Brain Rodent Analysis In Neuroimaging},
  author={Glerma, G. and ...},
  journal={...},
  year={2026}
}

@misc{glerma2026mbraindata,
  title={MBRAIN Test Data (12 subjects, BIDS + mrtrix_newmask derivatives)},
  author={Glerma, Garikoitz},
  year={2026},
  note={Open Science Framework},
  howpublished={\url{https://osf.io/XXXXX/}},
  doi={}  % Add OSF DOI if available
}
```

---

### Phase 5: Verify Everything Works

#### [ ] 5.1 Test from scratch on Mac
```bash
# Fresh directory on Mac
mkdir ~/test_mbrain && cd ~/test_mbrain

# Clone repo
git clone https://github.com/garikoitz/apallial_dmri.git
cd apallial_dmri

# Create symbolic links (simulating where data ends up)
ln -s ~/mbrain_data/data ./data
ln -s ~/mbrain_data/local ./local

# Run notebook (should auto-download on first cell)
jupyter notebook notebooks/pipeline_groups_new_mask.py

# Pre-Step -2 should detect missing local/ and auto-download
# Pre-Step -2 should then proceed to verify all data
```

Expected behavior:
```
[2026-07-02 14:30:00] ℹ️ LOCAL DATA NOT FOUND — AUTO-DOWNLOADING FROM OSF.io
[2026-07-02 14:30:02] ℹ️ Downloading from: https://osf.io/XXXXX/download
Progress: 45.2% (2.6 GB / 5.8 GB)
...
[2026-07-02 14:45:00] ✅ Setup complete! You can now run the notebook.
```

#### [ ] 5.2 Test skipping download (if data exists)
```bash
# Run again (should find local/ and skip download)
jupyter notebook notebooks/pipeline_groups_new_mask.py

# Should see:
# ✅ Local data directory found at ~/soft/mbrain/local
```

#### [ ] 5.3 Test Pre-Steps (optional, slow)
```bash
# Comment out osf url from config
# Run Pre-Step -1a, -1b, -1c to verify preprocessing pipeline works
# (This will take 12-36 hours for 12 subjects)
```

---

## 📋 Final Checklist

- [ ] ZIP created and verified (`/tmp/mbrain_test_data.zip`)
- [ ] Uploaded to OSF.io  
- [ ] OSF download URL obtained
- [ ] `utils/setup_data.py` updated with correct URL
- [ ] Pre-Steps -2 verified downloads work correctly
- [ ] README.md updated with link to SETUP.md
- [ ] Methods section of paper documents data availability with OSF link
- [ ] Paper includes data citation with DOI (if OSF generated one)
- [ ] Changes committed and pushed to GitHub
- [ ] Test on Mac confirmed auto-download works

---

## 🚀 You're Done!

Users can now:
1. Clone the repo: `git clone https://github.com/garikoitz/apallial_dmri.git`
2. Run the notebook: `jupyter notebook notebooks/pipeline_groups_new_mask.py`
3. Data auto-downloads automatically on first run
4. Pre-Steps verify everything is ready
5. Users choose: fast-track analysis OR full re-preprocessing

---

## 📞 Support

**If auto-download fails:**
- Manual download from: https://osf.io/XXXXX/
- Extract ZIP: `unzip mbrain_test_data.zip` in repo root
- Verify: `python utils/setup_data.py --check-only`

**For questions:**
- Contact: garikoitz@gmail.com
- Or: Open issue on GitHub
