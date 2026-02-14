# GitHub Copilot Instructions for MBRAIN Project

## ⚠️ CRITICAL: ALWAYS READ config1.yaml FIRST - NO EXCEPTIONS! ⚠️

**STOP! Before doing ANYTHING in this project:**

1. **READ `/Users/glerma/soft/mbrain/config1.yaml` FIRST** - This is MANDATORY before ANY task.

2. **WHEN YOU FEEL LOST OR NEED TO SEARCH FOR FILES → READ config1.yaml**
   - Do NOT randomly search for files or guess paths
   - Do NOT try to open files by trial and error
   - ALL file paths, subjects, sessions, and parameters are in config1.yaml
   - If you're unsure about anything, the answer is in config1.yaml

3. **config1.yaml contains:**
   - All subject IDs and their sessions (lines ~550-590)
   - All file paths and folder structures
   - All tractography parameters
   - All ROI naming conventions
   - All pipeline instructions and rules

Before performing ANY task in this project, you MUST:

1. **Read `/Users/glerma/soft/mbrain/config1.yaml`** - This file contains ALL configuration, parameters, and instructions for this project.

2. **Check the `functions:` section** (line ~1060) - It documents ALL available scripts and their usage. DO NOT create new scripts if functionality exists!

3. **Follow the MBRAIN Instructions** (lines 108-127 of config1.yaml), especially:
   - Rule #7: DO NOT create new Python scripts - modify existing ones in `functions/`
   - Rule #9: For visualization, use `t2w_highres_flip.mif` as background image
   - Rule #11: Everything must be in `mrtrix_preproc_flip1` pipeline (never the old pipeline)
   - Rule #12: All parameters must come from config1.yaml, no hardcoded values

4. **Use the flip1 pipeline** at the path specified in config1.yaml `paths.derivatives_dir`

5. **Use functions/config_loader.py** to read parameters in Python scripts

## Project Structure:

```
/Users/glerma/soft/mbrain/
├── config1.yaml              # ALL configuration and documentation
├── README.md                # Project overview
├── functions/               # ALL Python scripts (14 files)
│   ├── __init__.py
│   ├── config_loader.py     # Config loading utility
│   ├── atlas_coregistration.py
│   ├── build_group_templates.py
│   ├── clean_tracts.py
│   ├── dwi_flip1_fod.py
│   ├── generate_tract_report.py
│   ├── generate_tracts.py
│   ├── matlab_clean_standalone.py
│   ├── qc_viewer.py
│   ├── resample_rois_to_fod.py
│   ├── segment_graft_slicewise.py
│   ├── segment_isocortex_slicewise.py
│   ├── tractography_segments.py
│   └── view_tracts.py
├── for_checking/            # Archived old/temporary files
└── data/derivatives/mrtrix_preproc_flip1/  # All processed data
```

## Quick Reference from config1.yaml:

- **Background image for visualization**: `t2w_highres_flip.mif`
- **FOD image**: `wm_fod_flip1.mif`
- **Mask**: `mask_flip1.mif`
- **ROI folders**: 
  - AP/APTR: `manual_roi/fod_space/` (flip1 transformed, resampled to FOD)
  - CTRL: `atlas/rois/fod_space/` (Allen atlas, resampled to FOD)
- **Segment folders**:
  - AP/APTR: `dwi/graft_segments_v2/fod_space/`
  - CTRL: `dwi/isocortex_segments_v2/fod_space/`
- **Environment**: `micromamba activate mbrain`

## Tractography Parameters (from config1.yaml):
- algorithm: iFOD2
- select: 2000
- step_size: 0.02
- angle: 60
- cutoff: 0.1
- minlength: 2
- maxlength: 15

## Key Scripts Usage:
```bash
# All scripts are in functions/ folder
python functions/tractography_segments.py --group CTRL --target thalamus
python functions/clean_tracts.py --subject sub-129-921-R1F
python functions/view_tracts.py --subject sub-596-C-OO-NmM
```

## Visualization Command Template:
```bash
mrview {base}/t2w_highres_flip.mif \
  -overlay.load {roi} -overlay.opacity 0.5 \
  -tractography.load {tract_file} &
```

## Never:
- Create new Python scripts (modify existing ones in functions/)
- Use old pipeline (mrtrix_preproc or mrtrix_preproc_flip2)
- Hardcode parameters
- Use `wm_fod_flip1.mif` or `b0_flip1.mif` as visualization background (use `t2w_highres_flip.mif`)

## ⚠️ CRITICAL: ALWAYS CHECK FILE DATES BEFORE USING PIPELINE OUTPUTS ⚠️

**Before using ANY derived file, check its modification timestamp against its upstream dependencies.**

The pipeline has a strict dependency order. If a downstream file is OLDER than its upstream source, it is STALE and must be regenerated:

### Pipeline Dependency Chain:
```
1. DWI preprocessing (dwi_flip_gradcheck_denoise_bias_anat.mif)
   ↓
2. Mask (mask_anat.mif) — via dwi2mask on _anat DWI
   ↓
3. FOD estimation (wm_fod_norm_group.mif) — via dwi2fod + mtnormalise with mask_anat
   ↓
4. Atlas coregistration (annotation_warped.nii.gz) — ANTs to T2w
   ↓
5. ROI extraction (atlas/rois/roi_*.nii.gz) — from warped atlas
   ↓
6. ROI resampling to FOD space (atlas/rois/fod_space/, manual_roi/fod_space/)
   — uses FOD as template via mrgrid regrid
   ↓
7. Fixels (fixels_group/) — via fod2fixel on FOD + mask_anat
   ↓
8. Wholebrain tractograms (wholebrain/wholebrain_1000k.tck) — via tckgen on FOD + mask_anat
   ↓
9. SIFT2 weights, passthrough tracts, connectivity analysis
```

### Rules:
1. **Before running any step N, verify that step N-1 outputs are NEWER than step N-2 outputs**
2. **If you find step N output is OLDER than step N-1 output → ASK the user before regenerating**
3. **Use `ls -lahT` to check timestamps. Compare dates, not sizes (old files may exist with different sizes)**
4. **When a step is re-run (e.g., new FOD), ALL downstream steps that depend on it are potentially stale**
5. **The FOD template used for `mrgrid regrid` (ROI resampling) MUST match the FOD used for analysis. Currently: `wm_fod_norm_group.mif`**
6. **After re-running atlas coregistration, ROI extraction AND fod_space resampling must both be re-run**

### Example check before fixel analysis:
```bash
# Verify: mask_anat (step 2) < FOD (step 3) < fixels (step 7)
ls -lahT subject/dwi/mask_anat.mif           # Should be oldest
ls -lahT subject/dwi/wm_fod_norm_group.mif   # Should be after mask
ls -lahT subject/dwi/fixels_group/            # Should be after FOD
ls -lahT subject/atlas/rois/fod_space/*.nii.gz  # Should be after BOTH atlas coreg AND FOD
```
