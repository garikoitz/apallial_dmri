# GitHub Copilot Instructions for MBRAIN Project

## CRITICAL: Always Read config.yaml First!

Before performing ANY task in this project, you MUST:

1. **Read `/Users/glerma/soft/mbrain/config.yaml`** - This file contains ALL configuration, parameters, and instructions for this project.

2. **Check the `functions:` section** (line ~1060) - It documents ALL available scripts and their usage. DO NOT create new scripts if functionality exists!

3. **Follow the MBRAIN Instructions** (lines 10-27 of config.yaml), especially:
   - Rule #7: DO NOT create new Python scripts - modify existing ones in `functions/`
   - Rule #9: For visualization, use `t2w_highres_flip.mif` as background image
   - Rule #11: Everything must be in `mrtrix_preproc_flip2` pipeline (never the old pipeline)
   - Rule #12: All parameters must come from config.yaml, no hardcoded values

4. **Use the flip2 pipeline** at `/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2/`

5. **Use functions/config_loader.py** to read parameters in Python scripts

## Project Structure:

```
/Users/glerma/soft/mbrain/
├── config.yaml              # ALL configuration and documentation
├── README.md                # Project overview
├── functions/               # ALL Python scripts (14 files)
│   ├── __init__.py
│   ├── config_loader.py     # Config loading utility
│   ├── atlas_coregistration.py
│   ├── build_group_templates.py
│   ├── clean_tracts.py
│   ├── dwi_flip2_fod.py
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
└── data/derivatives/mrtrix_preproc_flip2/  # All processed data
```

## Quick Reference from config.yaml:

- **Background image for visualization**: `t2w_highres_flip.mif`
- **FOD image**: `wm_fod_flip2.mif`
- **Mask**: `mask_flip2.mif`
- **ROI folders**: 
  - AP/APTR: `manual_roi/fod_space/` (flip2 transformed, resampled to FOD)
  - CTRL: `atlas/rois/fod_space/` (Allen atlas, resampled to FOD)
- **Segment folders**:
  - AP/APTR: `dwi/graft_segments_v2/fod_space/`
  - CTRL: `dwi/isocortex_segments_v2/fod_space/`
- **Environment**: `micromamba activate mbrain`

## Tractography Parameters (from config.yaml):
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
- Use old pipeline (mrtrix_preproc)
- Hardcode parameters
- Use `wm_fod_flip2.mif` or `b0_flip2.mif` as visualization background (use `t2w_highres_flip.mif`)
