# Two-Stage Coregistration for ROI-based Tractography

## Overview

The pipeline implements a two-stage coregistration strategy optimized for high-resolution anatomical ROI drawing and tractography in native DWI space.

## Registration Strategy

### Stage 1: DWI → lowres T2w
- **Purpose**: Register diffusion to anatomical space
- **Rationale**: lowres T2w (200µm) matches DWI resolution for optimal alignment
- **Method**: ANTs rigid registration (6 DOF)
- **Output**: DWI in lowres T2w space

### Stage 2: lowres T2w → highres T2w
- **Purpose**: Link to high-resolution anatomical
- **Rationale**: highres T2w (100µm) provides best detail for ROI drawing
- **Method**: ANTs rigid registration (6 DOF)
- **Output**: Transform chain: highres → lowres → DWI

## Transform Chain

```
highres T2w (100µm) → lowres T2w (200µm) → DWI (200µm)
     [ROI drawn here]                      [Tractography here]
```

### Files Created

```
sub-XXX_ses-YYY_dwi2lowres.txt       # Stage 1: DWI → lowres
sub-XXX_ses-YYY_lowres2highres.txt   # Stage 2: lowres → highres
sub-XXX_ses-YYY_highres2lowres.txt   # Inverse of stage 2
sub-XXX_ses-YYY_highres2dwi.txt      # ⭐ Combined chain (use this!)
```

## Workflow

### 1. Run Coregistration

```bash
# Process all subjects
python main.py --preprocess --steps coreg

# Test on first subject
python main.py --preprocess --steps coreg --first-subject

# Dry-run to preview commands
python main.py --preprocess --steps coreg --first-subject --dry-run
```

**Output locations:**
- Transforms: `data/derivatives/mrtrix_preproc/sub-XXX/ses-YYY/dwi/`
- DWI in lowres space: `*_coreg.mif`
- T2w references: `*_t2w_lowres.mif`, `*_t2w_highres.mif`

### 2. Draw ROIs in Highres T2w

**Using mrview (MRtrix3):**
```bash
# Load highres T2w
mrview data/bids/sub-596-C-OO-NmM/ses-240911103349/anat/*_acq-highres_T2w.nii.gz

# In mrview:
# - Tool → ROI editor
# - Draw ROI
# - File → Save as... → roi_name_highres.mif
```

**Using ITK-SNAP or FSLeyes:**
- Load `*_acq-highres_T2w.nii.gz`
- Draw ROI using segmentation tools
- Save as NIfTI (will be auto-converted by pipeline)

**Tips:**
- Draw in highres for maximum precision
- Use anatomical landmarks visible in T2w
- ROI can be multi-label (different regions)
- Save with descriptive names: `motor_cortex_highres.mif`

### 3. Transform ROI to DWI Space

**Python API:**
```python
from utils import DWIPreprocessor
from pathlib import Path

preprocessor = DWIPreprocessor()

# Transform ROI using the highres → DWI chain
roi_dwi = preprocessor.transform_roi_to_dwi(
    subject_id='596-C-OO-NmM',
    session_id='240911103349',
    roi_highres=Path('motor_cortex_highres.mif')
)

print(f"ROI in DWI space: {roi_dwi}")
# Output: .../motor_cortex_dwi.mif
```

**Command line (manual):**
```bash
mrtransform \
  motor_cortex_highres.mif \
  motor_cortex_dwi.mif \
  -linear sub-XXX_ses-YYY_highres2dwi.txt \
  -template sub-XXX_ses-YYY_raw.mif \
  -interp nearest
```

### 4. Use ROI for Tractography

**Seed-based tractography:**
```bash
tckgen \
  sub-XXX_ses-YYY_fod_wm.mif \
  tracks_from_motor.tck \
  -seed_image motor_cortex_dwi.mif \
  -mask sub-XXX_ses-YYY_mask.mif \
  -select 100000 \
  -step 0.02 \
  -angle 60 \
  -minlength 1 \
  -maxlength 50
```

**Include/exclude ROIs:**
```bash
tckgen \
  sub-XXX_ses-YYY_fod_wm.mif \
  motor_to_cortex.tck \
  -seed_image motor_cortex_dwi.mif \
  -include cortex_dwi.mif \
  -exclude ventricles_dwi.mif \
  -select 10000
```

## Quality Control

### Check Registration Quality

**Overlay DWI on T2w:**
```bash
# Check DWI → lowres alignment
mrview sub-XXX_ses-YYY_t2w_lowres.mif \
  -overlay.load sub-XXX_ses-YYY_coreg.mif \
  -overlay.opacity 0.5
```

**Visual checks:**
- Brain boundaries aligned
- No major rotations or translations
- Tissue boundaries match
- Ventricles aligned

### Verify ROI Transformation

**Overlay transformed ROI on DWI:**
```bash
mrview sub-XXX_ses-YYY_raw.mif \
  -overlay.load motor_cortex_dwi.mif \
  -overlay.colourmap 2
```

**Check:**
- ROI in correct anatomical location
- No gaps or artifacts
- Labels preserved (for multi-label ROIs)

## Technical Details

### Registration Parameters

**Stage 1 (DWI → lowres):**
- Algorithm: `antsRegistrationSyN.sh`
- Transform: rigid (`-t r`)
- Dimensions: 3D (`-d 3`)
- Fixed: lowres T2w
- Moving: DWI b0 (mean)

**Stage 2 (lowres → highres):**
- Algorithm: `antsRegistrationSyN.sh`
- Transform: rigid (`-t r`)
- Dimensions: 3D (`-d 3`)
- Fixed: highres T2w
- Moving: lowres T2w

### Transform Composition

```python
# Inverse: lowres → highres becomes highres → lowres
transformcalc lowres2highres.txt invert highres2lowres.txt

# Compose chain: highres → lowres → DWI
transformcompose highres2lowres.txt dwi2lowres.txt highres2dwi.txt
```

### ROI Interpolation

- **Method**: Nearest-neighbor
- **Rationale**: Preserves discrete ROI labels
- **Template**: Native DWI space (ensures correct grid)

## File Organization

```
data/derivatives/mrtrix_preproc/
└── sub-596-C-OO-NmM/
    └── ses-240911103349/
        └── dwi/
            ├── sub-XXX_ses-YYY_raw.mif              # Native DWI
            ├── sub-XXX_ses-YYY_coreg.mif            # DWI in lowres space
            ├── sub-XXX_ses-YYY_t2w_lowres.mif       # Lowres T2w reference
            ├── sub-XXX_ses-YYY_t2w_highres.mif      # Highres T2w reference
            ├── sub-XXX_ses-YYY_dwi2lowres.txt       # Transform: DWI → lowres
            ├── sub-XXX_ses-YYY_lowres2highres.txt   # Transform: lowres → highres
            ├── sub-XXX_ses-YYY_highres2lowres.txt   # Transform: highres → lowres
            └── sub-XXX_ses-YYY_highres2dwi.txt      # ⭐ Transform chain
```

## Example: Motor Cortex to Spinal Cord

```python
from pathlib import Path
from utils import DWIPreprocessor

# Setup
preprocessor = DWIPreprocessor()
subject = '596-C-OO-NmM'
session = '240911103349'

# 1. Draw ROIs in highres (use mrview)
#    - motor_cortex_highres.mif
#    - spinal_cord_highres.mif

# 2. Transform both ROIs to DWI space
motor_dwi = preprocessor.transform_roi_to_dwi(
    subject, session, Path('motor_cortex_highres.mif')
)

spinal_dwi = preprocessor.transform_roi_to_dwi(
    subject, session, Path('spinal_cord_highres.mif')
)

# 3. Run tractography connecting the two regions
import subprocess
subprocess.run([
    'tckgen',
    f'data/derivatives/mrtrix_preproc/sub-{subject}/ses-{session}/dwi/*_fod_wm.mif',
    'motor_to_spinal.tck',
    '-seed_image', str(motor_dwi),
    '-include', str(spinal_dwi),
    '-select', '10000',
    '-step', '0.02',
    '-angle', '60'
])
```

## Advantages

### vs. Single-stage (DWI → highres)
- ✅ Better DWI alignment (matching resolutions)
- ✅ Less interpolation artifacts
- ✅ Smaller transform error

### vs. Drawing ROIs in lowres
- ✅ Higher anatomical precision (100µm vs 200µm)
- ✅ Better tissue contrast
- ✅ Easier to identify small structures

### vs. DWI space ROIs
- ✅ Superior anatomical detail
- ✅ Better tissue boundaries
- ✅ Independent of DWI quality

## Troubleshooting

### Registration Failed
```bash
# Check alignment visually
mrview *_t2w_lowres.mif -overlay.load *_b0.nii.gz

# Try with more iterations (edit dwi_preprocessor.py)
# Add: '-n', '8' to antsRegistrationSyN.sh calls
```

### ROI Not Aligned
```bash
# Verify transform chain
transforminfo *_highres2dwi.txt

# Check intermediate steps
mrtransform roi_highres.mif roi_lowres.mif -linear *_highres2lowres.txt
mrtransform roi_lowres.mif roi_dwi.mif -linear *_dwi2lowres.txt
```

### ROI Labels Lost
```bash
# Ensure nearest-neighbor interpolation
mrtransform ... -interp nearest

# Check ROI is integer type
mrinfo roi_highres.mif  # should show integer datatype
```

## References

- **rtp2-preproc**: Similar approach for human brain DWI
- **ANTs**: Registration tools (antsRegistrationSyN.sh)
- **MRtrix3**: Transform tools (transformconvert, transformcompose, mrtransform)
