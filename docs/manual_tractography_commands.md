# Manual Tractography Commands for Mouse Brain

This document provides step-by-step commands for performing tractography on preprocessed mouse brain DWI data.
These commands are designed for copy-paste execution in the terminal.

**Prerequisites:** Complete all steps in `manual_preprocessing_commands.md` first.

---

## PART 1: SETUP AND VARIABLES

### 1.1 Define paths (copy from preprocessing or re-run)

```bash
# Base paths
BIDS_DIR=/Users/glerma/soft/mbrain/data/bids
DERIV_DIR=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc
ATLAS_DIR=/Users/glerma/soft/mbrain/local/atlas/allen_ccf

# Subject info
SUBJECT=sub-596-C-OO-NmM
SESSION=ses-240911103349

# Output directory
OUT_DIR=$DERIV_DIR/$SUBJECT/$SESSION/dwi

# Verify preprocessed files exist
ls -la $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif
ls -la $OUT_DIR/mask.mif
```

---

## PART 2: TENSOR FITTING AND FA/MD MAPS

### 2.1 Fit diffusion tensor

```bash
# Fit tensor model to preprocessed DWI
dwi2tensor \
    -mask $OUT_DIR/mask.mif \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
    $OUT_DIR/tensor.mif -force

mrinfo $OUT_DIR/tensor.mif
```

### 2.2 Extract tensor metrics (FA, MD, AD, RD)

```bash
# Extract FA (Fractional Anisotropy)
tensor2metric $OUT_DIR/tensor.mif \
    -fa $OUT_DIR/fa.mif \
    -mask $OUT_DIR/mask.mif -force

# Extract MD (Mean Diffusivity)
tensor2metric $OUT_DIR/tensor.mif \
    -adc $OUT_DIR/md.mif \
    -mask $OUT_DIR/mask.mif -force

# Extract AD (Axial Diffusivity)
tensor2metric $OUT_DIR/tensor.mif \
    -ad $OUT_DIR/ad.mif \
    -mask $OUT_DIR/mask.mif -force

# Extract RD (Radial Diffusivity)
tensor2metric $OUT_DIR/tensor.mif \
    -rd $OUT_DIR/rd.mif \
    -mask $OUT_DIR/mask.mif -force

# List all metric files
ls -la $OUT_DIR/*.mif | grep -E "(fa|md|ad|rd|tensor)"
```

### 2.3 View tensor ellipsoids

```bash
# View tensor as ellipsoids overlaid on FA
mrview $OUT_DIR/fa.mif \
    -odf.load_tensor $OUT_DIR/tensor.mif &
```

### 2.4 View FA map

```bash
# View FA overlaid on T2w
mrview $OUT_DIR/t2w_highres_flip.nii.gz \
    -overlay.load $OUT_DIR/fa.mif \
    -overlay.opacity 0.5 &
```

### 2.5 View FA with ROIs

```bash
# View FA with atlas ROIs
mrview $OUT_DIR/fa.mif \
    -overlay.load $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_striatum.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_cerebellum.nii.gz &
```

---

## PART 3: CREATE FA-BASED WHITE MATTER MASK

### 3.1 Create WM mask from FA

```bash
# Threshold FA to create WM mask (FA > 0.2 for mouse brain)
mrthreshold $OUT_DIR/fa.mif \
    -abs 0.2 \
    $OUT_DIR/wm_mask.mif -force

# View WM mask
mrview $OUT_DIR/fa.mif \
    -overlay.load $OUT_DIR/wm_mask.mif \
    -overlay.opacity 0.5 &
```

### 3.2 Create GM/WM interface mask (for ACT if using)

```bash
# For mouse data without segmentation, use FA-based approach
# Higher FA threshold for core WM
mrthreshold $OUT_DIR/fa.mif -abs 0.3 $OUT_DIR/wm_core.mif -force

# Lower FA threshold for GM/WM boundary
mrthreshold $OUT_DIR/fa.mif -abs 0.15 $OUT_DIR/gm_wm_boundary.mif -force
```

---

## PART 4: RESPONSE FUNCTION ESTIMATION

### 4.1 Estimate response function (single-shell, single-tissue)

For single-shell data (like most mouse DWI), use dwi2response with tournier algorithm:

```bash
# Estimate WM response function
dwi2response tournier \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
    $OUT_DIR/response_wm.txt \
    -mask $OUT_DIR/mask.mif \
    -voxels $OUT_DIR/response_voxels.mif -force

# View response function voxels (should be in WM)
mrview $OUT_DIR/fa.mif \
    -overlay.load $OUT_DIR/response_voxels.mif \
    -overlay.opacity 0.8 &

# Display response function values
cat $OUT_DIR/response_wm.txt
```

### 4.2 Alternative: Manual response function (if automatic fails)

```bash
# For ex-vivo mouse brain, you may need to manually specify
# Check your b-value first
mrinfo $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif -shell_bvalues

# Typical single-shell response for mouse brain (adjust based on your data)
# Format: l=0, l=2, l=4 coefficients
echo "800 -200 20" > $OUT_DIR/response_wm_manual.txt
```

---

## PART 5: FOD ESTIMATION (FIBER ORIENTATION DISTRIBUTION)

### 5.1 Compute FODs (single-shell, single-tissue CSD)

```bash
# Constrained Spherical Deconvolution
dwi2fod csd \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
    $OUT_DIR/response_wm.txt \
    $OUT_DIR/fod_wm.mif \
    -mask $OUT_DIR/mask.mif -force

mrinfo $OUT_DIR/fod_wm.mif
```

### 5.2 View FODs

```bash
# View FODs overlaid on FA
mrview $OUT_DIR/fa.mif \
    -odf.load_sh $OUT_DIR/fod_wm.mif &
```

### 5.3 Normalize FODs (recommended before tractography)

```bash
# Normalize FOD intensities for better tractography
mtnormalise \
    $OUT_DIR/fod_wm.mif $OUT_DIR/fod_wm_norm.mif \
    -mask $OUT_DIR/mask.mif -force
```

---

## PART 6: TRACTOGRAPHY

### 6.1 Whole-brain tractography

```bash
# Generate whole-brain tractogram
# Adjust parameters for mouse brain size
tckgen \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/tracks_wholebrain.tck \
    -algorithm iFOD2 \
    -seed_image $OUT_DIR/wm_mask.mif \
    -mask $OUT_DIR/mask.mif \
    -select 100000 \
    -cutoff 0.05 \
    -minlength 1 \
    -maxlength 50 \
    -step 0.05 \
    -angle 45 \
    -force

# Check tractogram info
tckinfo $OUT_DIR/tracks_wholebrain.tck
```

### 6.2 View whole-brain tractogram

```bash
# View tracks on FA
mrview $OUT_DIR/fa.mif \
    -tractography.load $OUT_DIR/tracks_wholebrain.tck &
```

### 6.3 Subsample for visualization (optional)

```bash
# If too many tracks, subsample for faster visualization
tckedit $OUT_DIR/tracks_wholebrain.tck \
    $OUT_DIR/tracks_wholebrain_10k.tck \
    -number 10000 -force
```

---

## PART 7: ROI-BASED TRACTOGRAPHY

### 7.1 Seed from specific ROI (e.g., thalamus)

```bash
# Tractography seeded from thalamus
tckgen \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/tracks_from_thalamus.tck \
    -algorithm iFOD2 \
    -seed_image $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -mask $OUT_DIR/mask.mif \
    -select 10000 \
    -cutoff 0.05 \
    -minlength 1 \
    -maxlength 50 \
    -step 0.05 \
    -angle 45 \
    -force

mrview $OUT_DIR/fa.mif \
    -tractography.load $OUT_DIR/tracks_from_thalamus.tck \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz &
```

### 7.2 Connectivity between ROIs (e.g., thalamus to cortex)

```bash
# Tracks connecting thalamus and cortex (generate more initially for SIFT)
tckgen \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/tracks_thalamus_cortex_raw.tck \
    -algorithm iFOD2 \
    -seed_image $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -include $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
    -mask $OUT_DIR/mask.mif \
    -select 2000 \
    -cutoff 0.05 \
    -minlength 1 \
    -maxlength 20 \
    -step 0.05 \
    -angle 45 \
    -force

tckinfo $OUT_DIR/tracks_thalamus_cortex_raw.tck
```

### 7.3 Apply SIFT to thalamus-cortex tracks

```bash
# SIFT filtering - keep only most biologically plausible tracks
tcksift \
    $OUT_DIR/tracks_thalamus_cortex_raw.tck \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/tracks_thalamus_cortex_sift50.tck \
    -term_number 50 -force

tckinfo $OUT_DIR/tracks_thalamus_cortex_sift50.tck

# View SIFT-filtered thalamus-cortex tracks
mrview $OUT_DIR/fa.mif \
    -tractography.load $OUT_DIR/tracks_thalamus_cortex_sift50.tck \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_isocortex.nii.gz &
```

### 7.4 Extract tracks passing through hippocampus

```bash
# Filter existing tracks to those passing through hippocampus
tckedit $OUT_DIR/tracks_wholebrain.tck \
    $OUT_DIR/tracks_through_hippocampus.tck \
    -include $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz -force

tckinfo $OUT_DIR/tracks_through_hippocampus.tck

mrview $OUT_DIR/fa.mif \
    -tractography.load $OUT_DIR/tracks_through_hippocampus.tck \
    -overlay.load $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz &
```

---

## PART 8: SIFT FILTERING (SPHERICAL-DECONVOLUTION INFORMED FILTERING)

### 8.1 Apply SIFT to improve biological plausibility

```bash
# SIFT filtering - reduces tracks to biologically plausible subset
tcksift \
    $OUT_DIR/tracks_wholebrain.tck \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/tracks_wholebrain_sift.tck \
    -term_number 50000 -force

tckinfo $OUT_DIR/tracks_wholebrain_sift.tck

mrview $OUT_DIR/fa.mif \
    -tractography.load $OUT_DIR/tracks_wholebrain_sift.tck &
```

### 8.2 Alternative: SIFT2 for quantitative analysis

```bash
# SIFT2 - produces per-streamline weights instead of removing tracks
tcksift2 \
    $OUT_DIR/tracks_wholebrain.tck \
    $OUT_DIR/fod_wm_norm.mif \
    $OUT_DIR/sift2_weights.txt -force

# Check weights
head -20 $OUT_DIR/sift2_weights.txt
```

---

## PART 9: CONNECTOME GENERATION

### 9.1 Create connectome from atlas parcellation

```bash
# First, ensure annotation is in subject space
# (should already be done in preprocessing Part 3)

# Generate connectome matrix
tck2connectome \
    $OUT_DIR/tracks_wholebrain_sift.tck \
    $OUT_DIR/atlas/annotation_to_t2w.nii.gz \
    $OUT_DIR/connectome.csv \
    -symmetric \
    -zero_diagonal \
    -out_assignments $OUT_DIR/track_assignments.txt -force

# View connectome (if mrview supports it, or use Python)
echo "Connectome saved to $OUT_DIR/connectome.csv"
```

### 9.2 Visualize connectome with Python

```bash
# Quick visualization script
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

# Load connectome
conn = np.loadtxt("$OUT_DIR/connectome.csv", delimiter=",")

# Plot
plt.figure(figsize=(10, 8))
plt.imshow(np.log1p(conn), cmap='hot')
plt.colorbar(label='log(connections + 1)')
plt.title('Structural Connectome')
plt.savefig("$OUT_DIR/connectome.png", dpi=150)
plt.show()
EOF
```

---

## PART 10: TRACK DENSITY IMAGING (TDI)

### 10.1 Create track density image

```bash
# Create TDI from tractogram
tckmap \
    $OUT_DIR/tracks_wholebrain_sift.tck \
    $OUT_DIR/tdi.mif \
    -template $OUT_DIR/fa.mif -force

mrconvert $OUT_DIR/tdi.mif $OUT_DIR/tdi.nii.gz -force

# View TDI
mrview $OUT_DIR/tdi.nii.gz &
```

### 10.2 Create directionally-encoded color (DEC) TDI

```bash
# DEC-TDI - color-coded by direction
tckmap \
    $OUT_DIR/tracks_wholebrain_sift.tck \
    $OUT_DIR/tdi_dec.mif \
    -dec \
    -template $OUT_DIR/fa.mif -force

mrconvert $OUT_DIR/tdi_dec.mif $OUT_DIR/tdi_dec.nii.gz -force

mrview $OUT_DIR/tdi_dec.nii.gz &
```

---

## PART 11: VERIFY ALL OUTPUTS

### 11.1 List all tractography files

```bash
echo "=== Tensor/Metric Files ==="
ls -la $OUT_DIR/*.mif | grep -E "(tensor|fa|md|ad|rd|wm_mask|fod)"

echo ""
echo "=== Track Files ==="
ls -la $OUT_DIR/*.tck

echo ""
echo "=== Derived Images ==="
ls -la $OUT_DIR/*.nii.gz | grep -E "(fa|md|tdi|connectome)"
```

### 11.2 Summary visualization

```bash
# Final visualization: FA + FODs + tracks
mrview $OUT_DIR/fa.mif \
    -odf.load_sh $OUT_DIR/fod_wm_norm.mif \
    -tractography.load $OUT_DIR/tracks_wholebrain_sift.tck &
```

---

## Quick Reference: Tractography Parameters for Mouse Brain

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| `-step` | 0.05 mm | Smaller step for small brain |
| `-angle` | 45° | Max turning angle per step |
| `-cutoff` | 0.05-0.1 | FOD amplitude cutoff |
| `-minlength` | 1 mm | Minimum track length |
| `-maxlength` | 50 mm | Maximum track length |
| `-select` | 100k-1M | Number of tracks to generate |
| FA threshold | 0.2 | For WM mask in mouse |

---

## Troubleshooting

### If dwi2response fails:
```bash
# Check if there are enough high-FA voxels
mrstats $OUT_DIR/fa.mif -mask $OUT_DIR/mask.mif
# Try manual response function (see 4.2)
```

### If tractography produces no tracks:
```bash
# Check FOD amplitudes
mrstats $OUT_DIR/fod_wm_norm.mif -mask $OUT_DIR/mask.mif
# Lower cutoff value or check mask coverage
```

### If tracks don't follow expected anatomy:
```bash
# View FODs to check fiber orientations
mrview $OUT_DIR/fa.mif -odf.load_sh $OUT_DIR/fod_wm.mif &
# Check b-vector orientations
mrinfo $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif -dwgrad
```

---

## PART 12: MULTI-SUBJECT TRACT COMPARISON

This section describes how to compare individual tracts (e.g., thalamus-cortex) across multiple mice
without requiring group-average tractography. All subjects are registered to the Allen CCF template space.

### 12.1 Define template paths

```bash
# Allen CCF template (common space for all subjects)
TEMPLATE_DIR=$ATLAS_DIR
TEMPLATE=$TEMPLATE_DIR/average_template_25.nii.gz
PARCELLATION=$TEMPLATE_DIR/annotation_25.nii.gz

# List of subjects to process
SUBJECTS="sub-596-C-OO-NmM sub-769-C-O5-R1M sub-769-C-OO-NmM sub-875-C-OO-NM"
```

### 12.2 Project individual tracts to template space

For each subject, create tract density maps in template space:

```bash
# For current subject (already processed)
tckmap \
    $OUT_DIR/tracks_thalamus_cortex_sift50.tck \
    $OUT_DIR/thalamus_cortex_density.mif \
    -template $TEMPLATE -force

mrconvert $OUT_DIR/thalamus_cortex_density.mif $OUT_DIR/thalamus_cortex_density.nii.gz -force
```

### 12.3 Process all subjects (loop)

```bash
# Loop through all subjects
for SUBJ in $SUBJECTS; do
    echo "Processing $SUBJ..."
    
    # Find session (first available)
    SESS=$(ls $DERIV_DIR/$SUBJ/ | head -1)
    SUBJ_DIR=$DERIV_DIR/$SUBJ/$SESS/dwi
    
    # Create tract density map in template space
    if [ -f "$SUBJ_DIR/tracks_thalamus_cortex_sift50.tck" ]; then
        tckmap \
            $SUBJ_DIR/tracks_thalamus_cortex_sift50.tck \
            $SUBJ_DIR/thalamus_cortex_density.mif \
            -template $TEMPLATE -force
        
        mrconvert $SUBJ_DIR/thalamus_cortex_density.mif \
            $SUBJ_DIR/thalamus_cortex_density.nii.gz -force
        
        echo "  Created density map for $SUBJ"
    else
        echo "  WARNING: No thalamus-cortex tracks found for $SUBJ"
    fi
done
```

### 12.4 Stack densities for group analysis

```bash
# Create output directory for group analysis
GROUP_DIR=$DERIV_DIR/group_analysis
mkdir -p $GROUP_DIR

# Collect all density maps
DENSITY_FILES=""
for SUBJ in $SUBJECTS; do
    SESS=$(ls $DERIV_DIR/$SUBJ/ | head -1)
    DENSITY_FILES="$DENSITY_FILES $DERIV_DIR/$SUBJ/$SESS/dwi/thalamus_cortex_density.nii.gz"
done

# Stack into 4D volume for group statistics
mrcat $DENSITY_FILES $GROUP_DIR/all_thalamus_density.mif -axis 3 -force

mrinfo $GROUP_DIR/all_thalamus_density.mif
```

### 12.5 Compute group means

```bash
# Compute mean tract density across all subjects
mrmath $GROUP_DIR/all_thalamus_density.mif mean \
    $GROUP_DIR/mean_thalamus_density.mif -axis 3 -force

# Compute standard deviation
mrmath $GROUP_DIR/all_thalamus_density.mif std \
    $GROUP_DIR/std_thalamus_density.mif -axis 3 -force

# Convert for visualization
mrconvert $GROUP_DIR/mean_thalamus_density.mif $GROUP_DIR/mean_thalamus_density.nii.gz -force
mrconvert $GROUP_DIR/std_thalamus_density.mif $GROUP_DIR/std_thalamus_density.nii.gz -force
```

### 12.6 Visualize group mean tracts

```bash
# View mean tract density on template
mrview $TEMPLATE \
    -overlay.load $GROUP_DIR/mean_thalamus_density.nii.gz \
    -overlay.opacity 0.7 &
```

---

## PART 13: CONNECTOME-BASED TRACT COMPARISON

Alternative approach using connectome weights for more robust comparison of sparse mouse tracts.

### 13.1 Generate connectome for individual tracts

```bash
# Create connectome from thalamus-cortex tracts
# Scale by streamline length and inverse node volume for comparability
tck2connectome \
    $OUT_DIR/tracks_thalamus_cortex_sift50.tck \
    $OUT_DIR/atlas/annotation_to_t2w.nii.gz \
    $OUT_DIR/thalamus_cortex_connectome.csv \
    -scale_length \
    -scale_invnodevol \
    -symmetric \
    -zero_diagonal -force

echo "Connectome saved to $OUT_DIR/thalamus_cortex_connectome.csv"
```

### 13.2 Process connectomes for all subjects

```bash
for SUBJ in $SUBJECTS; do
    echo "Processing connectome for $SUBJ..."
    
    SESS=$(ls $DERIV_DIR/$SUBJ/ | head -1)
    SUBJ_DIR=$DERIV_DIR/$SUBJ/$SESS/dwi
    
    if [ -f "$SUBJ_DIR/tracks_thalamus_cortex_sift50.tck" ]; then
        tck2connectome \
            $SUBJ_DIR/tracks_thalamus_cortex_sift50.tck \
            $SUBJ_DIR/atlas/annotation_to_t2w.nii.gz \
            $SUBJ_DIR/thalamus_cortex_connectome.csv \
            -scale_length \
            -scale_invnodevol \
            -symmetric \
            -zero_diagonal -force
        
        echo "  Created connectome for $SUBJ"
    fi
done
```

### 13.3 Compare connectomes with Python

```bash
python3 << 'EOF'
import numpy as np
import os
import glob

# Find all connectome files
deriv_dir = os.environ.get('DERIV_DIR', '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc')
connectome_files = glob.glob(f"{deriv_dir}/sub-*/ses-*/dwi/thalamus_cortex_connectome.csv")

print(f"Found {len(connectome_files)} connectome files")

# Load and stack connectomes
connectomes = []
subjects = []
for f in sorted(connectome_files):
    subj = f.split('/')[-4]
    conn = np.loadtxt(f, delimiter=',')
    connectomes.append(conn)
    subjects.append(subj)
    print(f"  {subj}: shape {conn.shape}, sum {conn.sum():.2f}")

if len(connectomes) > 1:
    # Stack and compute statistics
    conn_stack = np.stack(connectomes, axis=0)
    mean_conn = np.mean(conn_stack, axis=0)
    std_conn = np.std(conn_stack, axis=0)
    
    print(f"\nGroup mean connectome: sum = {mean_conn.sum():.2f}")
    print(f"Group std connectome: mean std = {std_conn.mean():.4f}")
    
    # Save group results
    group_dir = f"{deriv_dir}/group_analysis"
    os.makedirs(group_dir, exist_ok=True)
    np.savetxt(f"{group_dir}/mean_thalamus_connectome.csv", mean_conn, delimiter=',')
    np.savetxt(f"{group_dir}/std_thalamus_connectome.csv", std_conn, delimiter=',')
    print(f"\nSaved to {group_dir}/")
EOF
```

### 13.4 Visualize connectome comparison

```bash
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt
import os

deriv_dir = os.environ.get('DERIV_DIR', '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc')
group_dir = f"{deriv_dir}/group_analysis"

# Load group connectomes
mean_conn = np.loadtxt(f"{group_dir}/mean_thalamus_connectome.csv", delimiter=',')

# Plot
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(np.log1p(mean_conn), cmap='hot', aspect='equal')
plt.colorbar(im, label='log(connection strength + 1)')
ax.set_title('Mean Thalamus-Cortex Connectome')
ax.set_xlabel('Region')
ax.set_ylabel('Region')

plt.tight_layout()
plt.savefig(f"{group_dir}/mean_thalamus_connectome.png", dpi=150)
print(f"Saved to {group_dir}/mean_thalamus_connectome.png")
plt.show()
EOF
```

---

## PART 14: STATISTICAL COMPARISON BETWEEN GROUPS

For comparing two groups of mice (e.g., control vs. treatment).

### 14.1 Define groups

```bash
# Define subject groups
GROUP1="sub-596-C-OO-NmM sub-769-C-OO-NmM"  # e.g., Control
GROUP2="sub-769-C-O5-R1M sub-875-C-OO-NM"   # e.g., Treatment

GROUP_DIR=$DERIV_DIR/group_analysis
mkdir -p $GROUP_DIR
```

### 14.2 Create group density stacks

```bash
# Group 1 density stack
G1_FILES=""
for SUBJ in $GROUP1; do
    SESS=$(ls $DERIV_DIR/$SUBJ/ | head -1)
    G1_FILES="$G1_FILES $DERIV_DIR/$SUBJ/$SESS/dwi/thalamus_cortex_density.nii.gz"
done
mrcat $G1_FILES $GROUP_DIR/group1_thalamus_density.mif -axis 3 -force

# Group 2 density stack
G2_FILES=""
for SUBJ in $GROUP2; do
    SESS=$(ls $DERIV_DIR/$SUBJ/ | head -1)
    G2_FILES="$G2_FILES $DERIV_DIR/$SUBJ/$SESS/dwi/thalamus_cortex_density.nii.gz"
done
mrcat $G2_FILES $GROUP_DIR/group2_thalamus_density.mif -axis 3 -force
```

### 14.3 Compute group means for visualization

```bash
# Group 1 mean
mrmath $GROUP_DIR/group1_thalamus_density.mif mean \
    $GROUP_DIR/group1_mean_density.mif -axis 3 -force

# Group 2 mean
mrmath $GROUP_DIR/group2_thalamus_density.mif mean \
    $GROUP_DIR/group2_mean_density.mif -axis 3 -force

# Convert to NIfTI
mrconvert $GROUP_DIR/group1_mean_density.mif $GROUP_DIR/group1_mean_density.nii.gz -force
mrconvert $GROUP_DIR/group2_mean_density.mif $GROUP_DIR/group2_mean_density.nii.gz -force
```

### 14.4 Visualize group differences

```bash
# View both group means on template
mrview $TEMPLATE \
    -overlay.load $GROUP_DIR/group1_mean_density.nii.gz \
    -overlay.load $GROUP_DIR/group2_mean_density.nii.gz &

# Compute difference map
mrcalc $GROUP_DIR/group2_mean_density.mif \
    $GROUP_DIR/group1_mean_density.mif \
    -subtract $GROUP_DIR/group_difference.mif -force

mrconvert $GROUP_DIR/group_difference.mif $GROUP_DIR/group_difference.nii.gz -force

mrview $TEMPLATE \
    -overlay.load $GROUP_DIR/group_difference.nii.gz \
    -overlay.opacity 0.7 &
```

### 14.5 Statistical testing with Python

```bash
python3 << 'EOF'
import numpy as np
from scipy import stats
import nibabel as nib
import os

deriv_dir = os.environ.get('DERIV_DIR', '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc')
group_dir = f"{deriv_dir}/group_analysis"

# Load group density maps
g1 = nib.load(f"{group_dir}/group1_thalamus_density.nii.gz")
g2 = nib.load(f"{group_dir}/group2_thalamus_density.nii.gz")

g1_data = g1.get_fdata()
g2_data = g2.get_fdata()

print(f"Group 1 shape: {g1_data.shape}")
print(f"Group 2 shape: {g2_data.shape}")

# Voxelwise t-test (where both groups have non-zero values)
mask = (g1_data.mean(axis=3) > 0) & (g2_data.mean(axis=3) > 0)
print(f"Testing {mask.sum()} voxels")

t_map = np.zeros(g1_data.shape[:3])
p_map = np.ones(g1_data.shape[:3])

for idx in zip(*np.where(mask)):
    i, j, k = idx
    t, p = stats.ttest_ind(g1_data[i,j,k,:], g2_data[i,j,k,:])
    t_map[i,j,k] = t
    p_map[i,j,k] = p

# Save results
t_img = nib.Nifti1Image(t_map, g1.affine)
p_img = nib.Nifti1Image(p_map, g1.affine)
nib.save(t_img, f"{group_dir}/tstat_map.nii.gz")
nib.save(p_img, f"{group_dir}/pval_map.nii.gz")

# Report significant voxels
sig_voxels = (p_map < 0.05) & mask
print(f"\nSignificant voxels (p<0.05): {sig_voxels.sum()}")
print(f"Saved t-stat and p-value maps to {group_dir}/")
EOF
```

### 14.6 View statistical results

```bash
# View t-statistic map
mrview $TEMPLATE \
    -overlay.load $GROUP_DIR/tstat_map.nii.gz \
    -overlay.opacity 0.7 &
```
