# Manual Preprocessing Commands for Mouse Brain DTI

This document provides copy-paste commands to test the preprocessing pipeline step by step.

## Setup Variables

```bash
# Set paths
MBRAIN=/Users/glerma/soft/mbrain
BIDS_DIR=$MBRAIN/data/bids
DERIV_DIR=$MBRAIN/data/derivatives/mrtrix_preproc
ATLAS_DIR=$MBRAIN/local/atlas/allen_ccf

# Choose subject/session
SUBJECT=sub-596-C-OO-NmM
SESSION=ses-240911103349

# Create output directory
OUT_DIR=$DERIV_DIR/$SUBJECT/$SESSION/dwi
mkdir -p $OUT_DIR
cd $OUT_DIR
```

---

## PART 1: ANATOMICAL T2w PROCESSING

### 1.1 Check available T2w files

```bash
ls -la $BIDS_DIR/$SUBJECT/$SESSION/anat/
```

### 1.2 Define T2w file paths

```bash
T2W_LOWRES=$BIDS_DIR/$SUBJECT/$SESSION/anat/${SUBJECT}_${SESSION}_acq-lowres_T2w.nii.gz
T2W_HIGHRES=$BIDS_DIR/$SUBJECT/$SESSION/anat/${SUBJECT}_${SESSION}_acq-highres_T2w.nii.gz

# Verify files exist
ls -la $T2W_LOWRES
ls -la $T2W_HIGHRES
```

### 1.3 Convert T2w LOWRES to MIF

```bash
mrconvert $T2W_LOWRES $OUT_DIR/t2w_lowres.mif -force
mrinfo $OUT_DIR/t2w_lowres.mif
```

### 1.4 Flip T2w LOWRES (axis permutation + flip)

```bash
# Step 1: Axis permutation (for 3D: 2,1,0)
mrconvert $OUT_DIR/t2w_lowres.mif $OUT_DIR/t2w_lowres_perm.mif \
    -axes 2,1,0 -strides 1,2,3 -force

# Step 2: Axis flip
mrtransform $OUT_DIR/t2w_lowres_perm.mif $OUT_DIR/t2w_lowres_flip.mif \
    -flip 1 -force

# Clean up
rm $OUT_DIR/t2w_lowres_perm.mif

# Check result
mrinfo $OUT_DIR/t2w_lowres_flip.mif
```

### 1.5 Convert T2w HIGHRES to MIF

```bash
mrconvert $T2W_HIGHRES $OUT_DIR/t2w_highres.mif -force
mrinfo $OUT_DIR/t2w_highres.mif
```

### 1.6 Flip T2w HIGHRES (axis permutation + flip)

```bash
# Step 1: Axis permutation
mrconvert $OUT_DIR/t2w_highres.mif $OUT_DIR/t2w_highres_perm.mif \
    -axes 2,1,0 -strides 1,2,3 -force

# Step 2: Axis flip
mrtransform $OUT_DIR/t2w_highres_perm.mif $OUT_DIR/t2w_highres_flip.mif \
    -flip 1 -force

# Clean up
rm $OUT_DIR/t2w_highres_perm.mif

# Check result
mrinfo $OUT_DIR/t2w_highres_flip.mif
```

### 1.7 View T2w images to verify orientation

```bash
# Open in mrview to check orientation
mrview $OUT_DIR/t2w_highres_flip.mif &
mrview $OUT_DIR/t2w_lowres_flip.mif &
```

---

## PART 2: ATLAS TO T2w COREGISTRATION

### 2.1 Check atlas files

```bash
ls -la $ATLAS_DIR/
ATLAS_TEMPLATE=$ATLAS_DIR/average_template_25.nii.gz
ATLAS_ANNOTATION=$ATLAS_DIR/annotation_25.nii.gz

# If atlas not downloaded yet, download it:
# python -c "from utils.atlas_manager import AllenAtlasManager; AllenAtlasManager().download()"
```

### 2.2 View atlas template and annotation

```bash
# View atlas template with annotation overlay
mrview $ATLAS_TEMPLATE -overlay.load $ATLAS_ANNOTATION &
```

### 2.3 Extract complete ROIs from atlas (with ALL child structures)

**IMPORTANT:** Allen CCF uses hierarchical labels. A single ID (e.g., 549) is just one 
subnucleus, not the full structure. We must extract ALL descendant IDs to get complete ROIs.

```bash
# Create master ROI directory in atlas folder
mkdir -p $ATLAS_DIR/rois

# Extract complete ROIs including all child structures
python3 << 'EOF'
import json
import numpy as np
import nibabel as nib
from pathlib import Path

atlas_dir = Path("/Users/glerma/soft/mbrain/local/atlas/allen_ccf")
structures_file = atlas_dir / "structures.json"
annotation_file = atlas_dir / "annotation_25.nii.gz"
rois_dir = atlas_dir / "rois"
rois_dir.mkdir(exist_ok=True)

# ROIs to extract with their parent IDs
# From Allen CCF ontology:
# - Isocortex: 315
# - Striatum (Caudoputamen + others): 477
# - Thalamus: 549 is just VPM! Use 549's parent or search for "TH"
# - Hippocampal formation: 1089

ROIS = {
    'isocortex': 315,        # Isocortex (all cortical areas)
    'striatum': 477,         # Striatum  
    'thalamus': 549,         # This is wrong - need to find parent
    'hippocampus': 1089,     # Hippocampal formation
    'cerebellum': 512,       # Cerebellum
}

# Load structures
print("Loading structures.json...")
with open(structures_file) as f:
    data = json.load(f)

def get_all_descendant_ids(structures_data, parent_id):
    """Recursively get all descendant structure IDs including parent."""
    ids = set()
    
    def find_and_collect(node, target_id, collecting=False):
        current_id = node.get('id')
        
        if current_id == target_id:
            collecting = True
        
        if collecting:
            ids.add(current_id)
            
        for child in node.get('children', []):
            find_and_collect(child, target_id, collecting)
    
    if 'msg' in structures_data:
        for root in structures_data['msg']:
            find_and_collect(root, parent_id)
    
    return ids

def find_structure_by_acronym(structures_data, acronym):
    """Find structure ID by acronym."""
    def search(node):
        if node.get('acronym') == acronym:
            return node.get('id')
        for child in node.get('children', []):
            result = search(child)
            if result:
                return result
        return None
    
    if 'msg' in structures_data:
        for root in structures_data['msg']:
            result = search(root)
            if result:
                return result
    return None

# Find correct parent IDs
# Thalamus parent is "TH" not 549 (which is VPM - ventral posteromedial)
th_id = find_structure_by_acronym(data, 'TH')
print(f"Found Thalamus (TH) ID: {th_id}")

# Update ROIs with correct IDs
ROIS_CORRECTED = {
    'isocortex': 315,
    'striatum': 477,
    'thalamus': th_id if th_id else 549,  # Use TH parent, fallback to 549
    'hippocampus': 1089,
    'cerebellum': 512,
}

# Load annotation volume
print(f"Loading annotation: {annotation_file}")
annot_img = nib.load(annotation_file)
annot_data = annot_img.get_fdata().astype(int)

print(f"Annotation shape: {annot_data.shape}")
print(f"Unique labels in annotation: {len(np.unique(annot_data))}")

# Extract each ROI
for roi_name, parent_id in ROIS_CORRECTED.items():
    print(f"\nExtracting {roi_name} (parent ID: {parent_id})...")
    
    # Get all descendant IDs
    all_ids = get_all_descendant_ids(data, parent_id)
    print(f"  Found {len(all_ids)} structure IDs")
    
    if len(all_ids) == 0:
        print(f"  WARNING: No IDs found for {roi_name}!")
        continue
    
    # Create binary mask
    mask = np.isin(annot_data, list(all_ids)).astype(np.uint8)
    n_voxels = np.sum(mask)
    print(f"  Mask contains {n_voxels:,} voxels")
    
    if n_voxels == 0:
        print(f"  WARNING: Empty mask for {roi_name}!")
        continue
    
    # Save mask
    output_file = rois_dir / f"roi_{roi_name}.nii.gz"
    mask_img = nib.Nifti1Image(mask, annot_img.affine, annot_img.header)
    nib.save(mask_img, output_file)
    print(f"  Saved: {output_file}")

print("\n✓ All ROIs extracted to:", rois_dir)
print("\nThese master ROIs can now be warped to each subject's space.")
EOF
```

### 2.4 View complete ROIs on original atlas

```bash
# View the complete ROIs (with all child structures)
mrview $ATLAS_TEMPLATE \
    -overlay.load $ATLAS_DIR/rois/roi_isocortex.nii.gz \
    -overlay.load $ATLAS_DIR/rois/roi_striatum.nii.gz \
    -overlay.load $ATLAS_DIR/rois/roi_thalamus.nii.gz \
    -overlay.load $ATLAS_DIR/rois/roi_hippocampus.nii.gz &
```

### 2.5 Create atlas output directory

```bash
mkdir -p $OUT_DIR/atlas
```

### 2.3 Convert T2w to NIfTI for ANTs

```bash
# Use highres T2w as target for atlas registration
mrconvert $OUT_DIR/t2w_highres_flip.mif $OUT_DIR/t2w_highres_flip.nii.gz -force
```

### 2.4 Run ANTs registration: Atlas → T2w

```bash
# SyN registration (deformable) - atlas template to subject T2w
# Fixed = subject T2w (target space)
# Moving = atlas template (to be warped)

antsRegistrationSyN.sh \
    -d 3 \
    -t s \
    -f $OUT_DIR/t2w_highres_flip.nii.gz \
    -m $ATLAS_TEMPLATE \
    -o $OUT_DIR/atlas/atlas2t2w_
```

### 2.5 Apply transform to atlas template

```bash
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_TEMPLATE \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -o $OUT_DIR/atlas/atlas_template_warped.nii.gz
```

### 2.6 View atlas overlay on T2w to CHECK REGISTRATION QUALITY

```bash
# View warped atlas template overlaid on subject T2w
mrview $OUT_DIR/t2w_highres_flip.nii.gz \
    -overlay.load $OUT_DIR/atlas/atlas_template_warped.nii.gz \
    -overlay.opacity 0.5 &
```

---

## PART 3: WARP ATLAS ROIs TO SUBJECT SPACE

**NOTE:** We use the complete ROIs from `$ATLAS_DIR/rois/` (extracted with all child 
structures) and warp them to subject space. This ensures accurate ROI boundaries.

### 3.1 Create subject ROI directory

```bash
mkdir -p $OUT_DIR/atlas/rois
```

### 3.2 Warp master ROIs to subject T2w space

```bash
# Warp Isocortex ROI to subject space
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_DIR/rois/roi_isocortex.nii.gz \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -n NearestNeighbor \
    -o $OUT_DIR/atlas/rois/roi_isocortex.nii.gz

# Warp Striatum ROI to subject space
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_DIR/rois/roi_striatum.nii.gz \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -n NearestNeighbor \
    -o $OUT_DIR/atlas/rois/roi_striatum.nii.gz

# Warp Thalamus ROI to subject space
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_DIR/rois/roi_thalamus.nii.gz \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -n NearestNeighbor \
    -o $OUT_DIR/atlas/rois/roi_thalamus.nii.gz

# Warp Hippocampus ROI to subject space
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_DIR/rois/roi_hippocampus.nii.gz \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -n NearestNeighbor \
    -o $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz

# Warp Cerebellum ROI to subject space
antsApplyTransforms \
    -d 3 \
    -i $ATLAS_DIR/rois/roi_cerebellum.nii.gz \
    -r $OUT_DIR/t2w_highres_flip.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
    -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat \
    -n NearestNeighbor \
    -o $OUT_DIR/atlas/rois/roi_cerebellum.nii.gz
```

### 3.3 Check ROI voxel counts

```bash
echo "=== ROI Voxel Counts ==="
for roi in isocortex striatum thalamus hippocampus cerebellum; do
    count=$(mrstats $OUT_DIR/atlas/rois/roi_${roi}.nii.gz -mask $OUT_DIR/atlas/rois/roi_${roi}.nii.gz -output count 2>/dev/null | tail -1)
    echo "$roi: $count voxels"
done
```

### 3.4 View ROIs overlaid on T2w - VERIFY PLACEMENT

```bash
# View all ROIs on T2w
mrview $OUT_DIR/t2w_highres_flip.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_striatum.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz &
```

---

## PART 4: DWI PREPROCESSING

### 4.0 Set MRtrix configuration for b=0 threshold

```bash
# Mouse DWI often has b-values ~10-20 instead of exactly 0
# This ensures MRtrix treats b-values < 20 as b=0 volumes
echo "BZeroThreshold: 20" >> ~/.mrtrix.conf
```

### 4.1 Define DWI files

```bash
DWI_NII=$BIDS_DIR/$SUBJECT/$SESSION/dwi/${SUBJECT}_${SESSION}_dwi.nii.gz
DWI_BVEC=$BIDS_DIR/$SUBJECT/$SESSION/dwi/${SUBJECT}_${SESSION}_dwi.bvec
DWI_BVAL=$BIDS_DIR/$SUBJECT/$SESSION/dwi/${SUBJECT}_${SESSION}_dwi.bval

# Verify files exist
ls -la $DWI_NII $DWI_BVEC $DWI_BVAL
```

### 4.2 Convert DWI to MIF with gradients

```bash
mrconvert $DWI_NII $OUT_DIR/dwi.mif \
    -fslgrad $DWI_BVEC $DWI_BVAL -force
mrinfo $OUT_DIR/dwi.mif
```

### 4.3 Flip DWI (axis permutation + flip)

```bash
# Step 1: Axis permutation (for 4D: 2,1,0,3)
mrconvert $OUT_DIR/dwi.mif $OUT_DIR/dwi_perm.mif \
    -axes 2,1,0,3 -strides 1,2,3,4 -force

# Step 2: Axis flip
mrtransform $OUT_DIR/dwi_perm.mif $OUT_DIR/dwi_flip.mif \
    -flip 1 -reorient_fod no -force

# Clean up
rm $OUT_DIR/dwi_perm.mif

# Check result
mrinfo $OUT_DIR/dwi_flip.mif
```

### 4.4 Gradient check

```bash
dwigradcheck $OUT_DIR/dwi_flip.mif \
    -export_grad_mrtrix $OUT_DIR/corrected.b -force

mrconvert $OUT_DIR/dwi_flip.mif $OUT_DIR/dwi_flip_gradcheck.mif \
    -grad $OUT_DIR/corrected.b -force
```

### 4.5 Denoise

```bash
dwidenoise -extent 5,5,5 -noise $OUT_DIR/noise.mif -estimator Exp2 \
    $OUT_DIR/dwi_flip_gradcheck.mif $OUT_DIR/dwi_flip_gradcheck_denoise.mif -force
```

### 4.6 Create mask

```bash
dwi2mask $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -force

# Convert mask to NIfTI
mrconvert $OUT_DIR/mask.mif $OUT_DIR/mask.nii.gz -force
```

### 4.7 Create b0 reference images

```bash
# Extract mean b0 (use -config BZeroThreshold 20 for b-values < 20 to be treated as b=0)
dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise.mif - -bzero -config BZeroThreshold 20 -force | \
    mrmath - mean $OUT_DIR/b0_dwi.mif -axis 3 -force

# Convert to NIfTI
mrconvert $OUT_DIR/b0_dwi.mif $OUT_DIR/b0_dwi.nii.gz -force

# Convert mask to NIfTI
mrconvert $OUT_DIR/mask.mif $OUT_DIR/mask.nii.gz -force

# Apply mask to create brain-extracted b0
mrcalc $OUT_DIR/b0_dwi.nii.gz $OUT_DIR/mask.nii.gz -mult $OUT_DIR/b0_dwi_brain.nii.gz -force
```

### 4.8 Bias correction (ANTs)

```bash
# Note: -ants.c is [numberOfIterations,convergenceThreshold]
dwibiascorrect ants \
    -ants.b "[200,3]" \
    -ants.c "[10,1e-6]" \
    -ants.s 4 \
    -mask $OUT_DIR/mask.mif \
    $OUT_DIR/dwi_flip_gradcheck_denoise.mif \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
    -scratch $OUT_DIR/tmp_bias -force
```

---

## PART 5: DWI TO T2w COREGISTRATION (FINAL STEP)

### 5.1 Create b0 from bias-corrected DWI

```bash
# Extract mean b0 from bias-corrected DWI
dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif - -bzero -force | \
    mrmath - mean $OUT_DIR/b0_dwi_bias.mif -axis 3 -force

# Convert to NIfTI
mrconvert $OUT_DIR/b0_dwi_bias.mif $OUT_DIR/b0_dwi_bias.nii.gz -force

# Apply mask to create brain-extracted b0
mrcalc $OUT_DIR/b0_dwi_bias.nii.gz $OUT_DIR/mask.nii.gz -mult $OUT_DIR/b0_dwi_bias_brain.nii.gz -force
```

### 5.2 Run ANTs registration: DWI b0 → T2w

```bash
# Rigid registration - uses masked b0 from bias-corrected DWI
antsRegistrationSyN.sh \
    -d 3 \
    -t r \
    -f $OUT_DIR/t2w_highres_flip.nii.gz \
    -m $OUT_DIR/b0_dwi_bias_brain.nii.gz \
    -o $OUT_DIR/dwi2anat_
```

### 5.3 Convert ANTs transform to MRtrix format

```bash
ConvertTransformFile 3 \
    $OUT_DIR/dwi2anat_0GenericAffine.mat \
    $OUT_DIR/dwi2anat_affine.txt

transformconvert \
    $OUT_DIR/dwi2anat_affine.txt \
    itk_import \
    $OUT_DIR/dwi2anat_mrtrix.txt -force
```

### 5.4 Apply transform to full DWI

```bash
mrtransform \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
    $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
    -linear $OUT_DIR/dwi2anat_mrtrix.txt \
    -reorient_fod no -force
```

### 5.5 View final DWI overlaid with T2w and ROIs

```bash
# Extract b0 from final DWI
dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif - -bzero -force | \
    mrmath - mean $OUT_DIR/b0_final.mif -axis 3 -force
mrconvert $OUT_DIR/b0_final.mif $OUT_DIR/b0_final.nii.gz -force

# View overlay
mrview $OUT_DIR/t2w_highres_flip.nii.gz \
    -overlay.load $OUT_DIR/b0_final.nii.gz \
    -overlay.opacity 0.5 &
```

### 5.6 Verify ROI alignment with final DWI

The ROIs warped in Part 3 are in T2w space. Since DWI is now also in T2w space,
they should align. View final DWI with subject ROIs:

```bash
# View final DWI b0 with all subject ROIs
mrview $OUT_DIR/b0_final.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_striatum.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_cerebellum.nii.gz &
```

---

## PART 6: VERIFY EVERYTHING

### 6.1 List all output files

```bash
echo "=== Output Files ==="
ls -la $OUT_DIR/

echo ""
echo "=== Atlas Files ==="
ls -la $OUT_DIR/atlas/

echo ""
echo "=== ROI Files ==="
ls -la $OUT_DIR/atlas/rois/
```

### 6.2 Final visualization: DWI + T2w + ROIs

```bash
mrview $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
    -overlay.load $OUT_DIR/t2w_highres_flip.nii.gz \
    -overlay.opacity 0.3 \
    -overlay.load $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_striatum.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_hippocampus.nii.gz \
    -overlay.load $OUT_DIR/atlas/rois/roi_cerebellum.nii.gz &
```

---

## Quick Reference: File Naming Convention

| Step | Suffix | Example |
|------|--------|---------|
| Convert | (base) | `dwi.mif` |
| Flip | `_flip` | `dwi_flip.mif` |
| Gradcheck | `_gradcheck` | `dwi_flip_gradcheck.mif` |
| Denoise | `_denoise` | `dwi_flip_gradcheck_denoise.mif` |
| Degibbs | `_degibbs` | (optional) |
| Bias | `_bias` | `dwi_flip_gradcheck_denoise_bias.mif` |
| Norm | `_norm` | (optional) |
| Anat coreg | `_anat` | `dwi_flip_gradcheck_denoise_bias_anat.mif` |

---

## Troubleshooting

### If dwi2mask fails (common for ex-vivo data):
```bash
# Use threshold-based mask instead
mrthreshold $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -abs 10 -force
```

### If ANTs registration fails:
```bash
# Try with different options
antsRegistrationSyN.sh -d 3 -t r -n 4 -f $FIXED -m $MOVING -o $PREFIX
```

### Check image orientations:
```bash
mrinfo $OUT_DIR/t2w_highres_flip.mif -transform
mrinfo $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif -transform
```
