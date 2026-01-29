# Mouse Brain DWI Pipeline - All Subjects

This document provides commands to run the complete preprocessing and tractography pipeline
on all subjects with intermediate quality checks.

---

## QUICK START

```bash
# Run full pipeline on all subjects
cd /Users/glerma/soft/mbrain
./run_pipeline.sh

# Or run on specific subjects
./run_pipeline.sh sub-596-C-OO-NmM sub-769-C-O5-R1M
```

---

## PART 1: SETUP

### 1.1 Define global paths

```bash
# Base paths
export MBRAIN_DIR=/Users/glerma/soft/mbrain
export BIDS_DIR=$MBRAIN_DIR/data/bids
export DERIV_DIR=$MBRAIN_DIR/data/derivatives/mrtrix_preproc
export ATLAS_DIR=$MBRAIN_DIR/local/atlas/allen_ccf
export LOG_DIR=$MBRAIN_DIR/logs

# Create directories
mkdir -p $DERIV_DIR $LOG_DIR

# Atlas files
export ATLAS_TEMPLATE=$ATLAS_DIR/average_template_25.nii.gz
export ATLAS_ANNOTATION=$ATLAS_DIR/annotation_25.nii.gz
```

### 1.2 List all available subjects

```bash
# Find all subjects in BIDS directory
ls -d $BIDS_DIR/sub-* | xargs -n1 basename

# Count subjects
echo "Total subjects: $(ls -d $BIDS_DIR/sub-* | wc -l)"
```

### 1.3 Set MRtrix configuration

```bash
# Ensure BZeroThreshold is set (run once)
grep -q "BZeroThreshold" ~/.mrtrix.conf 2>/dev/null || \
    echo "BZeroThreshold: 20" >> ~/.mrtrix.conf

cat ~/.mrtrix.conf
```

---

## PART 2: SINGLE SUBJECT PIPELINE FUNCTION

### 2.1 Define the processing function

```bash
process_subject() {
    local SUBJECT=$1
    local SESSION=$2
    
    echo "========================================"
    echo "Processing: $SUBJECT / $SESSION"
    echo "Started: $(date)"
    echo "========================================"
    
    # Paths
    local BIDS_SUBJ=$BIDS_DIR/$SUBJECT/$SESSION
    local OUT_DIR=$DERIV_DIR/$SUBJECT/$SESSION/dwi
    local LOG_FILE=$LOG_DIR/${SUBJECT}_${SESSION}_$(date +%Y%m%d_%H%M%S).log
    
    mkdir -p $OUT_DIR/atlas/rois
    
    # Start logging
    exec > >(tee -a "$LOG_FILE") 2>&1
    
    #===========================================
    # STEP 1: T2w ANATOMICAL PROCESSING
    #===========================================
    echo "[STEP 1/10] T2w anatomical processing..."
    
    # Find T2w files
    T2W_HIGHRES=$(find $BIDS_SUBJ/anat -name "*T2w*.nii*" | grep -i high | head -1)
    T2W_LOWRES=$(find $BIDS_SUBJ/anat -name "*T2w*.nii*" | grep -iv high | head -1)
    
    if [ -z "$T2W_HIGHRES" ]; then
        T2W_HIGHRES=$(find $BIDS_SUBJ/anat -name "*T2w*.nii*" | head -1)
    fi
    
    if [ -n "$T2W_HIGHRES" ]; then
        mrconvert $T2W_HIGHRES $OUT_DIR/t2w_highres.mif -force
        mrtransform $OUT_DIR/t2w_highres.mif $OUT_DIR/t2w_highres_flip.mif \
            -flip 0,1 -force
        mrconvert $OUT_DIR/t2w_highres_flip.mif $OUT_DIR/t2w_highres_flip.nii.gz -force
        echo "  ✓ T2w highres processed"
    else
        echo "  ✗ WARNING: No T2w highres found"
    fi
    
    if [ -n "$T2W_LOWRES" ]; then
        mrconvert $T2W_LOWRES $OUT_DIR/t2w_lowres.mif -force
        mrtransform $OUT_DIR/t2w_lowres.mif $OUT_DIR/t2w_lowres_flip.mif \
            -flip 0,1 -force
        echo "  ✓ T2w lowres processed"
    fi
    
    #===========================================
    # STEP 2: ATLAS TO T2w COREGISTRATION
    #===========================================
    echo "[STEP 2/10] Atlas to T2w coregistration..."
    
    if [ -f "$OUT_DIR/t2w_highres_flip.nii.gz" ]; then
        antsRegistrationSyN.sh \
            -d 3 \
            -t s \
            -f $OUT_DIR/t2w_highres_flip.nii.gz \
            -m $ATLAS_TEMPLATE \
            -o $OUT_DIR/atlas/atlas2t2w_ \
            -n 4
        
        # Warp annotation to subject space
        antsApplyTransforms \
            -d 3 \
            -i $ATLAS_ANNOTATION \
            -r $OUT_DIR/t2w_highres_flip.nii.gz \
            -o $OUT_DIR/atlas/annotation_to_t2w.nii.gz \
            -n NearestNeighbor \
            -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
            -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat
        
        echo "  ✓ Atlas coregistration complete"
    else
        echo "  ✗ SKIPPING: No T2w for atlas registration"
    fi
    
    #===========================================
    # STEP 3: EXTRACT ROIs
    #===========================================
    echo "[STEP 3/10] Extracting ROIs..."
    
    if [ -f "$OUT_DIR/atlas/annotation_to_t2w.nii.gz" ]; then
        # Extract major ROIs using Python
        python3 << PYEOF
import nibabel as nib
import numpy as np
import json
import os

# Load annotation
ann = nib.load("$OUT_DIR/atlas/annotation_to_t2w.nii.gz")
ann_data = ann.get_fdata()

# Load structures
with open("$ATLAS_DIR/structures.json") as f:
    structures = json.load(f)

def get_all_children(structures, parent_acronym):
    """Get all descendant structure IDs for a parent acronym."""
    id_to_children = {}
    acronym_to_id = {}
    for s in structures:
        sid = s['id']
        parent = s.get('parent_structure_id')
        acronym_to_id[s['acronym']] = sid
        if parent not in id_to_children:
            id_to_children[parent] = []
        id_to_children[parent].append(sid)
    
    if parent_acronym not in acronym_to_id:
        return []
    
    parent_id = acronym_to_id[parent_acronym]
    all_ids = [parent_id]
    queue = [parent_id]
    while queue:
        current = queue.pop(0)
        children = id_to_children.get(current, [])
        all_ids.extend(children)
        queue.extend(children)
    return all_ids

# Define ROIs
rois = {
    'isocortex': 'Isocortex',
    'striatum': 'STR', 
    'thalamus': 'TH',
    'hippocampus': 'HIP',
    'cerebellum': 'CB'
}

out_dir = "$OUT_DIR/atlas/rois"
os.makedirs(out_dir, exist_ok=True)

for name, acronym in rois.items():
    ids = get_all_children(structures, acronym)
    if ids:
        mask = np.isin(ann_data, ids).astype(np.uint8)
        roi_img = nib.Nifti1Image(mask, ann.affine)
        nib.save(roi_img, f"{out_dir}/roi_{name}.nii.gz")
        print(f"  ✓ roi_{name}: {len(ids)} structures, {mask.sum()} voxels")
    else:
        print(f"  ✗ {name}: acronym '{acronym}' not found")
PYEOF
    else
        echo "  ✗ SKIPPING: No annotation for ROI extraction"
    fi
    
    #===========================================
    # STEP 4: DWI PREPROCESSING
    #===========================================
    echo "[STEP 4/10] DWI preprocessing..."
    
    DWI_NII=$(find $BIDS_SUBJ/dwi -name "*dwi.nii*" | head -1)
    DWI_BVEC=$(find $BIDS_SUBJ/dwi -name "*dwi.bvec" | head -1)
    DWI_BVAL=$(find $BIDS_SUBJ/dwi -name "*dwi.bval" | head -1)
    
    if [ -z "$DWI_NII" ]; then
        echo "  ✗ ERROR: No DWI file found"
        return 1
    fi
    
    # Convert
    mrconvert $DWI_NII $OUT_DIR/dwi.mif -fslgrad $DWI_BVEC $DWI_BVAL -force
    
    # Flip
    mrtransform $OUT_DIR/dwi.mif $OUT_DIR/dwi_flip.mif -flip 0,1 -force
    
    # Gradient check
    dwigradcheck $OUT_DIR/dwi_flip.mif \
        -export_grad_mrtrix $OUT_DIR/corrected_grads.txt -force
    mrconvert $OUT_DIR/dwi_flip.mif $OUT_DIR/dwi_flip_gradcheck.mif \
        -grad $OUT_DIR/corrected_grads.txt -force
    
    # Denoise
    dwidenoise -extent 5,5,5 \
        $OUT_DIR/dwi_flip_gradcheck.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise.mif -force
    
    echo "  ✓ DWI preprocessing complete"
    
    #===========================================
    # STEP 5: MASKING
    #===========================================
    echo "[STEP 5/10] Creating brain mask..."
    
    dwi2mask $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -force || \
        mrthreshold $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -abs 10 -force
    
    mrconvert $OUT_DIR/mask.mif $OUT_DIR/mask.nii.gz -force
    echo "  ✓ Mask created"
    
    #===========================================
    # STEP 6: B0 EXTRACTION AND BIAS CORRECTION
    #===========================================
    echo "[STEP 6/10] B0 extraction and bias correction..."
    
    # Extract b0
    dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise.mif - -bzero \
        -config BZeroThreshold 20 -force | \
        mrmath - mean $OUT_DIR/b0_dwi.mif -axis 3 -force
    
    mrconvert $OUT_DIR/b0_dwi.mif $OUT_DIR/b0_dwi.nii.gz -force
    
    # Bias correction
    dwibiascorrect ants \
        -ants.b "[200,3]" \
        -ants.c "[10,1e-6]" \
        -ants.s 4 \
        -mask $OUT_DIR/mask.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
        -scratch $OUT_DIR/tmp_bias -force
    
    rm -rf $OUT_DIR/tmp_bias
    echo "  ✓ Bias correction complete"
    
    #===========================================
    # STEP 7: DWI TO T2w COREGISTRATION
    #===========================================
    echo "[STEP 7/10] DWI to T2w coregistration..."
    
    if [ -f "$OUT_DIR/t2w_highres_flip.nii.gz" ]; then
        # Extract b0 from bias-corrected
        dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif - -bzero -force | \
            mrmath - mean $OUT_DIR/b0_dwi_bias.mif -axis 3 -force
        mrconvert $OUT_DIR/b0_dwi_bias.mif $OUT_DIR/b0_dwi_bias.nii.gz -force
        mrcalc $OUT_DIR/b0_dwi_bias.nii.gz $OUT_DIR/mask.nii.gz -mult \
            $OUT_DIR/b0_dwi_bias_brain.nii.gz -force
        
        # Register
        antsRegistrationSyN.sh \
            -d 3 \
            -t r \
            -f $OUT_DIR/t2w_highres_flip.nii.gz \
            -m $OUT_DIR/b0_dwi_bias_brain.nii.gz \
            -o $OUT_DIR/dwi2anat_
        
        # Convert and apply transform
        ConvertTransformFile 3 \
            $OUT_DIR/dwi2anat_0GenericAffine.mat \
            $OUT_DIR/dwi2anat_affine.txt
        
        transformconvert \
            $OUT_DIR/dwi2anat_affine.txt \
            itk_import \
            $OUT_DIR/dwi2anat_mrtrix.txt -force
        
        mrtransform \
            $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
            $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
            -linear $OUT_DIR/dwi2anat_mrtrix.txt \
            -reorient_fod no -force
        
        echo "  ✓ DWI to T2w coregistration complete"
    else
        echo "  ✗ SKIPPING: No T2w for DWI coregistration"
        cp $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
           $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif
    fi
    
    #===========================================
    # STEP 8: TENSOR AND FA
    #===========================================
    echo "[STEP 8/10] Tensor fitting and FA..."
    
    dwi2tensor \
        -mask $OUT_DIR/mask.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
        $OUT_DIR/tensor.mif -force
    
    tensor2metric $OUT_DIR/tensor.mif \
        -fa $OUT_DIR/fa.mif \
        -adc $OUT_DIR/md.mif \
        -mask $OUT_DIR/mask.mif -force
    
    echo "  ✓ Tensor and FA computed"
    
    #===========================================
    # STEP 9: FOD AND TRACTOGRAPHY
    #===========================================
    echo "[STEP 9/10] FOD estimation and tractography..."
    
    # WM mask from FA
    mrthreshold $OUT_DIR/fa.mif -abs 0.2 $OUT_DIR/wm_mask.mif -force
    
    # Response function
    dwi2response tournier \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
        $OUT_DIR/response_wm.txt \
        -mask $OUT_DIR/mask.mif \
        -voxels $OUT_DIR/response_voxels.mif -force
    
    # FOD
    dwi2fod csd \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
        $OUT_DIR/response_wm.txt \
        $OUT_DIR/fod_wm.mif \
        -mask $OUT_DIR/mask.mif -force
    
    # Normalize FOD
    mtnormalise \
        $OUT_DIR/fod_wm.mif $OUT_DIR/fod_wm_norm.mif \
        -mask $OUT_DIR/mask.mif -force
    
    # Thalamus-cortex tractography
    if [ -f "$OUT_DIR/atlas/rois/roi_thalamus.nii.gz" ]; then
        tckgen \
            $OUT_DIR/fod_wm_norm.mif \
            $OUT_DIR/tracks_thalamus_cortex_raw.tck \
            -algorithm iFOD2 \
            -seed_image $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
            -include $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
            -mask $OUT_DIR/mask.mif \
            -select 2000 \
            -cutoff 0.1 \
            -minlength 5 \
            -maxlength 25 \
            -step 0.05 \
            -angle 45 \
            -force
        
        # SIFT
        tcksift \
            $OUT_DIR/tracks_thalamus_cortex_raw.tck \
            $OUT_DIR/fod_wm_norm.mif \
            $OUT_DIR/tracks_thalamus_cortex_sift50.tck \
            -term_number 50 -force
        
        echo "  ✓ Tractography complete"
    else
        echo "  ✗ SKIPPING tractography: No ROIs"
    fi
    
    #===========================================
    # STEP 10: QUALITY CHECK REPORT
    #===========================================
    echo "[STEP 10/10] Generating QC report..."
    
    QC_FILE=$OUT_DIR/qc_report.txt
    
    cat > $QC_FILE << QCEOF
================================================================================
QUALITY CHECK REPORT
Subject: $SUBJECT
Session: $SESSION
Date: $(date)
================================================================================

1. INPUT FILES
--------------
T2w highres: $T2W_HIGHRES
DWI: $DWI_NII
b-vectors: $DWI_BVEC
b-values: $DWI_BVAL

2. DWI INFO
-----------
$(mrinfo $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif 2>/dev/null | head -20)

3. MASK STATISTICS
------------------
$(mrstats $OUT_DIR/mask.mif 2>/dev/null)

4. FA STATISTICS
----------------
$(mrstats $OUT_DIR/fa.mif -mask $OUT_DIR/mask.mif 2>/dev/null)

5. TRACTOGRAPHY
---------------
$(tckinfo $OUT_DIR/tracks_thalamus_cortex_sift50.tck 2>/dev/null || echo "No tracks generated")

6. OUTPUT FILES
---------------
$(ls -la $OUT_DIR/*.mif $OUT_DIR/*.nii.gz 2>/dev/null | wc -l) files generated

================================================================================
QCEOF
    
    cat $QC_FILE
    
    echo ""
    echo "========================================"
    echo "Completed: $SUBJECT / $SESSION"
    echo "Finished: $(date)"
    echo "Log: $LOG_FILE"
    echo "========================================"
}
```

---

## PART 3: RUN ALL SUBJECTS

### 3.1 Process all subjects sequentially

```bash
# Find and process all subjects
for SUBJ_PATH in $BIDS_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    
    # Find first session
    SESSION=$(ls $SUBJ_PATH | head -1)
    
    if [ -n "$SESSION" ]; then
        process_subject $SUBJECT $SESSION
    else
        echo "WARNING: No session found for $SUBJECT"
    fi
done
```

### 3.2 Process specific subjects

```bash
# List of subjects to process
SUBJECTS_TO_RUN="sub-596-C-OO-NmM sub-769-C-O5-R1M sub-769-C-OO-NmM sub-875-C-OO-NM"

for SUBJECT in $SUBJECTS_TO_RUN; do
    if [ -d "$BIDS_DIR/$SUBJECT" ]; then
        SESSION=$(ls $BIDS_DIR/$SUBJECT | head -1)
        process_subject $SUBJECT $SESSION
    else
        echo "WARNING: Subject $SUBJECT not found"
    fi
done
```

### 3.3 Process subjects in parallel (use with caution)

```bash
# Process up to 2 subjects at a time
MAX_JOBS=2

for SUBJ_PATH in $BIDS_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH | head -1)
    
    # Wait if too many jobs running
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 10
    done
    
    # Run in background
    process_subject $SUBJECT $SESSION &
done

# Wait for all to complete
wait
echo "All subjects processed!"
```

---

## PART 4: QUALITY CHECK ALL SUBJECTS

### 4.1 Summary of all QC reports

```bash
echo "==============================================="
echo "QC SUMMARY - ALL SUBJECTS"
echo "==============================================="

for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    QC_FILE=$SUBJ_PATH/$SESSION/dwi/qc_report.txt
    
    if [ -f "$QC_FILE" ]; then
        echo ""
        echo "--- $SUBJECT ---"
        
        # Check key files exist
        DWI_OK="✗"
        FA_OK="✗"
        TRACKS_OK="✗"
        
        [ -f "$SUBJ_PATH/$SESSION/dwi/dwi_flip_gradcheck_denoise_bias_anat.mif" ] && DWI_OK="✓"
        [ -f "$SUBJ_PATH/$SESSION/dwi/fa.mif" ] && FA_OK="✓"
        [ -f "$SUBJ_PATH/$SESSION/dwi/tracks_thalamus_cortex_sift50.tck" ] && TRACKS_OK="✓"
        
        echo "  DWI processed: $DWI_OK"
        echo "  FA computed: $FA_OK"
        echo "  Tractography: $TRACKS_OK"
        
        # FA stats
        if [ -f "$SUBJ_PATH/$SESSION/dwi/fa.mif" ]; then
            FA_MEAN=$(mrstats $SUBJ_PATH/$SESSION/dwi/fa.mif -mask $SUBJ_PATH/$SESSION/dwi/mask.mif -output mean 2>/dev/null)
            echo "  Mean FA: $FA_MEAN"
        fi
    else
        echo ""
        echo "--- $SUBJECT ---"
        echo "  ✗ Not processed"
    fi
done
```

### 4.2 Visual QC for all subjects

```bash
# Generate QC screenshots (requires mrview)
for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    OUT_DIR=$SUBJ_PATH/$SESSION/dwi
    
    if [ -f "$OUT_DIR/fa.mif" ]; then
        echo "Opening mrview for $SUBJECT..."
        mrview $OUT_DIR/fa.mif \
            -overlay.load $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
            -tractography.load $OUT_DIR/tracks_thalamus_cortex_sift50.tck &
        
        echo "Press Enter to continue to next subject..."
        read
        killall mrview 2>/dev/null
    fi
done
```

---

## PART 5: GROUP ANALYSIS

### 5.1 Create tract density maps for all subjects

```bash
TEMPLATE=$ATLAS_DIR/average_template_25.nii.gz
GROUP_DIR=$DERIV_DIR/group_analysis
mkdir -p $GROUP_DIR

for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    SUBJ_DIR=$SUBJ_PATH/$SESSION/dwi
    
    if [ -f "$SUBJ_DIR/tracks_thalamus_cortex_sift50.tck" ]; then
        echo "Creating density map for $SUBJECT..."
        
        tckmap \
            $SUBJ_DIR/tracks_thalamus_cortex_sift50.tck \
            $SUBJ_DIR/thalamus_cortex_density.mif \
            -template $TEMPLATE -force
        
        mrconvert $SUBJ_DIR/thalamus_cortex_density.mif \
            $SUBJ_DIR/thalamus_cortex_density.nii.gz -force
    fi
done
```

### 5.2 Stack and compute group statistics

```bash
# Collect all density maps
DENSITY_FILES=$(find $DERIV_DIR -name "thalamus_cortex_density.nii.gz" | sort | tr '\n' ' ')
N_SUBJ=$(echo $DENSITY_FILES | wc -w)

echo "Found $N_SUBJ subjects with density maps"

if [ $N_SUBJ -gt 1 ]; then
    # Stack
    mrcat $DENSITY_FILES $GROUP_DIR/all_thalamus_density.mif -axis 3 -force
    
    # Mean
    mrmath $GROUP_DIR/all_thalamus_density.mif mean \
        $GROUP_DIR/mean_thalamus_density.mif -axis 3 -force
    
    # Std
    mrmath $GROUP_DIR/all_thalamus_density.mif std \
        $GROUP_DIR/std_thalamus_density.mif -axis 3 -force
    
    # Convert
    mrconvert $GROUP_DIR/mean_thalamus_density.mif $GROUP_DIR/mean_thalamus_density.nii.gz -force
    mrconvert $GROUP_DIR/std_thalamus_density.mif $GROUP_DIR/std_thalamus_density.nii.gz -force
    
    echo "Group analysis complete!"
    echo "Mean density: $GROUP_DIR/mean_thalamus_density.nii.gz"
fi
```

### 5.3 View group results

```bash
mrview $ATLAS_TEMPLATE \
    -overlay.load $GROUP_DIR/mean_thalamus_density.nii.gz \
    -overlay.opacity 0.7 &
```

---

## PART 6: EXPORT RESULTS

### 6.1 Create summary CSV

```bash
OUTPUT_CSV=$DERIV_DIR/pipeline_summary.csv

echo "subject,session,dwi_ok,fa_ok,tracks_ok,fa_mean,fa_std,n_tracks" > $OUTPUT_CSV

for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    OUT_DIR=$SUBJ_PATH/$SESSION/dwi
    
    DWI_OK=0
    FA_OK=0
    TRACKS_OK=0
    FA_MEAN="NA"
    FA_STD="NA"
    N_TRACKS="NA"
    
    [ -f "$OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif" ] && DWI_OK=1
    
    if [ -f "$OUT_DIR/fa.mif" ]; then
        FA_OK=1
        FA_MEAN=$(mrstats $OUT_DIR/fa.mif -mask $OUT_DIR/mask.mif -output mean 2>/dev/null)
        FA_STD=$(mrstats $OUT_DIR/fa.mif -mask $OUT_DIR/mask.mif -output std 2>/dev/null)
    fi
    
    if [ -f "$OUT_DIR/tracks_thalamus_cortex_sift50.tck" ]; then
        TRACKS_OK=1
        N_TRACKS=$(tckinfo $OUT_DIR/tracks_thalamus_cortex_sift50.tck 2>/dev/null | grep "count" | awk '{print $2}')
    fi
    
    echo "$SUBJECT,$SESSION,$DWI_OK,$FA_OK,$TRACKS_OK,$FA_MEAN,$FA_STD,$N_TRACKS" >> $OUTPUT_CSV
done

echo "Summary saved to $OUTPUT_CSV"
cat $OUTPUT_CSV
```

---

## TROUBLESHOOTING

### Check failed subjects

```bash
for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    
    if [ ! -f "$SUBJ_PATH/$SESSION/dwi/fa.mif" ]; then
        echo "FAILED: $SUBJECT"
        
        # Check log
        LOG=$(ls -t $LOG_DIR/${SUBJECT}_*.log 2>/dev/null | head -1)
        if [ -f "$LOG" ]; then
            echo "  Last error:"
            grep -i "error\|failed\|warning" $LOG | tail -5
        fi
    fi
done
```

### Re-run failed subjects

```bash
FAILED_SUBJECTS=""

for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    
    if [ ! -f "$SUBJ_PATH/$SESSION/dwi/fa.mif" ]; then
        FAILED_SUBJECTS="$FAILED_SUBJECTS $SUBJECT"
    fi
done

echo "Re-running failed subjects: $FAILED_SUBJECTS"

for SUBJECT in $FAILED_SUBJECTS; do
    SESSION=$(ls $BIDS_DIR/$SUBJECT | head -1)
    process_subject $SUBJECT $SESSION
done
```
