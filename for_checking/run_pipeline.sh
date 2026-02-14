#!/bin/bash
#
# Mouse Brain DWI Pipeline
# Run preprocessing and tractography on all subjects
#
# Usage:
#   ./run_pipeline.sh                    # Process all subjects
#   ./run_pipeline.sh sub-001 sub-002    # Process specific subjects
#
# Author: mbrain pipeline
# Date: 2026-01-03
#

set -e  # Exit on error

#===========================================
# CONFIGURATION
#===========================================
export MBRAIN_DIR=/Users/glerma/soft/mbrain
export BIDS_DIR=$MBRAIN_DIR/data/bids
export DERIV_DIR=$MBRAIN_DIR/data/derivatives/mrtrix_preproc
export ATLAS_DIR=$MBRAIN_DIR/local/atlas/allen_ccf
export LOG_DIR=$MBRAIN_DIR/logs

export ATLAS_TEMPLATE=$ATLAS_DIR/average_template_25.nii.gz
export ATLAS_ANNOTATION=$ATLAS_DIR/annotation_25.nii.gz

# Create directories
mkdir -p $DERIV_DIR $LOG_DIR

# Ensure MRtrix config
grep -q "BZeroThreshold" ~/.mrtrix.conf 2>/dev/null || \
    echo "BZeroThreshold: 20" >> ~/.mrtrix.conf

#===========================================
# PROCESSING FUNCTION
#===========================================
process_subject() {
    local SUBJECT=$1
    local SESSION=$2
    
    echo ""
    echo "========================================"
    echo "Processing: $SUBJECT / $SESSION"
    echo "Started: $(date)"
    echo "========================================"
    
    # Paths
    local BIDS_SUBJ=$BIDS_DIR/$SUBJECT/$SESSION
    local OUT_DIR=$DERIV_DIR/$SUBJECT/$SESSION/dwi
    local LOG_FILE=$LOG_DIR/${SUBJECT}_${SESSION}_$(date +%Y%m%d_%H%M%S).log
    
    mkdir -p $OUT_DIR/atlas/rois
    
    # Log everything
    exec > >(tee -a "$LOG_FILE") 2>&1
    
    #-------------------------------------------
    # STEP 1: T2w ANATOMICAL
    #-------------------------------------------
    echo "[1/10] T2w anatomical processing..."
    
    T2W_HIGHRES=$(find $BIDS_SUBJ/anat -name "*T2w*.nii*" 2>/dev/null | grep -i high | head -1)
    [ -z "$T2W_HIGHRES" ] && T2W_HIGHRES=$(find $BIDS_SUBJ/anat -name "*T2w*.nii*" 2>/dev/null | head -1)
    
    if [ -n "$T2W_HIGHRES" ] && [ -f "$T2W_HIGHRES" ]; then
        mrconvert $T2W_HIGHRES $OUT_DIR/t2w_highres.mif -force
        mrtransform $OUT_DIR/t2w_highres.mif $OUT_DIR/t2w_highres_flip.mif -flip 0,1 -force
        mrconvert $OUT_DIR/t2w_highres_flip.mif $OUT_DIR/t2w_highres_flip.nii.gz -force
        echo "  ✓ T2w processed"
    else
        echo "  ⚠ No T2w found"
    fi
    
    #-------------------------------------------
    # STEP 2: ATLAS COREGISTRATION
    #-------------------------------------------
    echo "[2/10] Atlas coregistration..."
    
    if [ -f "$OUT_DIR/t2w_highres_flip.nii.gz" ]; then
        antsRegistrationSyN.sh -d 3 -t s \
            -f $OUT_DIR/t2w_highres_flip.nii.gz \
            -m $ATLAS_TEMPLATE \
            -o $OUT_DIR/atlas/atlas2t2w_ -n 4
        
        antsApplyTransforms -d 3 \
            -i $ATLAS_ANNOTATION \
            -r $OUT_DIR/t2w_highres_flip.nii.gz \
            -o $OUT_DIR/atlas/annotation_to_t2w.nii.gz \
            -n NearestNeighbor \
            -t $OUT_DIR/atlas/atlas2t2w_1Warp.nii.gz \
            -t $OUT_DIR/atlas/atlas2t2w_0GenericAffine.mat
        echo "  ✓ Atlas registered"
    else
        echo "  ⚠ Skipping atlas registration"
    fi
    
    #-------------------------------------------
    # STEP 3: EXTRACT ROIs
    #-------------------------------------------
    echo "[3/10] Extracting ROIs..."
    
    if [ -f "$OUT_DIR/atlas/annotation_to_t2w.nii.gz" ]; then
        python3 << PYEOF
import nibabel as nib
import numpy as np
import json
import os

ann = nib.load("$OUT_DIR/atlas/annotation_to_t2w.nii.gz")
ann_data = ann.get_fdata()

with open("$ATLAS_DIR/structures.json") as f:
    structures = json.load(f)

def get_children(structures, acronym):
    id_map = {s['acronym']: s['id'] for s in structures}
    children = {s.get('parent_structure_id'): [] for s in structures}
    for s in structures:
        p = s.get('parent_structure_id')
        if p in children:
            children[p].append(s['id'])
    
    if acronym not in id_map:
        return []
    
    ids = [id_map[acronym]]
    queue = list(ids)
    while queue:
        c = queue.pop(0)
        kids = children.get(c, [])
        ids.extend(kids)
        queue.extend(kids)
    return ids

rois = {'isocortex': 'Isocortex', 'striatum': 'STR', 'thalamus': 'TH', 
        'hippocampus': 'HIP', 'cerebellum': 'CB'}

for name, acr in rois.items():
    ids = get_children(structures, acr)
    if ids:
        mask = np.isin(ann_data, ids).astype(np.uint8)
        nib.save(nib.Nifti1Image(mask, ann.affine), 
                 f"$OUT_DIR/atlas/rois/roi_{name}.nii.gz")
        print(f"  ✓ roi_{name}: {mask.sum()} voxels")
PYEOF
    fi
    
    #-------------------------------------------
    # STEP 4: DWI PREPROCESSING
    #-------------------------------------------
    echo "[4/10] DWI preprocessing..."
    
    DWI_NII=$(find $BIDS_SUBJ/dwi -name "*dwi.nii*" 2>/dev/null | head -1)
    DWI_BVEC=$(find $BIDS_SUBJ/dwi -name "*dwi.bvec" 2>/dev/null | head -1)
    DWI_BVAL=$(find $BIDS_SUBJ/dwi -name "*dwi.bval" 2>/dev/null | head -1)
    
    if [ -z "$DWI_NII" ]; then
        echo "  ✗ ERROR: No DWI found"
        return 1
    fi
    
    mrconvert $DWI_NII $OUT_DIR/dwi.mif -fslgrad $DWI_BVEC $DWI_BVAL -force
    mrtransform $OUT_DIR/dwi.mif $OUT_DIR/dwi_flip.mif -flip 0,1 -force
    dwigradcheck $OUT_DIR/dwi_flip.mif -export_grad_mrtrix $OUT_DIR/corrected_grads.txt -force
    mrconvert $OUT_DIR/dwi_flip.mif $OUT_DIR/dwi_flip_gradcheck.mif -grad $OUT_DIR/corrected_grads.txt -force
    dwidenoise -extent 5,5,5 $OUT_DIR/dwi_flip_gradcheck.mif $OUT_DIR/dwi_flip_gradcheck_denoise.mif -force
    echo "  ✓ DWI preprocessed"
    
    #-------------------------------------------
    # STEP 5: MASKING
    #-------------------------------------------
    echo "[5/10] Creating mask..."
    
    dwi2mask $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -force 2>/dev/null || \
        mrthreshold $OUT_DIR/dwi_flip_gradcheck_denoise.mif $OUT_DIR/mask.mif -abs 10 -force
    mrconvert $OUT_DIR/mask.mif $OUT_DIR/mask.nii.gz -force
    echo "  ✓ Mask created"
    
    #-------------------------------------------
    # STEP 6: B0 AND BIAS CORRECTION
    #-------------------------------------------
    echo "[6/10] Bias correction..."
    
    dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise.mif - -bzero -config BZeroThreshold 20 -force | \
        mrmath - mean $OUT_DIR/b0_dwi.mif -axis 3 -force
    mrconvert $OUT_DIR/b0_dwi.mif $OUT_DIR/b0_dwi.nii.gz -force
    
    dwibiascorrect ants -ants.b "[200,3]" -ants.c "[10,1e-6]" -ants.s 4 \
        -mask $OUT_DIR/mask.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
        -scratch $OUT_DIR/tmp_bias -force
    rm -rf $OUT_DIR/tmp_bias
    echo "  ✓ Bias corrected"
    
    #-------------------------------------------
    # STEP 7: DWI TO T2w REGISTRATION
    #-------------------------------------------
    echo "[7/10] DWI to T2w registration..."
    
    if [ -f "$OUT_DIR/t2w_highres_flip.nii.gz" ]; then
        dwiextract $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif - -bzero -force | \
            mrmath - mean $OUT_DIR/b0_dwi_bias.mif -axis 3 -force
        mrconvert $OUT_DIR/b0_dwi_bias.mif $OUT_DIR/b0_dwi_bias.nii.gz -force
        mrcalc $OUT_DIR/b0_dwi_bias.nii.gz $OUT_DIR/mask.nii.gz -mult $OUT_DIR/b0_dwi_bias_brain.nii.gz -force
        
        antsRegistrationSyN.sh -d 3 -t r \
            -f $OUT_DIR/t2w_highres_flip.nii.gz \
            -m $OUT_DIR/b0_dwi_bias_brain.nii.gz \
            -o $OUT_DIR/dwi2anat_
        
        ConvertTransformFile 3 $OUT_DIR/dwi2anat_0GenericAffine.mat $OUT_DIR/dwi2anat_affine.txt
        transformconvert $OUT_DIR/dwi2anat_affine.txt itk_import $OUT_DIR/dwi2anat_mrtrix.txt -force
        mrtransform $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif \
            $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
            -linear $OUT_DIR/dwi2anat_mrtrix.txt -reorient_fod no -force
        echo "  ✓ DWI registered to T2w"
    else
        cp $OUT_DIR/dwi_flip_gradcheck_denoise_bias.mif $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif
        echo "  ⚠ No T2w - skipping registration"
    fi
    
    #-------------------------------------------
    # STEP 8: TENSOR AND FA
    #-------------------------------------------
    echo "[8/10] Computing tensor and FA..."
    
    dwi2tensor -mask $OUT_DIR/mask.mif \
        $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif $OUT_DIR/tensor.mif -force
    tensor2metric $OUT_DIR/tensor.mif -fa $OUT_DIR/fa.mif -adc $OUT_DIR/md.mif \
        -mask $OUT_DIR/mask.mif -force
    echo "  ✓ FA computed"
    
    #-------------------------------------------
    # STEP 9: FOD AND TRACTOGRAPHY
    #-------------------------------------------
    echo "[9/10] Tractography..."
    
    mrthreshold $OUT_DIR/fa.mif -abs 0.2 $OUT_DIR/wm_mask.mif -force
    dwi2response tournier $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
        $OUT_DIR/response_wm.txt -mask $OUT_DIR/mask.mif -force
    dwi2fod csd $OUT_DIR/dwi_flip_gradcheck_denoise_bias_anat.mif \
        $OUT_DIR/response_wm.txt $OUT_DIR/fod_wm.mif -mask $OUT_DIR/mask.mif -force
    mtnormalise $OUT_DIR/fod_wm.mif $OUT_DIR/fod_wm_norm.mif -mask $OUT_DIR/mask.mif -force
    
    if [ -f "$OUT_DIR/atlas/rois/roi_thalamus.nii.gz" ]; then
        tckgen $OUT_DIR/fod_wm_norm.mif $OUT_DIR/tracks_thalamus_cortex_raw.tck \
            -algorithm iFOD2 \
            -seed_image $OUT_DIR/atlas/rois/roi_thalamus.nii.gz \
            -include $OUT_DIR/atlas/rois/roi_isocortex.nii.gz \
            -mask $OUT_DIR/mask.mif \
            -select 2000 -cutoff 0.1 -minlength 5 -maxlength 25 -step 0.05 -angle 45 -force
        
        tcksift $OUT_DIR/tracks_thalamus_cortex_raw.tck $OUT_DIR/fod_wm_norm.mif \
            $OUT_DIR/tracks_thalamus_cortex_sift50.tck -term_number 50 -force
        echo "  ✓ Tractography complete"
    else
        echo "  ⚠ No ROIs - skipping tractography"
    fi
    
    #-------------------------------------------
    # STEP 10: QC REPORT
    #-------------------------------------------
    echo "[10/10] Generating QC report..."
    
    FA_MEAN=$(mrstats $OUT_DIR/fa.mif -mask $OUT_DIR/mask.mif -output mean 2>/dev/null || echo "NA")
    N_TRACKS=$(tckinfo $OUT_DIR/tracks_thalamus_cortex_sift50.tck 2>/dev/null | grep "count" | awk '{print $2}' || echo "NA")
    
    cat > $OUT_DIR/qc_report.txt << QCEOF
Subject: $SUBJECT
Session: $SESSION
Date: $(date)
FA Mean: $FA_MEAN
N Tracks: $N_TRACKS
Status: COMPLETE
QCEOF
    
    echo ""
    echo "========================================"
    echo "✓ COMPLETE: $SUBJECT"
    echo "  FA Mean: $FA_MEAN"
    echo "  N Tracks: $N_TRACKS"
    echo "  Log: $LOG_FILE"
    echo "========================================"
}

#===========================================
# MAIN
#===========================================
echo ""
echo "╔════════════════════════════════════════╗"
echo "║   MOUSE BRAIN DWI PIPELINE             ║"
echo "║   $(date)              ║"
echo "╚════════════════════════════════════════╝"
echo ""

# Get subjects to process
if [ $# -gt 0 ]; then
    SUBJECTS="$@"
    echo "Processing specified subjects: $SUBJECTS"
else
    SUBJECTS=$(ls -d $BIDS_DIR/sub-* 2>/dev/null | xargs -n1 basename)
    echo "Processing all subjects: $(echo $SUBJECTS | wc -w) found"
fi

# Process each subject
PROCESSED=0
FAILED=0

for SUBJECT in $SUBJECTS; do
    if [ -d "$BIDS_DIR/$SUBJECT" ]; then
        SESSION=$(ls $BIDS_DIR/$SUBJECT 2>/dev/null | head -1)
        
        if [ -n "$SESSION" ]; then
            if process_subject $SUBJECT $SESSION; then
                ((PROCESSED++))
            else
                ((FAILED++))
                echo "✗ FAILED: $SUBJECT"
            fi
        else
            echo "⚠ No session for $SUBJECT"
            ((FAILED++))
        fi
    else
        echo "⚠ Subject not found: $SUBJECT"
        ((FAILED++))
    fi
done

# Summary
echo ""
echo "╔════════════════════════════════════════╗"
echo "║   PIPELINE COMPLETE                    ║"
echo "╠════════════════════════════════════════╣"
echo "║   Processed: $PROCESSED                         ║"
echo "║   Failed: $FAILED                            ║"
echo "╚════════════════════════════════════════╝"

# Generate summary CSV
echo "subject,session,fa_mean,n_tracks,status" > $DERIV_DIR/pipeline_summary.csv
for SUBJ_PATH in $DERIV_DIR/sub-*; do
    SUBJECT=$(basename $SUBJ_PATH)
    SESSION=$(ls $SUBJ_PATH 2>/dev/null | head -1)
    QC=$SUBJ_PATH/$SESSION/dwi/qc_report.txt
    
    if [ -f "$QC" ]; then
        FA=$(grep "FA Mean" $QC | cut -d: -f2 | tr -d ' ')
        NT=$(grep "N Tracks" $QC | cut -d: -f2 | tr -d ' ')
        echo "$SUBJECT,$SESSION,$FA,$NT,complete" >> $DERIV_DIR/pipeline_summary.csv
    else
        echo "$SUBJECT,$SESSION,NA,NA,failed" >> $DERIV_DIR/pipeline_summary.csv
    fi
done

echo ""
echo "Summary: $DERIV_DIR/pipeline_summary.csv"
echo "Logs: $LOG_DIR/"
