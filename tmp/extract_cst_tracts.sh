#!/usr/bin/env zsh
# Extract passthrough tracts from wholebrain using CST ROIs + same-hemisphere graft as seed
# CST ROIs are in T2w space → resample to FOD space first
set -e

ROI_SRC='/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/Shared drives/DTI Collab/mbrain_data/local/manual_roi_flip'
DERIV='/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1'

process_subject() {
    local SUBJ=$1 SES=$2 LEFT_ROI=$3 RIGHT_ROI=$4
    local BASE="$DERIV/$SUBJ/$SES"
    local DWI="$BASE/dwi"
    local FOD="$DWI/wm_fod_norm_group.mif"
    local WB_TCK="$DWI/wholebrain/wholebrain_1000k.tck"
    local OUT_DIR="$DWI/wholebrain/passthrough_cst"
    mkdir -p "$OUT_DIR"

    echo ""
    echo "============================================"
    echo "Subject: $SUBJ"
    echo "============================================"

    for HEMI in left right; do
        local ROI_FILE
        if [[ "$HEMI" == "left" ]]; then ROI_FILE="$LEFT_ROI"; else ROI_FILE="$RIGHT_ROI"; fi
        
        local ROI_T2W="$ROI_SRC/$ROI_FILE"
        local ROI_FOD="$OUT_DIR/roi_CST_${HEMI}_fod.mif"
        local GRAFT="$BASE/manual_roi/fod_space/roi_Graft_${HEMI}_flip_fod.mif"
        local OUT_TCK="$OUT_DIR/graft_${HEMI}_to_CST_${HEMI}.tck"

        echo ""
        echo "--- $HEMI hemisphere ---"

        # Step 1: Resample CST ROI from T2w to FOD space
        echo "  Resampling $ROI_FILE to FOD space..."
        micromamba run -n mbrain mrgrid "$ROI_T2W" regrid "$ROI_FOD" \
            -template "$FOD" -interp nearest -force -quiet
        
        # Verify resampled ROI has content
        local COUNT=$(micromamba run -n mbrain mrstats "$ROI_FOD" -output count -mask "$ROI_FOD" 2>/dev/null)
        echo "  Resampled ROI voxels: $COUNT"

        # Step 2: Check graft ROI exists
        if [[ ! -f "$GRAFT" ]]; then
            echo "  ERROR: Graft ROI not found: $GRAFT"
            continue
        fi

        # Step 3: Extract passthrough tract (include=graft AND include=CST)
        echo "  Extracting passthrough tract: graft_${HEMI} → CST_${HEMI}..."
        micromamba run -n mbrain tckedit "$WB_TCK" "$OUT_TCK" \
            -include "$GRAFT" -include "$ROI_FOD" -force -quiet
        
        # Count streamlines
        local N=$(micromamba run -n mbrain tckinfo "$OUT_TCK" -count 2>/dev/null | grep -o '[0-9]*' | tail -1)
        echo "  Result: $N streamlines in $OUT_TCK"
    done
}

# sub-769-C-O5-R1M (AP)
process_subject "sub-769-C-O5-R1M" "ses-240911115733" \
    "roi_CST_left_sub-769-C-O5-R1M_CST.mif" \
    "roi_CST_right_sub-769-C-O5-R1M_CST.mif"

# sub-Apallial20R1F (APTR)
process_subject "sub-Apallial20R1F" "ses-240909094232" \
    "roi_CST_left_sub-Apallial20R1F.mif" \
    "roi_CST_right_sub-Apallial20R1F.mif"

echo ""
echo "============================================"
echo "DONE - All 4 tracts extracted"
echo "============================================"
