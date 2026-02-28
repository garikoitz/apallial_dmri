#!/usr/bin/env zsh
# Run tckgen directly from graft→CST with RELAXED parameters
# Standard params: cutoff=0.1, angle=60, minlen=2, maxlen=15
# Relaxed params:  cutoff=0.05, angle=90, minlen=0.5, maxlen=30, no mask
set -e

DERIV='/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1'

run_tckgen() {
    local SUBJ=$1 SES=$2
    local BASE="$DERIV/$SUBJ/$SES"
    local DWI="$BASE/dwi"
    local FOD_NOMASK="$DWI/wm_fod_group_nomask.mif"
    local FOD_MASKED="$DWI/wm_fod_norm_group.mif"
    local OUT_DIR="$DWI/wholebrain/passthrough_cst"

    echo ""
    echo "============================================"
    echo "Subject: $SUBJ"
    echo "============================================"

    for HEMI in left right; do
        local GRAFT="$BASE/manual_roi/fod_space/roi_Graft_${HEMI}_flip_fod.mif"
        local CST="$OUT_DIR/roi_CST_${HEMI}_fod.mif"

        if [[ ! -f "$GRAFT" ]]; then echo "  SKIP $HEMI: no graft ROI"; continue; fi
        if [[ ! -f "$CST" ]]; then echo "  SKIP $HEMI: no CST ROI"; continue; fi

        echo ""
        echo "--- $HEMI hemisphere ---"

        # --- Attempt 1: nomask FOD, relaxed params ---
        local OUT1="$OUT_DIR/tckgen_graft_${HEMI}_to_CST_${HEMI}_relaxed_nomask.tck"
        echo "  [1] tckgen nomask FOD, cutoff=0.05, angle=90, minlen=0.5, maxlen=30, seeds=10k, select=1000"
        micromamba run -n mbrain tckgen "$FOD_NOMASK" "$OUT1" \
            -algorithm iFOD2 \
            -seed_image "$GRAFT" \
            -include "$CST" \
            -cutoff 0.05 \
            -angle 90 \
            -step 0.02 \
            -minlength 0.5 \
            -maxlength 30 \
            -seeds 10000 \
            -select 1000 \
            -force -quiet 2>&1 || true
        local N1=$(micromamba run -n mbrain tckinfo "$OUT1" 2>/dev/null | grep "count:" | head -1 | awk '{print $2}')
        echo "  → $N1 streamlines"

        # --- Attempt 2: masked FOD, relaxed params ---
        local OUT2="$OUT_DIR/tckgen_graft_${HEMI}_to_CST_${HEMI}_relaxed_masked.tck"
        echo "  [2] tckgen masked FOD, cutoff=0.05, angle=90, minlen=0.5, maxlen=30, seeds=10k, select=1000"
        micromamba run -n mbrain tckgen "$FOD_MASKED" "$OUT2" \
            -algorithm iFOD2 \
            -seed_image "$GRAFT" \
            -include "$CST" \
            -cutoff 0.05 \
            -angle 90 \
            -step 0.02 \
            -minlength 0.5 \
            -maxlength 30 \
            -seeds 10000 \
            -select 1000 \
            -force -quiet 2>&1 || true
        local N2=$(micromamba run -n mbrain tckinfo "$OUT2" 2>/dev/null | grep "count:" | head -1 | awk '{print $2}')
        echo "  → $N2 streamlines"

        # --- Attempt 3: nomask FOD, even more relaxed (cutoff=0.01) ---
        local OUT3="$OUT_DIR/tckgen_graft_${HEMI}_to_CST_${HEMI}_vrelaxed_nomask.tck"
        echo "  [3] tckgen nomask FOD, cutoff=0.01, angle=90, minlen=0.5, maxlen=30, seeds=10k, select=1000"
        micromamba run -n mbrain tckgen "$FOD_NOMASK" "$OUT3" \
            -algorithm iFOD2 \
            -seed_image "$GRAFT" \
            -include "$CST" \
            -cutoff 0.01 \
            -angle 90 \
            -step 0.02 \
            -minlength 0.5 \
            -maxlength 30 \
            -seeds 10000 \
            -select 1000 \
            -force -quiet 2>&1 || true
        local N3=$(micromamba run -n mbrain tckinfo "$OUT3" 2>/dev/null | grep "count:" | head -1 | awk '{print $2}')
        echo "  → $N3 streamlines"
    done
}

# Process both subjects
run_tckgen "sub-769-C-O5-R1M" "ses-240911115733"
run_tckgen "sub-Apallial20R1F" "ses-240909094232"

echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
for SUBJ_SES in "sub-769-C-O5-R1M/ses-240911115733" "sub-Apallial20R1F/ses-240909094232"; do
    echo "--- ${SUBJ_SES%%/*} ---"
    for f in "$DERIV/$SUBJ_SES/dwi/wholebrain/passthrough_cst/tckgen_"*.tck; do
        [[ -f "$f" ]] || continue
        local N=$(micromamba run -n mbrain tckinfo "$f" 2>/dev/null | grep "count:" | head -1 | awk '{print $2}')
        echo "  $(basename $f): $N streamlines"
    done
done
