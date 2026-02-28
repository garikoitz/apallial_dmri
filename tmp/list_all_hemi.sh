#!/bin/bash
# List all left/right hemisphere files for ALL AP and APTR subjects
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

AP_SUBJECTS=(
    "sub-596-C-OO-NmM/ses-240905102010"
    "sub-769-C-O5-R1M/ses-240911115733"
    "sub-769-C-OO-NmM/ses-241009130101"
    "sub-875-C-OO-NM/ses-241104111337"
    "sub-APtestmouseNmM/ses-240911111853"
)

APTR_SUBJECTS=(
    "sub-AP72-862-901-NmF/ses-241112122519"
    "sub-AP72-862-901-R1M/ses-241104123730"
    "sub-Apallial17NmF/ses-240816091633"
    "sub-Apallial20R1F/ses-240909094232"
)

echo "======= AP subjects: manual_roi/fod_space ======="
for s in "${AP_SUBJECTS[@]}"; do
    echo "--- $s ---"
    ls "$DERIV/$s/manual_roi/fod_space/" 2>/dev/null | grep -E "_left|_right" | sort
done

echo ""
echo "======= APTR subjects: manual_roi/fod_space ======="
for s in "${APTR_SUBJECTS[@]}"; do
    echo "--- $s ---"
    ls "$DERIV/$s/manual_roi/fod_space/" 2>/dev/null | grep -E "_left|_right" | sort
done

echo ""
echo "======= APTR subjects: graft_segments_fod ======="
for s in "${APTR_SUBJECTS[@]}"; do
    echo "--- $s ---"
    ls "$DERIV/$s/dwi/graft_segments_fod/" 2>/dev/null | grep -E "_left|_right" | sort
done

echo ""
echo "======= AP subjects: graft_segments_fod ======="
for s in "${AP_SUBJECTS[@]}"; do
    echo "--- $s ---"
    ls "$DERIV/$s/dwi/graft_segments_fod/" 2>/dev/null | grep -E "_left|_right" | sort
done

echo ""
echo "======= APTR: passthrough_tracts (left/right in filenames) ======="
for s in "${APTR_SUBJECTS[@]}"; do
    echo "--- $s ---"
    ls "$DERIV/$s/dwi/wholebrain/passthrough_tracts/" 2>/dev/null | grep -E "_left|_right" | head -5
    echo "  ..."
done
