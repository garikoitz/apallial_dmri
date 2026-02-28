#!/bin/bash
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

for s in \
    "sub-769-C-O5-R1M/ses-240911115733" \
    "sub-APtestmouseNmM/ses-240911111853" \
    "sub-Apallial17NmF/ses-240816091633" \
    "sub-Apallial20R1F/ses-240909094232"
do
    echo "=== $s ==="
    for d in \
        "manual_roi/fod_space" \
        "manual_roi/fod_space/shrink" \
        "dwi/graft_segments_fod" \
        "dwi/graft_segments_fod/shrink" \
        "dwi/wholebrain/passthrough_tracts" \
        "dwi/wholebrain/passthrough_tracts_shrink"
    do
        L=$(ls "$DERIV/$s/$d/" 2>/dev/null | grep -c "_left" || echo 0)
        R=$(ls "$DERIV/$s/$d/" 2>/dev/null | grep -c "_right" || echo 0)
        TOTAL=$(ls "$DERIV/$s/$d/" 2>/dev/null | grep -cE "_left|_right" || echo 0)
        echo "  $d: left=$L right=$R total=$TOTAL"
    done
done
