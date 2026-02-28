#!/bin/bash
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

echo "=== AP manual_roi/fod_space ==="
ls "$DERIV/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/" | grep -E "Graft|Thalamus"

echo ""
echo "=== APTR Apallial20R1F passthrough_tracts (Thalamus) ==="
ls "$DERIV/sub-Apallial20R1F/ses-240909094232/dwi/wholebrain/passthrough_tracts/" | grep "dorsal_to_Thalamus"
