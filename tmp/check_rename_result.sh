#!/bin/bash
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

echo "=== sub-769-C-O5-R1M manual_roi/fod_space ALL files ==="
ls "$DERIV/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/" | sort

echo ""
echo "=== sub-769-C-O5-R1M graft_segments_fod (first 6) ==="
ls "$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi/graft_segments_fod/" | sort | head -6

echo ""
echo "=== sub-769-C-O5-R1M passthrough_tracts (Thalamus only) ==="
ls "$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi/wholebrain/passthrough_tracts/" | grep "Thalamus" | sort
