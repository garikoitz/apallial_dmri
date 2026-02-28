#!/bin/bash
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

echo "=== CTRL (sub-138-1010-L1F) atlas/rois/fod_space ==="
ls "$DERIV/sub-138-1010-L1F/ses-260115101441/atlas/rois/fod_space/" 2>/dev/null | grep -E "left|right" | sort

echo ""
echo "=== AP (sub-769-C-O5-R1M) manual_roi/fod_space ==="
ls "$DERIV/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/" 2>/dev/null | grep -E "left|right" | sort

echo ""
echo "=== APTR (sub-Apallial20R1F) manual_roi/fod_space ==="
ls "$DERIV/sub-Apallial20R1F/ses-240909094232/manual_roi/fod_space/" 2>/dev/null | grep -E "left|right" | sort

echo ""
echo "=== APTR graft_segments_fod ==="
ls "$DERIV/sub-Apallial20R1F/ses-240909094232/dwi/graft_segments_fod/" 2>/dev/null | grep -E "left|right" | head -10 | sort
