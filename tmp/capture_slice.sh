#!/bin/bash
# Capture screenshot from mrview at central AP slice, then rotate 90° CW
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1
BASE="$DERIV/sub-Apallial20R1F/ses-240909094232"
OUT=/Users/glerma/soft/mbrain/tmp

micromamba run -n mbrain mrview "$BASE/dwi/t2w_highres_flip.mif" \
  -overlay.load "$BASE/manual_roi/fod_space/roi_Graft_left_flip_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/manual_roi/fod_space/roi_Graft_right_flip_fod.mif" -overlay.opacity 0.4 \
  -voxel 14,50,60 \
  -plane 0 \
  -capture.prefix "$OUT/mrview_cap" \
  -capture.grab \
  -exit
