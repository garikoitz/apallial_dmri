#!/bin/bash
# Open mrview for CTRL, AP and APTR to visually verify left/right hemisphere labeling
# Left ROI = RED (color 1), Right ROI = BLUE (color 2)
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

echo "Opening CTRL (sub-138-1010-L1F)..."
micromamba run -n mbrain mrview \
  "$DERIV/sub-138-1010-L1F/ses-260115101441/dwi/t2w_highres_flip.mif" \
  -overlay.load "$DERIV/sub-138-1010-L1F/ses-260115101441/atlas/rois/fod_space/roi_isocortex_left.nii.gz" \
  -overlay.opacity 0.6 \
  -overlay.load "$DERIV/sub-138-1010-L1F/ses-260115101441/atlas/rois/fod_space/roi_isocortex_right.nii.gz" \
  -overlay.opacity 0.6 \
  -plane 0 &

sleep 1

echo "Opening AP (sub-769-C-O5-R1M)..."
micromamba run -n mbrain mrview \
  "$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi/t2w_highres_flip.mif" \
  -overlay.load "$DERIV/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/roi_Graft_left_flip_fod.mif" \
  -overlay.opacity 0.6 \
  -overlay.load "$DERIV/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/roi_Graft_right_flip_fod.mif" \
  -overlay.opacity 0.6 \
  -plane 0 &

sleep 1

echo "Opening APTR (sub-Apallial20R1F)..."
micromamba run -n mbrain mrview \
  "$DERIV/sub-Apallial20R1F/ses-240909094232/dwi/t2w_highres_flip.mif" \
  -overlay.load "$DERIV/sub-Apallial20R1F/ses-240909094232/manual_roi/fod_space/roi_Graft_left_flip_fod.mif" \
  -overlay.opacity 0.6 \
  -overlay.load "$DERIV/sub-Apallial20R1F/ses-240909094232/manual_roi/fod_space/roi_Graft_right_flip_fod.mif" \
  -overlay.opacity 0.6 \
  -plane 0 &

echo "All 3 windows opened."
echo "RED = _left files, BLUE = _right files"
echo "Check if RED is on the anatomical LEFT and BLUE on the anatomical RIGHT."
