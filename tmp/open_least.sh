#!/bin/bash
# Open mrview for LEAST connected APTR mouse: sub-AP72-862-901-R1M
# NON-SHRINK tracts + destination ROIs (T2w space)
# Tracts grouped by destination (all 3 segments per target).
# Take 3 screenshots, toggling visibility per destination group:
#   reports/mrview_least_paleocortex.png
#   reports/mrview_least_striatum.png
#   reports/mrview_least_thalamus.png
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1
BASE="$DERIV/sub-AP72-862-901-R1M/ses-251210100043/dwi"
SEG="$BASE/graft_segments_fod"
TCK="$BASE/wholebrain/passthrough_tracts"
ROI="$DERIV/sub-AP72-862-901-R1M/ses-251210100043/manual_roi"

micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -overlay.load "$SEG/seg9_graft_left_anterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$SEG/seg9_graft_left_central_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$SEG/seg9_graft_left_posterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$ROI/roi_Palecortex_left_flip.mif" -overlay.opacity 0.3 \
  -overlay.load "$ROI/roi_Striatum_left_flip.mif" -overlay.opacity 0.3 \
  -overlay.load "$ROI/roi_Thalamus_left_flip.mif" -overlay.opacity 0.3 \
  -tractography.load "$TCK/seg9_graft_left_anterior_dorsal_to_Palecortex_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_central_dorsal_to_Palecortex_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_posterior_dorsal_to_Palecortex_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_anterior_dorsal_to_Striatum_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_central_dorsal_to_Striatum_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_posterior_dorsal_to_Striatum_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_anterior_dorsal_to_Thalamus_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_central_dorsal_to_Thalamus_left_cleaned.tck" \
  -tractography.load "$TCK/seg9_graft_left_posterior_dorsal_to_Thalamus_left_cleaned.tck" &
