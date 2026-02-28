#!/bin/bash
# Open mrview for all 4 APTR mice at the central AP graft slice
# Axis 0 = AP in flip1. mrview -voxel sets voxel coords as x,y,z = axis0,axis1,axis2
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

# Mouse 1: sub-AP72-862-901-NmF — central AP slice = 13
BASE="$DERIV/sub-AP72-862-901-NmF/ses-251210093342/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_anterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_central_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_posterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -voxel 13,50,53 \
  -mode 2 &

sleep 2

# Mouse 2: sub-AP72-862-901-R1M — central AP slice = 17
BASE="$DERIV/sub-AP72-862-901-R1M/ses-251210100043/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_anterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_central_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_posterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -voxel 17,50,53 \
  -mode 2 &

sleep 2

# Mouse 3: sub-Apallial17NmF — central AP slice = 15
BASE="$DERIV/sub-Apallial17NmF/ses-240816091633/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_anterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_central_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_posterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -voxel 15,50,60 \
  -mode 2 &

sleep 2

# Mouse 4: sub-Apallial20R1F — central AP slice = 14
BASE="$DERIV/sub-Apallial20R1F/ses-240909094232/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_anterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_central_dorsal_fod.mif" -overlay.opacity 0.4 \
  -overlay.load "$BASE/graft_segments_fod/seg9_graft_left_posterior_dorsal_fod.mif" -overlay.opacity 0.4 \
  -voxel 14,50,60 \
  -mode 2 &

echo "All 4 mrview windows opened at central AP graft slices"
echo "  NmF: slice 13 | R1M: slice 17 | 17NmF: slice 15 | 20R1F: slice 14"
