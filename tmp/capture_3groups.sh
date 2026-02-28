#!/bin/bash
# Capture FOD and fixel screenshots for 3 groups at central AP slices
# T2w voxel coordinates (converted from FOD central slices):
#   CTRL: (14, 50, 53)  -- FOD slice 13
#   AP:   (14, 49, 60)  -- FOD slice 14
#   APTR: (14, 50, 60)  -- FOD slice 15
set -e
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1
OUT=/Users/glerma/soft/mbrain/tmp

# ============================================================
# ROW 1: FOD screenshots
# ============================================================

echo "=== CTRL FOD ==="
BASE="$DERIV/sub-138-1010-L1F/ses-260115101441/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 14,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fod_ \
  -capture.grab \
  -exit

echo "=== AP FOD ==="
BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 14,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fod_ \
  -capture.grab \
  -exit

echo "=== APTR FOD ==="
BASE="$DERIV/sub-Apallial17NmF/ses-240816091633/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 14,50,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fod_ \
  -capture.grab \
  -exit

# ============================================================
# ROW 2: Fixel screenshots
# ============================================================

echo "=== CTRL Fixel ==="
BASE="$DERIV/sub-138-1010-L1F/ses-260115101441/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 14,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fixel_ \
  -capture.grab \
  -exit

echo "=== AP Fixel ==="
BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 14,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fixel_ \
  -capture.grab \
  -exit

echo "=== APTR Fixel ==="
BASE="$DERIV/sub-Apallial17NmF/ses-240816091633/dwi"
mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 14,50,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fixel_ \
  -capture.grab \
  -exit

echo "=== Done! ==="
ls -lah "$OUT"/{ctrl,ap,aptr}_{fod,fixel}_*
