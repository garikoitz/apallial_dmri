#!/bin/bash
# Capture FOD and fixel screenshots for 3 groups – version 02
# T2w voxel coordinates (first thalamus slice per user):
#   CTRL: (24, 50, 53)  -- T2W slice 24 → FOD slice 14
#   AP:   (25, 49, 60)  -- T2W slice 25 → FOD slice 15
#   APTR: (24, 50, 54)  -- T2W slice 24 → FOD slice 14  (sub-AP72-862-901-NmF)
set -e
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1
OUT=/Users/glerma/soft/mbrain/tmp

# ============================================================
# ROW 1: FOD screenshots
# ============================================================

echo "=== CTRL FOD ==="
BASE="$DERIV/sub-138-1010-L1F/ses-260115101441/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 24,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fod_v02_ \
  -capture.grab \
  -exit

echo "=== AP FOD ==="
BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 25,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fod_v02_ \
  -capture.grab \
  -exit

echo "=== APTR FOD ==="
BASE="$DERIV/sub-AP72-862-901-NmF/ses-251210093342/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -odf.load_sh "$BASE/wm_fod_norm_group.mif" \
  -voxel 24,50,54 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fod_v02_ \
  -capture.grab \
  -exit

# ============================================================
# ROW 2: Fixel screenshots
# ============================================================

echo "=== CTRL Fixel ==="
BASE="$DERIV/sub-138-1010-L1F/ses-260115101441/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 24,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fixel_v02_ \
  -capture.grab \
  -exit

echo "=== AP Fixel ==="
BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 25,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fixel_v02_ \
  -capture.grab \
  -exit

echo "=== APTR Fixel ==="
BASE="$DERIV/sub-AP72-862-901-NmF/ses-251210093342/dwi"
micromamba run -n mbrain mrview "$BASE/t2w_highres_flip.mif" \
  -fixel.load "$BASE/fixels_group/afd.mif" \
  -voxel 24,50,54 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fixel_v02_ \
  -capture.grab \
  -exit

echo "=== Done! ==="
ls -lah "$OUT"/{ctrl,ap,aptr}_{fod,fixel}_v02_*
