#!/bin/bash
# Capture FOD and fixel screenshots for 3 groups – version 03
# Same slices as v02, but with ROI edge overlay (isocortex for CTRL, graft for AP/APTR)
# Color: yellow (#ffff00 = 1,1,0) for high visibility
set -e
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1
OUT=/Users/glerma/soft/mbrain/tmp

CTRL_BASE="$DERIV/sub-138-1010-L1F/ses-260115101441"
AP_BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733"
APTR_BASE="$DERIV/sub-AP72-862-901-NmF/ses-251210093342"

CTRL_EDGE="$CTRL_BASE/dwi/roi_isocortex_left_edge.mif"
AP_EDGE="$AP_BASE/dwi/roi_Graft_left_edge.mif"
APTR_EDGE="$APTR_BASE/dwi/roi_Graft_left_edge.mif"

# ============================================================
# ROW 1: FOD screenshots with ROI outline
# ============================================================

echo "=== CTRL FOD ==="
mrview "$CTRL_BASE/dwi/t2w_highres_flip.mif" \
  -odf.load_sh "$CTRL_BASE/dwi/wm_fod_norm_group.mif" \
  -overlay.load "$CTRL_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 24,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fod_v03_ \
  -capture.grab \
  -exit

echo "=== AP FOD ==="
mrview "$AP_BASE/dwi/t2w_highres_flip.mif" \
  -odf.load_sh "$AP_BASE/dwi/wm_fod_norm_group.mif" \
  -overlay.load "$AP_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 25,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fod_v03_ \
  -capture.grab \
  -exit

echo "=== APTR FOD ==="
mrview "$APTR_BASE/dwi/t2w_highres_flip.mif" \
  -odf.load_sh "$APTR_BASE/dwi/wm_fod_norm_group.mif" \
  -overlay.load "$APTR_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 24,50,54 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fod_v03_ \
  -capture.grab \
  -exit

# ============================================================
# ROW 2: Fixel screenshots with ROI outline
# ============================================================

echo "=== CTRL Fixel ==="
mrview "$CTRL_BASE/dwi/t2w_highres_flip.mif" \
  -fixel.load "$CTRL_BASE/dwi/fixels_group/afd.mif" \
  -overlay.load "$CTRL_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 24,50,53 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ctrl_fixel_v03_ \
  -capture.grab \
  -exit

echo "=== AP Fixel ==="
mrview "$AP_BASE/dwi/t2w_highres_flip.mif" \
  -fixel.load "$AP_BASE/dwi/fixels_group/afd.mif" \
  -overlay.load "$AP_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 25,49,60 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix ap_fixel_v03_ \
  -capture.grab \
  -exit

echo "=== APTR Fixel ==="
mrview "$APTR_BASE/dwi/t2w_highres_flip.mif" \
  -fixel.load "$APTR_BASE/dwi/fixels_group/afd.mif" \
  -overlay.load "$APTR_EDGE" -overlay.opacity 0.7 -overlay.colour 1,1,0 \
  -voxel 24,50,54 \
  -plane 0 \
  -noannotations \
  -capture.folder "$OUT" \
  -capture.prefix aptr_fixel_v03_ \
  -capture.grab \
  -exit

echo "=== Done! ==="
ls -lah "$OUT"/{ctrl,ap,aptr}_{fod,fixel}_v03_*
