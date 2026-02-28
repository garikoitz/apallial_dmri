#!/bin/bash
# Create ROI edge masks for figure overlays
set -e
DERIV=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1

echo "=== CTRL: isocortex outline ==="
CTRL="$DERIV/sub-138-1010-L1F/ses-260115101441"
ROI_CTRL="$CTRL/atlas/rois/fod_space/roi_isocortex_left.nii.gz"
EDGE_CTRL="$CTRL/dwi/roi_isocortex_left_edge.mif"
maskfilter "$ROI_CTRL" erode /tmp/eroded_ctrl.mif -force
mrcalc "$ROI_CTRL" /tmp/eroded_ctrl.mif -subtract "$EDGE_CTRL" -force
echo "  Created $EDGE_CTRL"

echo "=== AP: graft outline ==="
AP="$DERIV/sub-769-C-O5-R1M/ses-240911115733"
ROI_AP="$AP/manual_roi/fod_space/roi_Graft_left_flip_fod.mif"
EDGE_AP="$AP/dwi/roi_Graft_left_edge.mif"
maskfilter "$ROI_AP" erode /tmp/eroded_ap.mif -force
mrcalc "$ROI_AP" /tmp/eroded_ap.mif -subtract "$EDGE_AP" -force
echo "  Created $EDGE_AP"

echo "=== APTR: graft outline ==="
APTR="$DERIV/sub-AP72-862-901-NmF/ses-251210093342"
ROI_APTR="$APTR/manual_roi/fod_space/roi_Graft_left_flip_fod.mif"
EDGE_APTR="$APTR/dwi/roi_Graft_left_edge.mif"
maskfilter "$ROI_APTR" erode /tmp/eroded_aptr.mif -force
mrcalc "$ROI_APTR" /tmp/eroded_aptr.mif -subtract "$EDGE_APTR" -force
echo "  Created $EDGE_APTR"

echo "=== Done ==="
