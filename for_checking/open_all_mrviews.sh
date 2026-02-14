#!/bin/bash
# Open mrview for all 12 subjects with colored ROIs and tracts
# Color scheme:
#   Green (0,255,0): Graft/Isocortex seed segment
#   Cyan (0,255,255): Thalamus ROI
#   Magenta (255,0,255): Striatum ROI
#   Blue (0,150,255): Thalamus tracts
#   Orange (255,100,0): Striatum tracts

DERIV="/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2"

# CTRL subjects
echo "Opening CTRL subjects..."

# sub-129-1014-L1M
BASE="$DERIV/sub-129-1014-L1M/ses-260115113814"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/isocortex_segments_v2/isocortex_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-129-921-R1F
BASE="$DERIV/sub-129-921-R1F/ses-260115121231"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/isocortex_segments_v2/isocortex_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-129-921-R1L1F
BASE="$DERIV/sub-129-921-R1L1F/ses-260116085433"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/isocortex_segments_v2/isocortex_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_thalamus_left.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/fod_space/roi_striatum_left.nii.gz -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_isocortex_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-138-1010-L1F (has empty ROIs - no tracts)
BASE="$DERIV/sub-138-1010-L1F/ses-260115101441"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
[ -f "$BASE/dwi/isocortex_segments_v2/isocortex_left_anterior_dorsal.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/isocortex_segments_v2/isocortex_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/dwi/atlas/rois/roi_thalamus_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/roi_thalamus_left.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/dwi/atlas/rois/roi_striatum_left.nii.gz" ] && cmd="$cmd -overlay.load $BASE/dwi/atlas/rois/roi_striatum_left.nii.gz -overlay.opacity 0.5 -overlay.colour 255,0,255"
eval "$cmd &"
sleep 0.5

echo "Opening AP subjects..."

# sub-596-C-OO-NmM
BASE="$DERIV/sub-596-C-OO-NmM/ses-240911103349"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-769-C-O5-R1M
BASE="$DERIV/sub-769-C-O5-R1M/ses-240911115733"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-769-C-OO-NmM
BASE="$DERIV/sub-769-C-OO-NmM/ses-240911111853"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-875-C-OO-NM
BASE="$DERIV/sub-875-C-OO-NM/ses-240911100253"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

echo "Opening APTR subjects..."

# sub-AP72-862-901-NmF
BASE="$DERIV/sub-AP72-862-901-NmF/ses-251210093342"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-AP72-862-901-R1M
BASE="$DERIV/sub-AP72-862-901-R1M/ses-251210100043"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-Apallial17NmF
BASE="$DERIV/sub-Apallial17NmF/ses-240816091633"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"
sleep 0.5

# sub-Apallial20R1F
BASE="$DERIV/sub-Apallial20R1F/ses-240909094232"
cmd="mrview $BASE/dwi/t2w_highres_flip.mif"
cmd="$cmd -overlay.load $BASE/dwi/graft_segments_v2/graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.5 -overlay.colour 0,255,0"
[ -f "$BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Thalamus_left_flip.mif -overlay.opacity 0.5 -overlay.colour 0,255,255"
[ -f "$BASE/manual_roi/flipped/roi_Striatum_left_flip.mif" ] && cmd="$cmd -overlay.load $BASE/manual_roi/flipped/roi_Striatum_left_flip.mif -overlay.opacity 0.5 -overlay.colour 255,0,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_thalamus_left_cleaned.tck -tractography.colour 0,150,255"
[ -f "$BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck" ] && cmd="$cmd -tractography.load $BASE/dwi/segment_tracts/tracks_graft_left_anterior_dorsal_to_striatum_left_cleaned.tck -tractography.colour 255,100,0"
eval "$cmd &"

echo "Opened all 12 subjects"
