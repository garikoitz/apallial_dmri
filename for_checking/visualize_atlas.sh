#!/bin/bash

echo "========================================================================="
echo "VISUALIZATION 1: Your warped brain overlaid on the atlas"
echo "========================================================================="
echo "This shows how well your brain deforms to match the atlas shape"
echo ""

mrview /Users/glerma/soft/mbrain/local/atlas/allen_ccf/average_template_25.mif \
  -overlay.load /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/sub-596-C-OO-NmM_ses-240911103349_subject2atlas_Warped.nii.gz \
  -overlay.opacity 0.5 &

echo ""
echo "Press Enter to continue to the next visualization..."
read

echo ""
echo "========================================================================="
echo "VISUALIZATION 2: Warped ROIs overlaid on your original brain"
echo "========================================================================="
echo "Check if ROIs align with the correct anatomical structures"
echo "(cortex, striatum, thalamus)"
echo ""

mrview /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_t2w_highres.mif \
  -overlay.load /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/sub-596-C-OO-NmM_ses-240911103349_subject2atlas_InverseWarped.nii.gz \
  -overlay.opacity 0.5 \
  -roi.load /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/rois/sub-596-C-OO-NmM_ses-240911103349_cortex_left.mif \
  -roi.load /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/rois/sub-596-C-OO-NmM_ses-240911103349_striatum_left.mif \
  -roi.load /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/rois/sub-596-C-OO-NmM_ses-240911103349_thalamus_left.mif
