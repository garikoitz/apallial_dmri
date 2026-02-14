#!/bin/bash
# Setup flip2 pipeline for the 6 missing AP/APTR subjects
# According to config.yaml - only applies flip 2 to old pipeline DWI (which already has permutation + flip1)

DERIV_OLD="/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc"
DERIV_NEW="/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip2"

# Subject list with sessions
declare -a SUBJECTS=(
  "sub-769-C-O5-R1M:ses-240911115733"
  "sub-769-C-OO-NmM:ses-240911111853"
  "sub-875-C-OO-NM:ses-240911100253"
  "sub-AP72-862-901-R1M:ses-251210100043"
  "sub-Apallial17NmF:ses-240816091633"
  "sub-Apallial20R1F:ses-240909094232"
)

for entry in "${SUBJECTS[@]}"; do
  sub="${entry%%:*}"
  ses="${entry##*:}"
  
  echo ""
  echo "=============================================="
  echo "Processing $sub/$ses"
  echo "=============================================="
  
  OLD="$DERIV_OLD/$sub/$ses"
  NEW="$DERIV_NEW/$sub/$ses"
  
  # Skip if already done
  if [ -f "$NEW/dwi/wm_fod_flip2.mif" ]; then
    echo "  SKIP: Already processed (FOD exists)"
    continue
  fi
  
  # Create folders
  mkdir -p "$NEW/dwi"
  mkdir -p "$NEW/manual_roi/flipped"
  mkdir -p "$NEW/dwi/graft_segments_v2"
  mkdir -p "$NEW/dwi/segment_tracts"
  
  # Find DWI in old pipeline
  OLD_DWI=$(ls "$OLD/dwi/dwi_flip_gradcheck_denoise_bias.mif" 2>/dev/null)
  if [ -z "$OLD_DWI" ]; then
    OLD_DWI=$(ls "$OLD/dwi/dwi_preproc"*.mif 2>/dev/null | head -1)
  fi
  
  if [ -z "$OLD_DWI" ]; then
    echo "  ERROR: No DWI found in old pipeline"
    continue
  fi
  echo "  DWI source: $(basename $OLD_DWI)"
  
  # Step 1: Apply flip2 to DWI
  echo "  Step 1: Applying flip2 to DWI..."
  mrtransform "$OLD_DWI" "$NEW/dwi/dwi_flip2.mif" -flip 2 -reorient_fod no -force
  
  # Step 2: Generate mask
  echo "  Step 2: Generating mask..."
  dwi2mask "$NEW/dwi/dwi_flip2.mif" "$NEW/dwi/mask_flip2.mif" -force
  
  # Step 3: Response functions
  echo "  Step 3: Estimating response functions..."
  dwi2response dhollander "$NEW/dwi/dwi_flip2.mif" \
    "$NEW/dwi/response_wm.txt" "$NEW/dwi/response_gm.txt" "$NEW/dwi/response_csf.txt" \
    -mask "$NEW/dwi/mask_flip2.mif" -force
  
  # Step 4: FOD
  echo "  Step 4: Computing FOD..."
  dwi2fod msmt_csd "$NEW/dwi/dwi_flip2.mif" \
    "$NEW/dwi/response_wm.txt" "$NEW/dwi/wm_fod_flip2.mif" \
    "$NEW/dwi/response_gm.txt" "$NEW/dwi/gm_fod_flip2.mif" \
    "$NEW/dwi/response_csf.txt" "$NEW/dwi/csf_fod_flip2.mif" \
    -mask "$NEW/dwi/mask_flip2.mif" -force
  
  # Step 5: Copy T2w
  echo "  Step 5: Copying T2w..."
  if [ -f "$OLD/dwi/t2w_highres_flip.mif" ]; then
    cp "$OLD/dwi/t2w_highres_flip.mif" "$NEW/dwi/"
  else
    echo "    WARNING: T2w not found"
  fi
  
  # Step 6: Copy manual ROIs (flipped)
  echo "  Step 6: Copying flipped manual ROIs..."
  if [ -d "$OLD/manual_roi/flipped" ]; then
    cp "$OLD/manual_roi/flipped"/*.mif "$NEW/manual_roi/flipped/" 2>/dev/null
  else
    echo "    WARNING: No flipped manual ROIs found"
  fi
  
  # Step 7: Copy graft segments
  echo "  Step 7: Copying graft segments..."
  if [ -d "$OLD/dwi/graft_segments_v2" ]; then
    cp "$OLD/dwi/graft_segments_v2"/*.nii.gz "$NEW/dwi/graft_segments_v2/" 2>/dev/null
  else
    echo "    WARNING: No graft_segments_v2 found"
  fi
  
  echo "  DONE: $sub"
done

echo ""
echo "=============================================="
echo "All subjects processed!"
echo "=============================================="
