#!/bin/bash
# Overnight passthrough tract regeneration for 4 renamed subjects
# Run: bash tmp/regen_passthrough_overnight.sh > tmp/regen_passthrough_log.txt 2>&1 &

set -e
cd /Users/glerma/soft/mbrain

SUBJECTS=(sub-769-C-O5-R1M sub-APtestmouseNmM sub-Apallial17NmF sub-Apallial20R1F)
LOG_PREFIX="$(date +%Y%m%d_%H%M%S)"

echo "=== Passthrough Tract Regeneration ==="
echo "Started: $(date)"
echo "Subjects: ${SUBJECTS[*]}"
echo ""

# Step 1: Extract passthrough (normal segments) - both hemispheres
echo "============================================"
echo "STEP 1: Extract passthrough (regular segments)"
echo "============================================"
for SUBJ in "${SUBJECTS[@]}"; do
  echo ""
  echo ">>> $SUBJ left hemisphere"
  micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
    --step extract --subject "$SUBJ" --hemi left

  echo ""
  echo ">>> $SUBJ right hemisphere"
  micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
    --step extract --subject "$SUBJ" --hemi right
done

echo ""
echo "============================================"
echo "STEP 2: Extract passthrough (shrink segments)"
echo "============================================"
for SUBJ in "${SUBJECTS[@]}"; do
  echo ""
  echo ">>> $SUBJ left hemisphere (shrink)"
  micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
    --step extract --shrink --subject "$SUBJ" --hemi left

  echo ""
  echo ">>> $SUBJ right hemisphere (shrink)"
  micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
    --step extract --shrink --subject "$SUBJ" --hemi right
done

echo ""
echo "============================================"
echo "STEP 3: MATLAB AFQ cleaning"
echo "============================================"
for SUBJ in "${SUBJECTS[@]}"; do
  for HEMI in left right; do
    echo ">>> $SUBJ $HEMI"
    micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
      --step clean --subject "$SUBJ" --hemi "$HEMI" || echo "WARNING: clean failed for $SUBJ $HEMI"
  done
done

echo ""
echo "============================================"
echo "STEP 4: SIFT2 weighting"
echo "============================================"
for SUBJ in "${SUBJECTS[@]}"; do
  for HEMI in left right; do
    echo ">>> $SUBJ $HEMI"
    micromamba run -n mbrain python functions/run_passthrough_pipeline.py \
      --step sift2 --subject "$SUBJ" --hemi "$HEMI" || echo "WARNING: sift2 failed for $SUBJ $HEMI"
  done
done

echo ""
echo "=== Done: $(date) ==="
