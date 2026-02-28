#!/bin/bash
# Recover all passthrough tracts for sub-769-C-O5-R1M
# Run from: /Users/glerma/soft/mbrain
set -e

cd /Users/glerma/soft/mbrain/functions
LOG=/Users/glerma/soft/mbrain/tmp/recover_769.log

echo "=== Starting sub-769 recovery: $(date) ===" | tee "$LOG"

# Step 1: Non-shrink extraction (both hemispheres)
echo "" | tee -a "$LOG"
echo "=== STEP 1a: Extract passthrough tracts (left, non-shrink) ===" | tee -a "$LOG"
micromamba run -n mbrain python extract_passthrough_tracts.py --subject sub-769-C-O5-R1M --hemi left 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== STEP 1b: Extract passthrough tracts (right, non-shrink) ===" | tee -a "$LOG"
micromamba run -n mbrain python extract_passthrough_tracts.py --subject sub-769-C-O5-R1M --hemi right 2>&1 | tee -a "$LOG"

# Step 2: Shrink extraction (both hemispheres)
echo "" | tee -a "$LOG"
echo "=== STEP 2a: Extract passthrough tracts (left, shrink) ===" | tee -a "$LOG"
micromamba run -n mbrain python extract_passthrough_tracts.py --subject sub-769-C-O5-R1M --hemi left --shrink 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== STEP 2b: Extract passthrough tracts (right, shrink) ===" | tee -a "$LOG"
micromamba run -n mbrain python extract_passthrough_tracts.py --subject sub-769-C-O5-R1M --hemi right --shrink 2>&1 | tee -a "$LOG"

# Step 3: MATLAB cleaning (both directories, both hemispheres)
echo "" | tee -a "$LOG"
echo "=== STEP 3a: MATLAB cleaning (non-shrink, left) ===" | tee -a "$LOG"
micromamba run -n mbrain python run_passthrough_pipeline.py --subject sub-769-C-O5-R1M --hemi left --step clean 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== STEP 3b: MATLAB cleaning (non-shrink, right) ===" | tee -a "$LOG"
micromamba run -n mbrain python run_passthrough_pipeline.py --subject sub-769-C-O5-R1M --hemi right --step clean 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== STEP 3c: MATLAB cleaning (shrink, left) ===" | tee -a "$LOG"
micromamba run -n mbrain python run_passthrough_pipeline.py --subject sub-769-C-O5-R1M --hemi left --step clean --shrink 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== STEP 3d: MATLAB cleaning (shrink, right) ===" | tee -a "$LOG"
micromamba run -n mbrain python run_passthrough_pipeline.py --subject sub-769-C-O5-R1M --hemi right --step clean --shrink 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "=== RECOVERY COMPLETE: $(date) ===" | tee -a "$LOG"

# Final summary
echo "" | tee -a "$LOG"
B=/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1/sub-769-C-O5-R1M/ses-240911115733/dwi/wholebrain
echo "passthrough_tracts: $(ls "$B/passthrough_tracts/"*.tck 2>/dev/null | wc -l | tr -d ' ') tck files" | tee -a "$LOG"
echo "passthrough_tracts_shrink: $(ls "$B/passthrough_tracts_shrink/"*.tck 2>/dev/null | wc -l | tr -d ' ') tck files" | tee -a "$LOG"
