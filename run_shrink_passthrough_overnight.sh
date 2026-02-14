#!/bin/bash
# =============================================================================
# SHRINK PASSTHROUGH TRACTOGRAPHY - OVERNIGHT BATCH SCRIPT
# =============================================================================
# This script runs passthrough tract extraction and cleaning using SHRUNK graft 
# SEGMENTS (eroded by 1 voxel) for AP and APTR subjects.
#
# Steps:
# 1. Extract passthrough tracts (LEFT hemisphere) - using shrunk dorsal segments
# 2. MATLAB AFQ cleaning (LEFT hemisphere)  
# 3. Extract passthrough tracts (RIGHT hemisphere) - using shrunk dorsal segments
# 4. MATLAB AFQ cleaning (RIGHT hemisphere)
#
# Segment source: dwi/graft_segments_fod/shrink/seg9_graft_{hemi}_{segment}_fod_shrink.mif
# Output folder: dwi/wholebrain/passthrough_tracts_shrink/
#
# Usage:
#   nohup bash run_shrink_passthrough_overnight.sh > logs/shrink_passthrough_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#
# Created: 2026-02-04
# =============================================================================

set -e  # Exit on error

# Setup
cd /Users/glerma/soft/mbrain
mkdir -p logs

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="logs/shrink_passthrough_${TIMESTAMP}.log"

echo "============================================================"
echo "SHRINK PASSTHROUGH TRACTOGRAPHY - OVERNIGHT BATCH"
echo "============================================================"
echo "Started: $(date)"
echo "Log file: ${LOGFILE}"
echo "============================================================"

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate mbrain

# Function to run a step with timing
run_step() {
    local step_name="$1"
    local hemi="$2"
    local step_type="$3"
    
    echo ""
    echo "============================================================"
    echo "STEP: ${step_name} - ${hemi} hemisphere"
    echo "Started: $(date)"
    echo "============================================================"
    
    local start_time=$(date +%s)
    
    python functions/run_passthrough_pipeline.py \
        --shrink \
        --step "${step_type}" \
        --hemi "${hemi}"
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    echo ""
    echo "Completed: $(date)"
    echo "Duration: $((duration / 60)) minutes $((duration % 60)) seconds"
    echo "============================================================"
}

# Run all steps
echo ""
echo "PROCESSING PLAN:"
echo "  1. Extract LEFT hemisphere passthrough tracts (shrunk segments)"
echo "  2. Clean LEFT hemisphere tracts (MATLAB AFQ)"
echo "  3. Extract RIGHT hemisphere passthrough tracts (shrunk segments)"
echo "  4. Clean RIGHT hemisphere tracts (MATLAB AFQ)"
echo ""
echo "Subjects: 8 (4 AP + 4 APTR)"
echo "Tracts per subject: 12 (3 segments × 4 targets)"
echo "Total tracts: 192 (8 subjects × 12 tracts × 2 hemispheres)"
echo ""

# Step 1: Extract LEFT
run_step "EXTRACTION" "left" "extract"

# Step 2: Clean LEFT  
run_step "MATLAB CLEANING" "left" "clean"

# Step 3: Extract RIGHT
run_step "EXTRACTION" "right" "extract"

# Step 4: Clean RIGHT
run_step "MATLAB CLEANING" "right" "clean"

# Final summary
echo ""
echo "============================================================"
echo "ALL STEPS COMPLETE"
echo "Finished: $(date)"
echo "============================================================"
echo ""
echo "Output folders:"
echo "  dwi/wholebrain/passthrough_tracts_shrink/"
echo ""
echo "Next steps:"
echo "  - Review logs for any errors"
echo "  - Generate reports with cleaned tracts"
echo "============================================================"
