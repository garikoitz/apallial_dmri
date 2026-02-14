#!/bin/bash
# Run MATLAB tract cleaning in background
# Usage: ./run_matlab_cleaning.sh

cd /Users/glerma/soft/mbrain

LOG_FILE="/tmp/matlab_cleaning_$(date +%Y%m%d_%H%M%S).log"

echo "Starting MATLAB tract cleaning..."
echo "Log file: $LOG_FILE"
echo "Started at: $(date)"
echo ""

# Run in background with nohup
/opt/homebrew/bin/micromamba run -n mbrain python clean_all_tracts_matlab.py >> "$LOG_FILE" 2>&1

echo ""
echo "Finished at: $(date)"
echo "Check log: $LOG_FILE"
