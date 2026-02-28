#!/bin/bash
set -e
BASE="/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1/sub-Apallial20R1F/ses-240909094232/dwi/wholebrain"

echo "Subsampling wholebrain_1000k.tck -> wholebrain_5k.tck ..."
micromamba run -n mbrain tckedit \
    "$BASE/wholebrain_1000k.tck" \
    "$BASE/wholebrain_5k.tck" \
    -number 5000 \
    -force

echo "Done."
ls -lah "$BASE/wholebrain_5k.tck"
