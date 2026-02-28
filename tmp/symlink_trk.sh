#!/bin/bash
TRK_DIR="/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/Shared drives/DTI Collab/mbrain_data/data/derivatives/mrtrix_preproc_flip1/sub-Apallial20R1F/ses-240909094232/dwi/wholebrain/passthrough_tracts/trk"
LOCAL="/Users/glerma/soft/mbrain/tmp/freeview_ap20"
mkdir -p "$LOCAL"
# Symlink everything
for f in "$TRK_DIR"/*.trk "$TRK_DIR/t2w_highres_flip.nii.gz"; do
    ln -sf "$f" "$LOCAL/"
done
echo "=== Contents of $LOCAL ==="
ls "$LOCAL/"
