#!/usr/bin/env python3
"""Debug: Show both K-means clusters for sub-AP72-862-901-NmF at its first-thalamus slice."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from sklearn.cluster import KMeans
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors, get_directional_mean
)
from functions.config_loader import load_config
from pathlib import Path

cfg = load_config('config1.yaml')
DERIV = Path(cfg['paths']['derivatives_flip1'])

sid = 'sub-AP72-862-901-NmF'
ses = 'ses-251210093342'
subj_dir = DERIV / sid / ses
fixel_dir = subj_dir / 'dwi' / 'fixels_group'
roi_path = subj_dir / 'manual_roi' / 'fod_space' / 'roi_Graft_left_flip_fod.mif'

fixel_data = load_fixel_dir(fixel_dir)
roi_mask = load_roi_mask(roi_path)

slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
    fixel_data, roi_mask, slice_axis=0)

print(f"Available slices: {sorted(slice_vecs.keys())}")

# Check multiple slices around the first-thalamus
for s in sorted(slice_vecs.keys()):
    vecs = slice_vecs[s]
    peaks = slice_peaks[s]
    n_vox = slice_voxels.get(s, len(vecs))
    
    if len(vecs) < 4:
        print(f"\nSlice {s}: only {len(vecs)} fixels, skipping")
        continue
    
    abs_v = np.abs(vecs)
    km = KMeans(n_clusters=2, n_init=10, random_state=42)
    labels = km.fit_predict(abs_v)
    
    print(f"\n=== Slice {s}: {len(vecs)} fixels, {n_vox} ROI voxels ===")
    
    for k in range(2):
        mask_k = labels == k
        cv = vecs[mask_k]
        pk = peaks[mask_k]
        mean_dir, coh = get_directional_mean(cv)
        if mean_dir[2] < 0:
            mean_dir *= -1
        signs = np.sign(cv @ mean_dir)
        signs[signs == 0] = 1.0
        aligned = cv * signs[:, None]
        weighted_sum = np.sum(aligned * pk[:, None], axis=0)
        scaled_dir = weighted_sum / n_vox
        mean_peak = float(np.sum(pk)) / n_vox
        
        # d[0]=AP, d[1]=DV, d[2]=ML
        print(f"  Cluster {k}: n={len(cv)}, coherence={coh:.3f}")
        print(f"    mean_dir  = [{mean_dir[0]:+.4f}, {mean_dir[1]:+.4f}, {mean_dir[2]:+.4f}]  "
              f"(AP={abs(mean_dir[0]):.3f}, DV={abs(mean_dir[1]):.3f}, ML={abs(mean_dir[2]):.3f})")
        print(f"    scaled_dir= [{scaled_dir[0]:+.4f}, {scaled_dir[1]:+.4f}, {scaled_dir[2]:+.4f}]  "
              f"(AP={abs(scaled_dir[0]):.4f}, DV={abs(scaled_dir[1]):.4f}, ML={abs(scaled_dir[2]):.4f})")
        print(f"    mean_peak = {mean_peak:.4f}")
        print(f"    DV dominance (|scaled[1]|) = {abs(scaled_dir[1]):.4f}")
