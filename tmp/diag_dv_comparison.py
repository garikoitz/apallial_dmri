#!/usr/bin/env python3
"""Diagnostic: compare ROI voxels, fixel counts, peaks, coherence between CTRL and APTR."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from functions.fixel_roi_analysis import load_fixel_dir, load_roi_mask, extract_slicewise_vectors, slicewise_kmeans
from functions.config_loader import load_config
from pathlib import Path

cfg = load_config('config1.yaml')
DERIV = Path(cfg['paths']['derivatives_flip1'])

subjects = {
    'APTR': [
        ('sub-AP72-862-901-NmF', 'ses-251210093342', 10, 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'),
        ('sub-Apallial17NmF',    'ses-240816091633', 10, 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'),
    ],
}

for group, mice in subjects.items():
    print(f'\n{"="*60}')
    print(f'  {group}')
    print(f'{"="*60}')
    for sid, ses, slc, roi_rel in mice:
        subj_dir = DERIV / sid / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        roi_path = subj_dir / roi_rel

        fd = load_fixel_dir(fixel_dir)
        roi = load_roi_mask(roi_path)
        sv, sp, svox = extract_slicewise_vectors(fd, roi, slice_axis=0)

        if slc not in sv:
            avail = sorted(sv.keys())
            nearest = min(avail, key=lambda s: abs(s - slc))
            print(f'  {sid} slice={slc} NOT AVAILABLE, using {nearest}')
            slc = nearest

        n_vox = svox.get(slc, 0)
        n_fix = len(sv.get(slc, []))
        peaks = sp.get(slc, np.array([]))

        km = slicewise_kmeans(
            {slc: sv[slc]}, n_clusters=2,
            slice_afds={slc: sp[slc]},
            slice_roi_voxels={slc: n_vox})
        ks = km[slc]

        # DV-dominant cluster
        dv_vals = [np.abs(sd[1]) for sd in ks['scaled_directions']]
        best_k = int(np.argmax(dv_vals))
        other_k = 1 - best_k

        # Compute peak stats per cluster
        labels_all = np.zeros(n_fix)  # not directly available, re-derive
        # Just use the cluster sizes and overall stats
        
        print(f'\n  {sid} slice={slc}:')
        print(f'    ROI voxels in slice:   {n_vox}')
        print(f'    Total fixels in slice: {n_fix}')
        print(f'    Fixels/voxel ratio:    {n_fix/max(n_vox,1):.2f}')
        print(f'    Mean peak (all fixels): {np.mean(peaks):.4f}')
        print(f'    Sum peaks (all fixels): {np.sum(peaks):.2f}')
        print(f'    --- DV-dom cluster (k={best_k}) ---')
        print(f'      n_fixels:    {ks["cluster_sizes"][best_k]}')
        print(f'      scaled_dir:  {ks["scaled_directions"][best_k]}')
        print(f'      |DV|:        {np.abs(ks["scaled_directions"][best_k][1]):.4f}')
        print(f'      coherence:   {ks["coherences"][best_k]:.4f}')
        print(f'      mean_peak:   {ks["mean_peaks"][best_k]:.4f}')
        print(f'    --- Other cluster (k={other_k}) ---')
        print(f'      n_fixels:    {ks["cluster_sizes"][other_k]}')
        print(f'      scaled_dir:  {ks["scaled_directions"][other_k]}')
        print(f'      coherence:   {ks["coherences"][other_k]:.4f}')
        print(f'      mean_peak:   {ks["mean_peaks"][other_k]:.4f}')
