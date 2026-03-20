#!/usr/bin/env python3
"""Calculate percentage of fixels by dominant direction (DV/AP/ML).

For each fixel in the bilateral ROI, the dominant direction is determined
by which component of the direction vector has the largest absolute value:
   dominant = argmax(|AP|, |DV|, |ML|)

Reported per mouse and per group, across all kept slices (1-11).

Usage:
    cd /home/glerma/glerma/soft/mbrain
    python3 tmp/diag_fixel_dominance.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors
)
from functions.config_loader import load_config

cfg = load_config('config1.yaml')
_data_base = Path(cfg['paths']['data_base'])
if not _data_base.exists() and str(_data_base).startswith('/export/home'):
    _data_base = Path(str(_data_base).replace('/export/home', '/home', 1))
DERIV = _data_base / 'data' / 'derivatives' / 'mrtrix_newmask'

SUBJECT_MAP = {
    'CTRL': [
        ('sub-129-921-R1F',   'ses-260115121231'),
        ('sub-129-921-R1L1F', 'ses-260116085433'),
        ('sub-129-1014-L1M',  'ses-260115113814'),
        ('sub-138-1010-L1F',  'ses-260115101441'),
    ],
    'AP': [
        ('sub-596-C-OO-NmM', 'ses-240911103349'),
        ('sub-769-C-O5-R1M', 'ses-240911115733'),
        ('sub-769-C-OO-NmM', 'ses-240911111853'),
        ('sub-875-C-OO-NM',  'ses-240911100253'),
    ],
    'APTR': [
        ('sub-AP72-862-901-NmF', 'ses-251210093342'),
        ('sub-AP72-862-901-R1M', 'ses-251210100043'),
        ('sub-Apallial17NmF',    'ses-240816091633'),
        ('sub-Apallial20R1F',    'ses-240909094232'),
    ],
}

COMP_NAMES = {0: 'AP', 1: 'DV', 2: 'ML'}

def get_roi_paths(subj_dir, group):
    if group == 'CTRL':
        base = subj_dir / 'atlas' / 'rois' / 'fod_space'
        return (base / 'roi_isocortex_left.nii.gz',
                base / 'roi_isocortex_right.nii.gz')
    else:
        base = subj_dir / 'manual_roi' / 'fod_space'
        return (base / 'roi_Graft_left_flip_fod.mif',
                base / 'roi_Graft_right_flip_fod.mif')


print('=' * 90)
print('Fixel Dominant Direction Analysis — Bilateral ROI')
print('  Method: For each fixel, dominant = argmax(|AP|, |DV|, |ML|)')
print('  ROI: bilateral (L|R), slices 1-11 (slice 0 dropped)')
print('=' * 90)

# Collect results for group-level summary
group_results = {}

for group in ['CTRL', 'AP', 'APTR']:
    print(f'\n{"="*90}')
    roi_label = 'isocortex' if group == 'CTRL' else 'graft'
    print(f'  {group} ({roi_label}, bilateral)')
    print(f'{"="*90}')

    group_counts = {'AP': 0, 'DV': 0, 'ML': 0}
    group_total = 0
    mouse_results = []

    for (subj, ses) in SUBJECT_MAP[group]:
        subj_dir = DERIV / subj / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        roi_left, roi_right = get_roi_paths(subj_dir, group)

        if not (fixel_dir / 'index.mif').exists():
            print(f'  ⚠ Skipping {subj}: no fixels')
            continue
        if not roi_left.exists() or not roi_right.exists():
            print(f'  ⚠ Skipping {subj}: ROI not found')
            continue

        fixel_data = load_fixel_dir(fixel_dir)
        mask_l = load_roi_mask(roi_left)
        mask_r = load_roi_mask(roi_right)
        roi_mask = mask_l | mask_r

        slice_vecs, _, _ = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)

        # Drop first slice, keep up to 11
        sorted_slices = sorted(slice_vecs.keys())
        kept_slices = sorted_slices[1:12]

        # Collect all fixel directions across kept slices
        all_dirs = []
        for s in kept_slices:
            vecs = slice_vecs[s]  # array of shape (n_fixels, 3)
            if len(vecs) > 0:
                all_dirs.append(np.array(vecs))

        if not all_dirs:
            print(f'  ⚠ {subj}: no fixels in kept slices')
            continue

        all_dirs = np.vstack(all_dirs)  # (N, 3) — columns: [AP, DV, ML]
        n_total = len(all_dirs)

        # Dominant direction = argmax of |component|
        dominant = np.argmax(np.abs(all_dirs), axis=1)  # 0=AP, 1=DV, 2=ML
        n_ap = int(np.sum(dominant == 0))
        n_dv = int(np.sum(dominant == 1))
        n_ml = int(np.sum(dominant == 2))

        pct_ap = 100.0 * n_ap / n_total
        pct_dv = 100.0 * n_dv / n_total
        pct_ml = 100.0 * n_ml / n_total

        short = subj.replace('sub-', '')
        print(f'\n  {short}:')
        print(f'    Total fixels: {n_total}')
        print(f'    AP-dominant:  {n_ap:6d}  ({pct_ap:5.1f}%)')
        print(f'    DV-dominant:  {n_dv:6d}  ({pct_dv:5.1f}%)')
        print(f'    ML-dominant:  {n_ml:6d}  ({pct_ml:5.1f}%)')

        group_counts['AP'] += n_ap
        group_counts['DV'] += n_dv
        group_counts['ML'] += n_ml
        group_total += n_total

        mouse_results.append({
            'subj': short, 'total': n_total,
            'AP': n_ap, 'DV': n_dv, 'ML': n_ml,
            'pct_AP': pct_ap, 'pct_DV': pct_dv, 'pct_ML': pct_ml,
        })

    # Group summary
    if group_total > 0:
        gpct_ap = 100.0 * group_counts['AP'] / group_total
        gpct_dv = 100.0 * group_counts['DV'] / group_total
        gpct_ml = 100.0 * group_counts['ML'] / group_total

        print(f'\n  --- {group} GROUP TOTAL ---')
        print(f'    Total fixels: {group_total}')
        print(f'    AP-dominant:  {group_counts["AP"]:6d}  ({gpct_ap:5.1f}%)')
        print(f'    DV-dominant:  {group_counts["DV"]:6d}  ({gpct_dv:5.1f}%)')
        print(f'    ML-dominant:  {group_counts["ML"]:6d}  ({gpct_ml:5.1f}%)')

        # Also show mean ± std across mice
        if mouse_results:
            pcts_ap = [m['pct_AP'] for m in mouse_results]
            pcts_dv = [m['pct_DV'] for m in mouse_results]
            pcts_ml = [m['pct_ML'] for m in mouse_results]
            print(f'\n    Mean ± SD across mice (n={len(mouse_results)}):')
            print(f'    AP: {np.mean(pcts_ap):5.1f}% ± {np.std(pcts_ap):4.1f}%')
            print(f'    DV: {np.mean(pcts_dv):5.1f}% ± {np.std(pcts_dv):4.1f}%')
            print(f'    ML: {np.mean(pcts_ml):5.1f}% ± {np.std(pcts_ml):4.1f}%')

    group_results[group] = {'counts': group_counts, 'total': group_total, 'mice': mouse_results}

# Final cross-group comparison table
print(f'\n\n{"="*90}')
print('CROSS-GROUP COMPARISON (mean ± SD across mice)')
print(f'{"="*90}')
print(f'{"Group":<8} {"n":>3} {"AP%":>12} {"DV%":>12} {"ML%":>12}')
print('-' * 50)
for group in ['CTRL', 'AP', 'APTR']:
    mice = group_results[group]['mice']
    n = len(mice)
    if n > 0:
        pcts_ap = [m['pct_AP'] for m in mice]
        pcts_dv = [m['pct_DV'] for m in mice]
        pcts_ml = [m['pct_ML'] for m in mice]
        print(f'{group:<8} {n:>3} {np.mean(pcts_ap):5.1f}±{np.std(pcts_ap):4.1f} '
              f'{np.mean(pcts_dv):5.1f}±{np.std(pcts_dv):4.1f} '
              f'{np.mean(pcts_ml):5.1f}±{np.std(pcts_ml):4.1f}')

print(f'\n✅ Done')
