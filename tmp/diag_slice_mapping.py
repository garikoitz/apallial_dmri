#!/usr/bin/env python3
"""Diagnostic: map first_thal_fod (4x3 v09) to v10 3group slice index."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import load_fixel_dir, load_roi_mask, extract_slicewise_vectors
from functions.config_loader import load_config

cfg = load_config('config1.yaml')
_data_base = Path(cfg['paths']['data_base'])
if not _data_base.exists() and str(_data_base).startswith('/export/home'):
    _data_base = Path(str(_data_base).replace('/export/home', '/home', 1))
DERIV = _data_base / 'data' / 'derivatives' / 'mrtrix_newmask'

FIRST_THAL = {
    'CTRL': {
        'sub-129-921-R1F': 10, 'sub-129-921-R1L1F': 10,
        'sub-129-1014-L1M': 9, 'sub-138-1010-L1F': 9,
    },
    'AP': {
        'sub-596-C-OO-NmM': 13, 'sub-769-C-O5-R1M': 12,
        'sub-769-C-OO-NmM': 11, 'sub-875-C-OO-NM': 12,
    },
    'APTR': {
        'sub-AP72-862-901-NmF': 10, 'sub-AP72-862-901-R1M': 10,
        'sub-Apallial17NmF': 10, 'sub-Apallial20R1F': 11,
    },
}

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

def get_roi_paths(subj_dir, group):
    if group == 'CTRL':
        base = subj_dir / 'atlas' / 'rois' / 'fod_space'
        return (base / 'roi_isocortex_left.nii.gz', base / 'roi_isocortex_right.nii.gz')
    else:
        base = subj_dir / 'manual_roi' / 'fod_space'
        return (base / 'roi_Graft_left_flip_fod.mif', base / 'roi_Graft_right_flip_fod.mif')

hdr = "{:<8} {:<28} {:>14} {:<55} {:>12}".format(
    'Group', 'Subject', 'first_thal_fod', 'Sorted slices (all)', 'v10 slice#')
print(hdr)
print('-' * len(hdr))

for group in ['CTRL', 'AP', 'APTR']:
    for (subj, ses) in SUBJECT_MAP[group]:
        subj_dir = DERIV / subj / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        roi_left, roi_right = get_roi_paths(subj_dir, group)

        fixel_data = load_fixel_dir(fixel_dir)
        mask_l = load_roi_mask(roi_left)
        mask_r = load_roi_mask(roi_right)
        roi_mask = mask_l | mask_r

        slice_vecs, _, _ = extract_slicewise_vectors(fixel_data, roi_mask, slice_axis=0)
        sorted_slices = sorted(slice_vecs.keys())

        ftf = FIRST_THAL[group][subj]

        # v10: drop first slice, re-index. Displayed as slice 1-11.
        kept_slices = sorted_slices[1:12]  # drop first, take up to 11
        if ftf in kept_slices:
            v10_idx = str(kept_slices.index(ftf) + 1)  # 1-based display
        elif ftf in sorted_slices:
            v10_idx = "(dropped)"
        else:
            v10_idx = "NOT FOUND"

        row = "{:<8} {:<28} {:>14} {:<55} {:>12}".format(
            group, subj, ftf, str(sorted_slices), v10_idx)
        print(row)
    print()
