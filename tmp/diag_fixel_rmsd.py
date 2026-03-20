#!/usr/bin/env python3
"""Calculate RMSD of fixel direction vectors across mice within each group.

RMSD = sqrt( (1/(N-1)) * sum_i || v_i - v_bar ||^2 )

where v_i is the 3D direction vector [|AP|, |DV|, |ML|] for mouse i,
and v_bar is the group mean vector.

Two methods:
  A) K-means Dir 1: |scaled_directions[0]| per slice, averaged across slices
  B) Raw fixel mean: mean |direction| across all fixels per slice, averaged

All components use absolute values (orientation-invariant).
Reported per-slice (1-11) and overall per group.

Usage:
    cd /home/glerma/glerma/soft/mbrain
    python3 tmp/diag_fixel_rmsd.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans
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


def get_roi_paths(subj_dir, group):
    if group == 'CTRL':
        base = subj_dir / 'atlas' / 'rois' / 'fod_space'
        return (base / 'roi_isocortex_left.nii.gz',
                base / 'roi_isocortex_right.nii.gz')
    else:
        base = subj_dir / 'manual_roi' / 'fod_space'
        return (base / 'roi_Graft_left_flip_fod.mif',
                base / 'roi_Graft_right_flip_fod.mif')


def compute_rmsd(vectors):
    """Compute RMSD across N vectors (N x 3 array).

    RMSD = sqrt( 1/(N-1) * sum || v_i - v_bar ||^2 )
    """
    vectors = np.array(vectors)
    N = len(vectors)
    if N < 2:
        return 0.0
    v_bar = np.mean(vectors, axis=0)
    diffs = vectors - v_bar
    sq_norms = np.sum(diffs ** 2, axis=1)  # ||v_i - v_bar||^2
    rmsd = np.sqrt(np.sum(sq_norms) / (N - 1))
    return rmsd


# ============================================================
# Step 1: Load all data
# ============================================================
print('=' * 90)
print('Fixel Direction RMSD Analysis — Bilateral ROI')
print('  RMSD = sqrt( 1/(N-1) * sum || v_i - v_bar ||^2 )')
print('  Components use absolute values (orientation-invariant)')
print('  Methods: A) K-means Dir 1  |  B) Raw fixel mean direction')
print('=' * 90)

# Storage: group -> mouse -> slice -> vectors
# Method A: K-means Dir 1 per slice
# Method B: Raw fixel mean per slice
all_data = {}

for group in ['CTRL', 'AP', 'APTR']:
    roi_label = 'isocortex' if group == 'CTRL' else 'graft'
    print(f'\n--- {group} ({roi_label}, bilateral) ---')
    all_data[group] = {}

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

        print(f'  Loading {subj}...')
        fixel_data = load_fixel_dir(fixel_dir)
        mask_l = load_roi_mask(roi_left)
        mask_r = load_roi_mask(roi_right)
        roi_mask = mask_l | mask_r

        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)

        # Drop first slice, keep up to 11
        sorted_slices = sorted(slice_vecs.keys())
        kept_slices = sorted_slices[1:12]

        # K-means for Dir 1
        slice_km = slicewise_kmeans(
            {s: slice_vecs[s] for s in kept_slices},
            n_clusters=2,
            slice_afds={s: slice_peaks[s] for s in kept_slices},
            slice_roi_voxels={s: slice_voxels.get(s, len(slice_vecs[s])) for s in kept_slices},
        )

        short = subj.replace('sub-', '')
        mouse_data = {'kmeans_per_slice': {}, 'rawmean_per_slice': {}}

        for i, s in enumerate(kept_slices):
            # Method A: K-means Dir 1 (absolute values)
            if s in slice_km and 'scaled_directions' in slice_km[s]:
                sd = slice_km[s]['scaled_directions'][0]
                mouse_data['kmeans_per_slice'][i] = np.abs(sd)  # [|AP|, |DV|, |ML|]

            # Method B: Raw fixel mean direction (absolute values)
            vecs = np.array(slice_vecs[s])
            if len(vecs) > 0:
                mouse_data['rawmean_per_slice'][i] = np.mean(np.abs(vecs), axis=0)

        all_data[group][short] = mouse_data
        n_km = len(mouse_data['kmeans_per_slice'])
        n_raw = len(mouse_data['rawmean_per_slice'])
        print(f'    {n_km} slices (kmeans), {n_raw} slices (raw)')

# Determine common slice count
target_n_slices = 11

# ============================================================
# Step 2: Compute RMSD per slice and overall
# ============================================================
print(f'\n\n{"="*90}')
print('RMSD RESULTS')
print(f'{"="*90}')

for method, key in [('A) K-means Dir 1', 'kmeans_per_slice'),
                     ('B) Raw fixel mean', 'rawmean_per_slice')]:
    print(f'\n{"─"*90}')
    print(f'Method {method}')
    print(f'{"─"*90}')

    # Per-slice RMSD table
    print(f'\nPer-slice RMSD:')
    header = f'{"Slice":>6}'
    for group in ['CTRL', 'AP', 'APTR']:
        header += f'  {group:>10}'
    print(header)
    print('-' * 40)

    group_overall_vecs = {g: [] for g in ['CTRL', 'AP', 'APTR']}

    for sl in range(target_n_slices):
        row = f'{sl+1:>6}'
        for group in ['CTRL', 'AP', 'APTR']:
            mice = all_data[group]
            slice_vectors = []
            for mouse, mdata in mice.items():
                if sl in mdata[key]:
                    slice_vectors.append(mdata[key][sl])
            if len(slice_vectors) >= 2:
                rmsd = compute_rmsd(slice_vectors)
                row += f'  {rmsd:10.4f}'
            else:
                row += f'  {"N/A":>10}'
        print(row)

    # Overall RMSD: average each mouse's vectors across slices, then RMSD across mice
    print(f'\nOverall RMSD (mouse = mean across slices):')
    row = f'{"OVERALL":>6}'
    for group in ['CTRL', 'AP', 'APTR']:
        mice = all_data[group]
        mouse_avg_vecs = []
        for mouse, mdata in mice.items():
            vecs = list(mdata[key].values())
            if vecs:
                mouse_avg_vecs.append(np.mean(vecs, axis=0))
        if len(mouse_avg_vecs) >= 2:
            rmsd = compute_rmsd(mouse_avg_vecs)
            row += f'  {rmsd:10.4f}'
        else:
            row += f'  {"N/A":>10}'
    print(row)

    # Also print each mouse's average vector + norm for inspection
    print(f'\nMouse average vectors (|AP|, |DV|, |ML|) and ||v||:')
    for group in ['CTRL', 'AP', 'APTR']:
        mice = all_data[group]
        print(f'  {group}:')
        norms = []
        for mouse, mdata in mice.items():
            vecs = list(mdata[key].values())
            if vecs:
                avg = np.mean(vecs, axis=0)
                norm = np.linalg.norm(avg)
                norms.append(norm)
                print(f'    {mouse:<25} [{avg[0]:.4f}, {avg[1]:.4f}, {avg[2]:.4f}]  ||v||={norm:.4f}')
        if norms:
            mean_norm = np.mean(norms)
            std_norm = np.std(norms)
            print(f'    {"mean ||v||":<25} {mean_norm:.4f} ± {std_norm:.4f}')

    # For Method A: print RMSD contextualized by vector magnitude
    if key == 'kmeans_per_slice':
        print(f'\nRMSD relative to mean vector magnitude (RMSD / mean||v||):')
        header_ctx = f'{"":>6}'
        for group in ['CTRL', 'AP', 'APTR']:
            header_ctx += f'  {group:>10}'
        print(header_ctx)
        print('-' * 40)
        row_rmsd = f'{"RMSD":>6}'
        row_norm = f'{"||v||":>6}'
        row_ratio = f'{"ratio":>6}'
        for group in ['CTRL', 'AP', 'APTR']:
            mice = all_data[group]
            mouse_avg_vecs = []
            for mouse, mdata in mice.items():
                vecs = list(mdata[key].values())
                if vecs:
                    mouse_avg_vecs.append(np.mean(vecs, axis=0))
            if len(mouse_avg_vecs) >= 2:
                rmsd = compute_rmsd(mouse_avg_vecs)
                mean_norm = np.mean([np.linalg.norm(v) for v in mouse_avg_vecs])
                ratio = rmsd / mean_norm if mean_norm > 1e-8 else float('inf')
                row_rmsd += f'  {rmsd:10.4f}'
                row_norm += f'  {mean_norm:10.4f}'
                row_ratio += f'  {ratio:10.4f}'
            else:
                row_rmsd += f'  {"N/A":>10}'
                row_norm += f'  {"N/A":>10}'
                row_ratio += f'  {"N/A":>10}'
        print(row_rmsd)
        print(row_norm)
        print(row_ratio)
        print('  (ratio < 1 = RMSD smaller than vector size = consistent group)')

# ============================================================
# Step 3: CSV-friendly output
# ============================================================
print(f'\n\n{"="*90}')
print('CSV OUTPUT — Per-slice RMSD')
print(f'{"="*90}')

for method, key in [('KmeansDir1', 'kmeans_per_slice'),
                     ('RawFixelMean', 'rawmean_per_slice')]:
    print(f'\nMethod: {method}')
    print('Slice,CTRL_RMSD,AP_RMSD,APTR_RMSD')
    for sl in range(target_n_slices):
        vals = []
        for group in ['CTRL', 'AP', 'APTR']:
            mice = all_data[group]
            slice_vectors = []
            for mouse, mdata in mice.items():
                if sl in mdata[key]:
                    slice_vectors.append(mdata[key][sl])
            if len(slice_vectors) >= 2:
                vals.append(f'{compute_rmsd(slice_vectors):.4f}')
            else:
                vals.append('NA')
        print(f'{sl+1},{",".join(vals)}')

    # Overall row
    vals = []
    for group in ['CTRL', 'AP', 'APTR']:
        mice = all_data[group]
        mouse_avg_vecs = []
        for mouse, mdata in mice.items():
            vecs = list(mdata[key].values())
            if vecs:
                mouse_avg_vecs.append(np.mean(vecs, axis=0))
        if len(mouse_avg_vecs) >= 2:
            vals.append(f'{compute_rmsd(mouse_avg_vecs):.4f}')
        else:
            vals.append('NA')
    print(f'OVERALL,{",".join(vals)}')

print(f'\n✅ Done')
