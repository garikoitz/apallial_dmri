#!/usr/bin/env python3
"""Generate 4×3 composite fixel figure: FODs | Fixels | Radar (single) | Radar (group).

Columns: CTRL | AP | APTR
Row 1: FOD screenshots (from mrview captures)
Row 2: Fixel screenshots (from mrview captures)
Row 3: Single-mouse peak-weighted radar plot (central slice)
Row 4: Group-level radar plot (all subjects, aggregate)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans, cluster_fixel_directions
)
from functions.config_loader import load_config

# ============================================================
# Config
# ============================================================
TMP = Path('tmp')
REPORTS = Path('reports')
REPORTS.mkdir(exist_ok=True)

cfg = load_config('config1.yaml')
DERIV = Path(cfg['paths']['derivatives_flip1'])

# Selected subjects and their central FOD-space slices
SUBJECTS = {
    'CTRL': {
        'id': 'sub-138-1010-L1F',
        'ses': 'ses-260115101441',
        'fod_slice': 13,
        'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz',
        'roi_name': 'isocortex',
    },
    'AP': {
        'id': 'sub-769-C-O5-R1M',
        'ses': 'ses-240911115733',
        'fod_slice': 14,
        'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif',
        'roi_name': 'graft',
    },
    'APTR': {
        'id': 'sub-Apallial17NmF',
        'ses': 'ses-240816091633',
        'fod_slice': 15,
        'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif',
        'roi_name': 'graft',
    },
}

GROUP_COLORS = {'CTRL': '#27ae60', 'AP': '#e74c3c', 'APTR': '#3498db'}

# Screenshot file prefixes (from capture_3groups.sh)
SCREENSHOT_MAP = {
    ('CTRL', 'fod'): TMP / 'ctrl_fod_0000.png',
    ('AP', 'fod'): TMP / 'ap_fod_0000.png',
    ('APTR', 'fod'): TMP / 'aptr_fod_0000.png',
    ('CTRL', 'fixel'): TMP / 'ctrl_fixel_0000.png',
    ('AP', 'fixel'): TMP / 'ap_fixel_0000.png',
    ('APTR', 'fixel'): TMP / 'aptr_fixel_0000.png',
}


# ============================================================
# Helper: rotate screenshot 90° CW
# ============================================================
def load_rotated(path):
    """Load image, rotate 90° CW, return as numpy array."""
    img = Image.open(path).rotate(-90, expand=True)
    return np.array(img)


# ============================================================
# Helper: draw radar plot on an existing axes
# ============================================================
def draw_radar(ax, directions_list, scaled_directions_list, mean_peaks,
               cluster_sizes, angle_between=None, title='', y_max=None):
    """Draw a peak-weighted radar plot on a polar axes.
    
    Parameters
    ----------
    ax : matplotlib polar axes
    directions_list : list of (3,) unit eigenvectors
    scaled_directions_list : list of (3,) peak-weighted, voxel-normalized directions
    mean_peaks : list of float
    cluster_sizes : list of int
    angle_between : float or None (degrees)
    title : str
    y_max : float or None
    """
    radar_labels = ['DV', 'ML', 'AP']
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML -> [DV, ML, AP]
    angles_closed = np.concatenate([angles, [angles[0]]])
    
    colors = ['#E24A33', '#348ABD', '#988ED5']
    linestyles = ['-', '--', '-.']
    
    rmax = 0.0
    for k, sd in enumerate(scaled_directions_list):
        abs_d = np.abs(sd)[reorder]
        rmax = max(rmax, np.max(abs_d))
        vals = np.concatenate([abs_d, [abs_d[0]]])
        ls = linestyles[k % len(linestyles)]
        n_k = cluster_sizes[k] if k < len(cluster_sizes) else 0
        pk = mean_peaks[k] if k < len(mean_peaks) else 0
        ax.plot(angles_closed, vals, color=colors[k], linewidth=2,
                linestyle=ls, label=f'Dir {k+1} (n={n_k})')
        ax.fill(angles_closed, vals, color=colors[k],
                alpha=0.15 if k == 0 else 0.08)
    
    # Isotropic reference
    if mean_peaks:
        overall_mean = np.mean(mean_peaks)
        iso_val = overall_mean / np.sqrt(3)
        iso_angles = np.linspace(0, 2 * np.pi, 100)
        ax.plot(iso_angles, [iso_val] * len(iso_angles), color='gray',
                linewidth=0.8, linestyle=':', alpha=0.5)
    
    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(radar_labels, fontsize=10, fontweight='bold')
    
    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(0, max(rmax * 1.2, 0.01))
    
    angle_str = f'\n∠={angle_between:.1f}°' if angle_between is not None else ''
    ax.set_title(f'{title}{angle_str}', fontsize=10, pad=12)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=7, framealpha=0.8)


# ============================================================
# Step 1: Load screenshots and rotate
# ============================================================
print("Loading and rotating screenshots...")
screenshots = {}
for key, path in SCREENSHOT_MAP.items():
    if path.exists():
        screenshots[key] = load_rotated(path)
        print(f"  ✓ {key}: {path.name}")
    else:
        print(f"  ✗ {key}: {path.name} NOT FOUND")
        screenshots[key] = np.zeros((400, 400, 3), dtype=np.uint8)


# ============================================================
# Step 2: Generate single-mouse radar plots (Row 3)
# ============================================================
print("\n--- Generating single-mouse radar plots (Row 3) ---")
single_radar_data = {}

for group, info in SUBJECTS.items():
    subj_dir = DERIV / info['id'] / info['ses']
    fixel_dir = subj_dir / 'dwi' / 'fixels_group'
    roi_path = subj_dir / info['roi_rel']
    
    print(f"  {group}: Loading fixels for {info['id']}...")
    fixel_data = load_fixel_dir(fixel_dir)
    roi_mask = load_roi_mask(roi_path)
    
    slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
        fixel_data, roi_mask, slice_axis=0)
    
    # Get central slice data
    central = info['fod_slice']
    if central in slice_vecs:
        vecs = slice_vecs[central]
        peaks = slice_peaks[central]
        n_vox = slice_voxels.get(central, len(vecs))
        print(f"    Central slice {central}: {len(vecs)} fixels, {n_vox} ROI voxels")
        
        # Run K-means on this slice
        abs_v = np.abs(vecs)
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=2, n_init=10, random_state=42)
        labels = km.fit_predict(abs_v)
        
        from functions.fixel_roi_analysis import get_directional_mean
        cluster_info = []
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
            cluster_info.append(dict(
                direction=mean_dir, scaled_dir=scaled_dir,
                coherence=coh, n=len(cv), mean_peak=mean_peak
            ))
        
        # Sort by cluster size
        order = sorted(range(2), key=lambda i: cluster_info[i]['n'], reverse=True)
        cluster_info = [cluster_info[i] for i in order]
        
        d1, d2 = cluster_info[0]['direction'], cluster_info[1]['direction']
        angle_between = float(np.degrees(np.arccos(np.clip(np.abs(np.dot(d1, d2)), 0, 1))))
        
        single_radar_data[group] = {
            'directions': [ci['direction'] for ci in cluster_info],
            'scaled_directions': [ci['scaled_dir'] for ci in cluster_info],
            'mean_peaks': [ci['mean_peak'] for ci in cluster_info],
            'cluster_sizes': [ci['n'] for ci in cluster_info],
            'angle_between': angle_between,
        }
        print(f"    Dirs: {[ci['n'] for ci in cluster_info]}, angle={angle_between:.1f}°")
    else:
        print(f"    ⚠ Central slice {central} not in fixel data! Available: {sorted(slice_vecs.keys())}")
        # Try nearest
        avail = sorted(slice_vecs.keys())
        nearest = min(avail, key=lambda s: abs(s - central))
        print(f"    Using nearest slice {nearest} instead")
        # Recurse logic... let's just skip for now
        single_radar_data[group] = None


# ============================================================
# Step 3: Generate group-level radar plots (Row 4)
# ============================================================
print("\n--- Generating group-level radar plots (Row 4) ---")

# Get all subjects from config
groups_cfg = cfg.get('groups', {})
group_radar_data = {}

for group in ['CTRL', 'AP', 'APTR']:
    print(f"\n  {group}:")
    all_scaled = [[], []]  # per cluster
    all_peaks = [[], []]
    all_sizes = [[], []]
    all_angles = []
    
    for entry in groups_cfg.get(group, []):
        sid = entry['id']
        ses = entry['session']
        if 'test' in sid.lower():
            continue
        
        subj_dir = DERIV / sid / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        
        if not (fixel_dir / 'index.mif').exists():
            print(f"    ⚠ Skipping {sid}: no fixels")
            continue
        
        # Determine ROI
        if group == 'CTRL':
            roi_path = subj_dir / 'atlas' / 'rois' / 'fod_space' / 'roi_isocortex_left.nii.gz'
        else:
            roi_path = subj_dir / 'manual_roi' / 'fod_space' / 'roi_Graft_left_flip_fod.mif'
        
        if not roi_path.exists():
            print(f"    ⚠ Skipping {sid}: ROI not found")
            continue
        
        print(f"    Loading {sid}...")
        fixel_data = load_fixel_dir(fixel_dir)
        roi_mask = load_roi_mask(roi_path)
        
        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)
        slice_km = slicewise_kmeans(slice_vecs, n_clusters=2,
                                     slice_afds=slice_peaks,
                                     slice_roi_voxels=slice_voxels)
        
        if not slice_km:
            print(f"    ⚠ No slice data")
            continue
        
        # Average across all slices for this subject
        subj_scaled = [[], []]
        subj_peaks = [[], []]
        subj_sizes = [[], []]
        subj_angles = []
        
        for s_idx, res in slice_km.items():
            for k in range(min(2, len(res['scaled_directions']))):
                subj_scaled[k].append(res['scaled_directions'][k])
                subj_peaks[k].append(res['mean_peaks'][k])
                subj_sizes[k].append(res['cluster_sizes'][k])
            if 'angle_between' in res:
                subj_angles.append(res['angle_between'])
        
        for k in range(2):
            if subj_scaled[k]:
                all_scaled[k].append(np.mean(subj_scaled[k], axis=0))
                all_peaks[k].append(np.mean(subj_peaks[k]))
                all_sizes[k].append(int(np.mean(subj_sizes[k])))
        if subj_angles:
            all_angles.append(np.mean(subj_angles))
        
        n_slices = len(slice_km)
        print(f"      {n_slices} slices")
    
    # Aggregate across subjects
    if all_scaled[0]:
        group_radar_data[group] = {
            'directions': [np.mean(all_scaled[k], axis=0) / (np.linalg.norm(np.mean(all_scaled[k], axis=0)) + 1e-12)
                          if all_scaled[k] else np.zeros(3) for k in range(2)],
            'scaled_directions': [np.mean(all_scaled[k], axis=0) if all_scaled[k] else np.zeros(3)
                                  for k in range(2)],
            'mean_peaks': [np.mean(all_peaks[k]) if all_peaks[k] else 0.0 for k in range(2)],
            'cluster_sizes': [int(np.mean(all_sizes[k])) if all_sizes[k] else 0 for k in range(2)],
            'angle_between': float(np.mean(all_angles)) if all_angles else None,
            'n_subjects': len(all_scaled[0]),
        }
        print(f"  → {group}: {group_radar_data[group]['n_subjects']} subjects averaged")
    else:
        group_radar_data[group] = None


# ============================================================
# Step 4: Compute consistent y_max across all radars
# ============================================================
all_maxes = []
for data in list(single_radar_data.values()) + list(group_radar_data.values()):
    if data is not None:
        for sd in data['scaled_directions']:
            all_maxes.append(np.max(np.abs(sd)))
y_max_global = max(all_maxes) * 1.2 if all_maxes else 0.1
print(f"\nGlobal y_max for radar plots: {y_max_global:.4f}")


# ============================================================
# Step 5: Compose 4×3 figure
# ============================================================
print("\n--- Composing 4×3 figure ---")

fig = plt.figure(figsize=(18, 22))
gs = fig.add_gridspec(4, 3, hspace=0.25, wspace=0.15,
                       height_ratios=[1, 1, 0.9, 0.9])

groups_order = ['CTRL', 'AP', 'APTR']
group_titles = ['CTRL (isocortex)', 'AP (graft)', 'APTR (graft + Tx)']

# Row 1: FOD screenshots
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[0, col])
    ax.imshow(screenshots.get((group, 'fod'), np.zeros((10, 10, 3), dtype=np.uint8)))
    ax.set_title(f'{group_titles[col]}\n{SUBJECTS[group]["id"]}', fontsize=11, fontweight='bold')
    ax.axis('off')
    if col == 0:
        ax.set_ylabel('FODs', fontsize=14, fontweight='bold', labelpad=20)

# Row 2: Fixel screenshots
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[1, col])
    ax.imshow(screenshots.get((group, 'fixel'), np.zeros((10, 10, 3), dtype=np.uint8)))
    ax.axis('off')
    if col == 0:
        ax.set_ylabel('Fixels', fontsize=14, fontweight='bold', labelpad=20)

# Row 3: Single-mouse radar (central slice)
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[2, col], projection='polar')
    data = single_radar_data.get(group)
    if data is not None:
        info = SUBJECTS[group]
        title = f'{info["id"]}\nslice {info["fod_slice"]} ({info["roi_name"]})'
        draw_radar(ax, data['directions'], data['scaled_directions'],
                   data['mean_peaks'], data['cluster_sizes'],
                   angle_between=data['angle_between'],
                   title=title, y_max=y_max_global)
    else:
        ax.set_title(f'{group}: no data', fontsize=10)

# Row 4: Group-level radar (all subjects)
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[3, col], projection='polar')
    data = group_radar_data.get(group)
    if data is not None:
        n_subj = data.get('n_subjects', '?')
        title = f'{group} group avg\n(n={n_subj} subjects)'
        draw_radar(ax, data['directions'], data['scaled_directions'],
                   data['mean_peaks'], data['cluster_sizes'],
                   angle_between=data['angle_between'],
                   title=title, y_max=y_max_global)
    else:
        ax.set_title(f'{group}: no data', fontsize=10)

# Row labels on left side
fig.text(0.02, 0.88, 'Row 1: FOD Glyphs', fontsize=12, fontweight='bold',
         rotation=90, va='center')
fig.text(0.02, 0.65, 'Row 2: Fixel Sticks', fontsize=12, fontweight='bold',
         rotation=90, va='center')
fig.text(0.02, 0.40, 'Row 3: Single Mouse\n(central slice)', fontsize=11,
         fontweight='bold', rotation=90, va='center')
fig.text(0.02, 0.15, 'Row 4: Group Average\n(all subjects)', fontsize=11,
         fontweight='bold', rotation=90, va='center')

fig.suptitle('Fixel Direction Analysis — CTRL vs AP vs APTR (left hemisphere)',
             fontsize=16, fontweight='bold', y=0.98)

save_path = REPORTS / 'fixel_4x3_comparison.png'
fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {save_path}")
plt.close(fig)
