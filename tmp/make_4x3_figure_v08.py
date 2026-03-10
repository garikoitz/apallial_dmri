#!/usr/bin/env python3
"""Generate 3×3 composite fixel figure v08: FODs | Fixels | 4-mice radar.

Changes from v07:
  - Uses mrtrix_newmask pipeline (new mask derivatives) instead of mrtrix_preproc_flip1.
  - Same layout: Row1=FOD, Row2=Fixel screenshots, Row3=per-mouse radar.
  - first_thal_fod slices kept from v07 with automatic nearest-slice fallback.

Columns: CTRL | AP | APTR
Row 1: FOD screenshots (from mrview captures, same as v02)
Row 2: Fixel screenshots (from mrview captures, same as v02)
Row 3: Per-mouse radar (4 mice × DV-dominant direction, first-thalamus slice)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from PIL import Image
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    get_directional_mean, slicewise_kmeans
)
from functions.config_loader import load_config

# ============================================================
# Config
# ============================================================
TMP = Path('tmp')
REPORTS = Path('reports')
REPORTS.mkdir(exist_ok=True)

cfg = load_config('config1.yaml')
# v08: use mrtrix_newmask pipeline instead of mrtrix_preproc_flip1
# Resolve data_base: /export/home → /home for non-lmc nodes
_data_base = Path(cfg['paths']['data_base'])
if not _data_base.exists() and str(_data_base).startswith('/export/home'):
    _data_base = Path(str(_data_base).replace('/export/home', '/home', 1))
DERIV = _data_base / 'data' / 'derivatives' / 'mrtrix_newmask'

# Representative subjects for FOD/Fixel screenshots (same as v02)
REPR_SUBJECTS = {
    'CTRL': {
        'id': 'sub-138-1010-L1F',
        'ses': 'ses-260115101441',
        'fod_slice': 14,
        'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz',
        'roi_name': 'isocortex',
    },
    'AP': {
        'id': 'sub-769-C-O5-R1M',
        'ses': 'ses-240911115733',
        'fod_slice': 15,
        'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif',
        'roi_name': 'graft',
    },
    'APTR': {
        'id': 'sub-AP72-862-901-NmF',
        'ses': 'ses-251210093342',
        'fod_slice': 14,
        'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif',
        'roi_name': 'graft',
    },
}

# All subjects per group with first-thalamus FOD slice
ALL_SUBJECTS = {
    'CTRL': [
        {'id': 'sub-129-921-R1F',   'ses': 'ses-260115121231', 'first_thal_fod': 10,
         'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz'},
        {'id': 'sub-129-921-R1L1F', 'ses': 'ses-260116085433', 'first_thal_fod': 10,
         'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz'},
        {'id': 'sub-129-1014-L1M',  'ses': 'ses-260115113814', 'first_thal_fod': 9,
         'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz'},
        {'id': 'sub-138-1010-L1F',  'ses': 'ses-260115101441', 'first_thal_fod': 9,
         'roi_rel': 'atlas/rois/fod_space/roi_isocortex_left.nii.gz'},
    ],
    'AP': [
        {'id': 'sub-596-C-OO-NmM', 'ses': 'ses-240911103349', 'first_thal_fod': 13,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-769-C-O5-R1M', 'ses': 'ses-240911115733', 'first_thal_fod': 12,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-769-C-OO-NmM', 'ses': 'ses-240911111853', 'first_thal_fod': 11,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-875-C-OO-NM',  'ses': 'ses-240911100253', 'first_thal_fod': 12,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
    ],
    'APTR': [
        {'id': 'sub-AP72-862-901-NmF', 'ses': 'ses-251210093342', 'first_thal_fod': 10,  # slice 9 is edge-only; DV cluster appears at 10
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-AP72-862-901-R1M', 'ses': 'ses-251210100043', 'first_thal_fod': 10,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-Apallial17NmF',    'ses': 'ses-240816091633', 'first_thal_fod': 10,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
        {'id': 'sub-Apallial20R1F',    'ses': 'ses-240909094232', 'first_thal_fod': 11,
         'roi_rel': 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'},
    ],
}

GROUP_COLORS = {'CTRL': '#333333', 'AP': '#800040', 'APTR': '#3377DD'}

# Screenshot file prefixes — v08 (re-captured from mrtrix_newmask)
SCREENSHOT_MAP = {
    ('CTRL', 'fod'):   TMP / 'ctrl_fod_v08_0000.png',
    ('AP', 'fod'):     TMP / 'ap_fod_v08_0000.png',
    ('APTR', 'fod'):   TMP / 'aptr_fod_v08_0000.png',
    ('CTRL', 'fixel'): TMP / 'ctrl_fixel_v08_0000.png',
    ('AP', 'fixel'):   TMP / 'ap_fixel_v08_0000.png',
    ('APTR', 'fixel'): TMP / 'aptr_fixel_v08_0000.png',
}

# Per-group mouse color shades (4 shades each)
GROUP_MOUSE_COLORS = {
    'CTRL': ['#1a1a1a', '#3d3d3d', '#5c5c5c', '#7a7a7a'],   # dark grey shades
    'AP':   ['#600030', '#800040', '#a0005c', '#c03080'],   # magenta shades
    'APTR': ['#1a3fa0', '#3366cc', '#5588dd', '#88aaee'],   # blue shades
}
MOUSE_MARKERS = ['o', 's', '^', 'D']


# ============================================================
# Helper: rotate screenshot 90° CW
# ============================================================
def load_rotated(path):
    img = Image.open(path).rotate(-90, expand=True)
    return np.array(img)


# ============================================================
# Helper: reduce radial grid — keep only ~1/3 of ticks
# ============================================================
def thin_radial_grid(ax, n_ticks=3):
    """Keep only n_ticks radial gridlines (removing ~2/3 of default)."""
    ylim = ax.get_ylim()
    ticks = np.linspace(0, ylim[1], n_ticks + 1)[1:]  # exclude 0
    ax.set_rticks(ticks)
    ax.set_yticklabels([f'{t:.2f}' for t in ticks], fontsize=7, color='gray')


# ============================================================
# Helper: draw radar for Row 3 — one line per mouse
# ============================================================
def draw_radar_mice(ax, mouse_data_list, title='', y_max=None, group='CTRL'):
    """Draw radar with two K-Means cluster traces per mouse.
    
    mouse_data_list: list of dicts with keys:
        'label', 'scaled_dir' (DV-dominant), 'scaled_dir_other', 'n_fixels'
    DV-dominant: solid line + marker.
    Other cluster: dashed line, same color, lower alpha, no marker.
    """
    colors = GROUP_MOUSE_COLORS.get(group, GROUP_MOUSE_COLORS['CTRL'])
    radar_labels = ['DV', 'ML', 'AP']
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML -> [DV, ML, AP]
    angles_closed = np.concatenate([angles, [angles[0]]])

    rmax = 0.0
    for k, md in enumerate(mouse_data_list):
        color = colors[k % len(colors)]
        # --- DV-dominant cluster (solid) ---
        sd = md['scaled_dir']
        abs_d = np.abs(sd)[reorder]
        rmax = max(rmax, np.max(abs_d))
        vals = np.concatenate([abs_d, [abs_d[0]]])
        ax.plot(angles_closed, vals, color=color, linewidth=1.3,
                linestyle='-', label=md['label'], alpha=0.75,
                marker=MOUSE_MARKERS[k % len(MOUSE_MARKERS)],
                markersize=5, zorder=len(mouse_data_list) - k + 2)
        ax.fill(angles_closed, vals, color=color, alpha=0.18)
        # --- Other cluster (dashed) ---
        sd2 = md.get('scaled_dir_other')
        if sd2 is not None:
            abs_d2 = np.abs(sd2)[reorder]
            rmax = max(rmax, np.max(abs_d2))
            vals2 = np.concatenate([abs_d2, [abs_d2[0]]])
            ax.plot(angles_closed, vals2, color=color, linewidth=1.0,
                    linestyle='--', alpha=0.45,
                    zorder=len(mouse_data_list) - k + 1)
            ax.fill(angles_closed, vals2, color=color, alpha=0.08)

    ax.set_theta_zero_location('S')
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(radar_labels, fontsize=10, fontweight='bold')

    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(0, max(rmax * 1.2, 0.01))

    thin_radial_grid(ax, n_ticks=3)

    ax.set_title(title, fontsize=10, pad=18)
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.10), fontsize=6.5,
              framealpha=0.8, handlelength=1.5)


# ============================================================
# Helper: draw radar for Row 4 — group avg (2 clusters, same as v02)
# ============================================================
def draw_radar_group(ax, directions_list, scaled_directions_list, mean_peaks,
                     cluster_sizes, angle_between=None, title='', y_max=None):
    radar_labels = ['DV', 'ML', 'AP']
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]
    angles_closed = np.concatenate([angles, [angles[0]]])

    colors = ['#E24A33', '#348ABD', '#988ED5']
    linestyles = ['-', '--', '-.']
    linewidths = [3.0, 1.5, 1.5]
    dir_labels = ['DV-dom', 'Other', 'Dir 3']

    rmax = 0.0
    for k, sd in enumerate(scaled_directions_list):
        abs_d = np.abs(sd)[reorder]
        rmax = max(rmax, np.max(abs_d))
        vals = np.concatenate([abs_d, [abs_d[0]]])
        ls = linestyles[k % len(linestyles)]
        lw = linewidths[k % len(linewidths)]
        n_k = cluster_sizes[k] if k < len(cluster_sizes) else 0
        ax.plot(angles_closed, vals, color=colors[k], linewidth=lw,
                linestyle=ls, label=f'{dir_labels[k]} (n={n_k})')
        ax.fill(angles_closed, vals, color=colors[k],
                alpha=0.20 if k == 0 else 0.08)

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

    thin_radial_grid(ax, n_ticks=3)

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
# Step 2: Generate per-mouse radar data (Row 3)
#   For each group, load all 4 mice.
#   For each mouse: extract fixels at first-thalamus slice,
#   K-means (2 clusters) → pick DV-dominant cluster's scaled_dir.
# ============================================================
print("\n--- Generating per-mouse radar data (Row 3) ---")
per_mouse_data = {}  # group → list of mouse dicts

for group in ['CTRL', 'AP', 'APTR']:
    mice_results = []
    print(f"\n  {group}:")
    for info in ALL_SUBJECTS[group]:
        sid = info['id']
        ses = info['ses']
        fod_slice = info['first_thal_fod']
        subj_dir = DERIV / sid / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        roi_path = subj_dir / info['roi_rel']

        if not (fixel_dir / 'index.mif').exists():
            print(f"    ⚠ Skipping {sid}: no fixels")
            continue
        if not roi_path.exists():
            print(f"    ⚠ Skipping {sid}: ROI not found")
            continue

        print(f"    {sid}: loading fixels, slice={fod_slice}...")
        fixel_data = load_fixel_dir(fixel_dir)
        roi_mask = load_roi_mask(roi_path)

        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)

        if fod_slice not in slice_vecs:
            # Try nearest available slice
            avail = sorted(slice_vecs.keys())
            if not avail:
                print(f"      ⚠ No slices with fixels in ROI!")
                continue
            nearest = min(avail, key=lambda s: abs(s - fod_slice))
            print(f"      ⚠ Slice {fod_slice} not available, using {nearest}")
            fod_slice = nearest

        vecs = slice_vecs[fod_slice]
        peaks = slice_peaks[fod_slice]
        n_vox = slice_voxels.get(fod_slice, len(vecs))

        if len(vecs) < 4:
            print(f"      ⚠ Only {len(vecs)} fixels, skipping")
            continue

        # K-Means (2 clusters) on all fixels in the graft/isocortex ROI at this slice.
        # Then pick only the DV-dominant cluster for plotting.
        km_results = slicewise_kmeans(
            {fod_slice: vecs}, n_clusters=2,
            slice_afds={fod_slice: peaks},
            slice_roi_voxels={fod_slice: n_vox},
        )
        km_slice = km_results[fod_slice]

        # Pick DV-dominant cluster: highest |scaled_dir[1]| (d[1]=DV)
        dv_vals = [np.abs(sd[1]) for sd in km_slice['scaled_directions']]
        best_k = int(np.argmax(dv_vals))
        other_k = 1 - best_k  # the other cluster
        scaled_dir = km_slice['scaled_directions'][best_k]
        scaled_dir_other = km_slice['scaled_directions'][other_k]
        n_cluster = km_slice['cluster_sizes'][best_k]
        n_other = km_slice['cluster_sizes'][other_k]

        # Short label: last meaningful part of subject ID
        short = sid.replace('sub-', '')
        mice_results.append({
            'label': short,
            'scaled_dir': scaled_dir,
            'scaled_dir_other': scaled_dir_other,
            'n_fixels': n_cluster,
        })
        print(f"      DV-dom scaled_dir={scaled_dir}, n={n_cluster}")
        print(f"      Other  scaled_dir={scaled_dir_other}, n={n_other}")

    per_mouse_data[group] = mice_results


# ============================================================
# Step 3: Compute consistent y_max across per-mouse radars
# ============================================================
all_maxes = []
for group, mice in per_mouse_data.items():
    for md in mice:
        all_maxes.append(np.max(np.abs(md['scaled_dir'])))
        if md.get('scaled_dir_other') is not None:
            all_maxes.append(np.max(np.abs(md['scaled_dir_other'])))
y_max_global = 1.0
print(f"\nGlobal y_max for radar plots: {y_max_global:.4f}")


# ============================================================
# Step 5: Compose 4×3 figure
# ============================================================
print("\n--- Composing 3×3 figure ---")

fig = plt.figure(figsize=(18, 20))
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.20,
                       height_ratios=[1, 1, 1.0],
                       top=0.93, bottom=0.04, left=0.08, right=0.95)

groups_order = ['CTRL', 'AP', 'APTR']
group_titles = ['Control (isocortex)', 'Apallial (graft)', '$t^{AP}$-hCO (graft + Tx)']

# Row 1: FOD screenshots
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[0, col])
    ax.imshow(screenshots.get((group, 'fod'), np.zeros((10, 10, 3), dtype=np.uint8)))
    ax.set_title(f'{group_titles[col]}\n{REPR_SUBJECTS[group]["id"]}', fontsize=11, fontweight='bold', pad=10)
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

# Row 3: Per-mouse radar (4 mice, DV-dominant direction)
GROUP_DISPLAY = {'CTRL': 'Control', 'AP': 'Apallial', 'APTR': '$t^{AP}$-hCO'}
for col, group in enumerate(groups_order):
    ax = fig.add_subplot(gs[2, col], projection='polar')
    mice = per_mouse_data.get(group, [])
    title = GROUP_DISPLAY[group]
    if mice:
        draw_radar_mice(ax, mice, title=title, y_max=y_max_global, group=group)
    else:
        ax.set_title(f'{group}: no data', fontsize=10)

# Row labels on left side
fig.text(0.02, 0.78, 'Row 1: FOD Glyphs', fontsize=12, fontweight='bold',
         rotation=90, va='center')
fig.text(0.02, 0.50, 'Row 2: Fixel Sticks', fontsize=12, fontweight='bold',
         rotation=90, va='center')
fig.text(0.02, 0.20, 'Row 3: Individual Mice\n(first-thalamus slice)', fontsize=11,
         fontweight='bold', rotation=90, va='center')

fig.suptitle('Fixel Direction Analysis — Control vs Apallial vs $t^{AP}$-hCO (left hemisphere) v08 [newmask]',
             fontsize=16, fontweight='bold', y=0.97)

save_path = REPORTS / 'fixel_4x3_comparison_v08.png'
fig.savefig(save_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"\n✅ Saved: {save_path}")

save_path_svg = REPORTS / 'fixel_4x3_comparison_v08.svg'
try:
    import cairo  # pycairo — enables matplotlib cairo backend
    from matplotlib.backends.backend_cairo import FigureCanvasCairo
    canvas_orig = fig.canvas
    FigureCanvasCairo(fig)        # attach cairo canvas
    fig.savefig(save_path_svg, format='svg', bbox_inches='tight', facecolor='white')
    fig.set_canvas(canvas_orig)   # restore original canvas
    print(f"✅ Saved (cairo SVG): {save_path_svg}")
except Exception as e:
    print(f"⚠️  Cairo SVG failed ({e}), using default SVG with fonttype=none")
    # Use fonttype='none' so text stays as <text> elements (not paths),
    # and skip bbox_inches='tight' which mis-calculates text bboxes in SVG.
    with matplotlib.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(save_path_svg, format='svg', facecolor='white')
    print(f"✅ Saved (SVG fonttype=none): {save_path_svg}")
plt.close(fig)
