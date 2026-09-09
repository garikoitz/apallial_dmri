#!/usr/bin/env python3
"""Generate fixel_3group v10 figures: bilateral ROIs, newmask pipeline.

Changes from v09:
  - Radar plots match the style from fixel_4x3_comparison_v09:
      * Per-group color shades (grey=CTRL, magenta=AP, blue=APTR)
      * Markers per mouse (o, s, ^, D)
      * White background for ALL radar plots
      * thin_radial_grid with 3 ticks
      * Larger axis labels (fontsize=10)
      * Adjusted line widths and alphas
  - Slice 0 removed: data starts at slice 1 (first slice dropped).
  - Outputs both PNG and SVG for each figure.

Usage:
    cd /home/glerma/glerma/soft/mbrain
    python3 tmp/make_3group_v10.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans
)
from functions.config_loader import load_config

# ============================================================
# Config
# ============================================================
REPORTS = Path('reports')
REPORTS.mkdir(exist_ok=True)

cfg = load_config('config1.yaml')
_data_base = Path(cfg['paths']['data_base'])
if not _data_base.exists() and str(_data_base).startswith('/export/home'):
    _data_base = Path(str(_data_base).replace('/export/home', '/home', 1))
DERIV = _data_base / 'data' / 'derivatives' / 'mrtrix_newmask'

# Subject definitions (same as config1.yaml groups, excluding test)
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

# Per-group mouse color shades (matching 4x3 v09 style)
GROUP_MOUSE_COLORS = {
    'CTRL': ['#1a1a1a', '#3d3d3d', '#5c5c5c', '#7a7a7a'],   # dark grey shades
    'AP':   ['#600030', '#800040', '#a0005c', '#c03080'],     # magenta shades
    'APTR': ['#1a3fa0', '#3366cc', '#5588dd', '#88aaee'],     # blue shades
}
MOUSE_MARKERS = ['o', 's', '^', 'D']

# ROI paths (bilateral): for each hemisphere, we load _left + _right
def get_roi_paths(subj_dir, group):
    """Return (left_path, right_path) for ROI files."""
    if group == 'CTRL':
        base = subj_dir / 'atlas' / 'rois' / 'fod_space'
        return (base / 'roi_isocortex_left.nii.gz',
                base / 'roi_isocortex_right.nii.gz')
    else:
        base = subj_dir / 'manual_roi' / 'fod_space'
        return (base / 'roi_Graft_left_flip_fod.mif',
                base / 'roi_Graft_right_flip_fod.mif')


# ============================================================
# Helper: thin_radial_grid (from 4x3 v09)
# ============================================================
def thin_radial_grid(ax, n_ticks=3):
    """Keep only n_ticks radial gridlines (removing ~2/3 of default)."""
    ylim = ax.get_ylim()
    ticks = np.linspace(0, ylim[1], n_ticks + 1)[1:]  # exclude 0
    ax.set_rticks(ticks)
    ax.set_yticklabels([f'{t:.2f}' for t in ticks], fontsize=7, color='gray')


# ============================================================
# Step 1: Load all fixel data with bilateral ROIs
# ============================================================
def load_all_group_data(verbose=True):
    """Load fixels for all 12 subjects using bilateral (L|R) ROI masks.

    Returns
    -------
    all_group_data : dict
        {group: {subj_id: {slice_idx: kmeans_result, ...}, ...}, ...}
    target_n_slices : int
        Minimum number of slices across all AP/APTR graft subjects,
        minus 1 (because we drop the first slice).
    """
    all_group_data = {}
    graft_slice_counts = []

    # Process AP and APTR first to determine target_n_slices
    for group in ['AP', 'APTR']:
        if verbose:
            print(f'\n--- {group} (bilateral, newmask) ---')
        all_group_data[group] = {}

        for (subj, ses) in SUBJECT_MAP[group]:
            subj_dir = DERIV / subj / ses
            fixel_dir = subj_dir / 'dwi' / 'fixels_group'

            if not (fixel_dir / 'index.mif').exists():
                if verbose:
                    print(f'  ⚠ Skipping {subj}: no fixels')
                continue

            roi_left, roi_right = get_roi_paths(subj_dir, group)
            if not roi_left.exists() or not roi_right.exists():
                if verbose:
                    print(f'  ⚠ Skipping {subj}: ROI L/R not found')
                continue

            if verbose:
                print(f'  Loading {subj}...')
            fixel_data = load_fixel_dir(fixel_dir)

            # Bilateral mask: L | R
            mask_l = load_roi_mask(roi_left)
            mask_r = load_roi_mask(roi_right)
            roi_mask = mask_l | mask_r
            n_l, n_r = int(mask_l.sum()), int(mask_r.sum())
            if verbose:
                print(f'    ROI voxels: L={n_l}, R={n_r}, L|R={int(roi_mask.sum())}')

            slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
                fixel_data, roi_mask, slice_axis=0)
            slice_km = slicewise_kmeans(slice_vecs, n_clusters=2,
                                        slice_afds=slice_peaks,
                                        slice_roi_voxels=slice_voxels)

            if slice_km:
                all_group_data[group][subj] = slice_km
                graft_slice_counts.append(len(slice_km))
                if verbose:
                    print(f'    {len(slice_km)} slices with fixel data')

    if not graft_slice_counts:
        print('  ⚠ No graft fixel data found!')
        return None, 0

    # target_n_slices: smallest graft minus 1 (we drop the first slice)
    raw_target = min(graft_slice_counts)
    target_n_slices = raw_target - 1
    if verbose:
        print(f'\n  Raw smallest graft slice count: {raw_target}')
        print(f'  After dropping slice 0: target_n_slices = {target_n_slices}')

    # Re-index AP/APTR: skip first slice (index 0), keep next target_n_slices
    for group in ['AP', 'APTR']:
        for subj, subj_data in all_group_data[group].items():
            sorted_slices = sorted(subj_data.keys())
            # Drop the first slice, then take target_n_slices
            kept_slices = sorted_slices[1:1 + target_n_slices]
            all_group_data[group][subj] = {i: subj_data[s]
                                            for i, s in enumerate(kept_slices)}

    # Process CTRL
    if verbose:
        print(f'\n--- CTRL (bilateral, newmask) ---')
    all_group_data['CTRL'] = {}

    for (subj, ses) in SUBJECT_MAP['CTRL']:
        subj_dir = DERIV / subj / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'

        if not (fixel_dir / 'index.mif').exists():
            if verbose:
                print(f'  ⚠ Skipping {subj}: no fixels')
            continue

        roi_left, roi_right = get_roi_paths(subj_dir, 'CTRL')
        if not roi_left.exists() or not roi_right.exists():
            if verbose:
                print(f'  ⚠ Skipping {subj}: ROI L/R not found')
            continue

        if verbose:
            print(f'  Loading {subj}...')
        fixel_data = load_fixel_dir(fixel_dir)

        mask_l = load_roi_mask(roi_left)
        mask_r = load_roi_mask(roi_right)
        roi_mask = mask_l | mask_r
        n_l, n_r = int(mask_l.sum()), int(mask_r.sum())
        if verbose:
            print(f'    ROI voxels: L={n_l}, R={n_r}, L|R={int(roi_mask.sum())}')

        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(
            fixel_data, roi_mask, slice_axis=0)
        slice_km = slicewise_kmeans(slice_vecs, n_clusters=2,
                                    slice_afds=slice_peaks,
                                    slice_roi_voxels=slice_voxels)

        if slice_km:
            sorted_slices = sorted(slice_km.keys())
            # Drop the first slice, then take target_n_slices
            kept_slices = sorted_slices[1:1 + target_n_slices]
            all_group_data['CTRL'][subj] = {i: slice_km[s]
                                             for i, s in enumerate(kept_slices)}
            if verbose:
                print(f'    {len(kept_slices)} slices with fixel data (after drop)')

    return all_group_data, target_n_slices


# ============================================================
# Compute global max from all_group_data
# ============================================================
def compute_global_max(all_group_data):
    """Find the maximum |scaled_dir| component across all groups/subjects/slices."""
    global_max = 0.0
    for group in ['CTRL', 'AP', 'APTR']:
        for subj, data in all_group_data.get(group, {}).items():
            for s_idx, res in data.items():
                for sd in res.get('scaled_directions', []):
                    val = np.max(np.abs(sd))
                    if val > global_max:
                        global_max = val
    return global_max


# ============================================================
# Plot: 3-group radar comparison v10 (4x3-style radars)
# ============================================================
def plot_3group_radar_v10(all_group_data, target_n_slices, y_max,
                          save_png=None, save_svg=None):
    """Generate the 3-group K-means radar comparison figure.

    v10 changes vs v09:
    - Per-group color shades (matching 4x3 style)
    - Markers per mouse
    - White background for all radars
    - thin_radial_grid (3 ticks)
    - Larger axis labels
    - Slice 0 already dropped in data loading
    """
    row_groups = ['CTRL', 'AP', 'APTR']
    group_labels = ['CTRL\n(isocortex)', 'AP\n(graft, no Tx)', 'XCX\n(graft + Tx)']

    radar_labels_list = ['DV', 'ML', 'AP']
    angles_radar = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML -> [DV, ML, AP]
    angles_closed = np.concatenate([angles_radar, [angles_radar[0]]])

    n_cols = target_n_slices + 3  # slices + AVG/subj + GROUP AVG + INFO
    n_rows = len(row_groups)

    fig = plt.figure(figsize=(2.4 * n_cols, 3.5 * n_rows))

    def _draw_radar(ax, directions_list, color, marker='o',
                    linewidth_solid=1.3,
                    alpha_solid=0.75,
                    fill_solid=0.18,
                    show_marker=True, zorder_base=2):
        # Only draw Dir 1 (solid/filled), skip Dir 2 (dashed)
        d = directions_list[0]
        abs_d = np.abs(d)[reorder]
        vals = np.concatenate([abs_d, [abs_d[0]]])
        kw = dict(color=color, linewidth=linewidth_solid,
                  linestyle='-', alpha=alpha_solid,
                  zorder=zorder_base + 1)
        if show_marker:
            kw['marker'] = marker
            kw['markersize'] = 5
        ax.plot(angles_closed, vals, **kw)
        ax.fill(angles_closed, vals, color=color, alpha=fill_solid)

    def _setup_radar_axes(ax, show_labels=False, show_rtick_numbers=False):
        ax.set_facecolor('white')
        ax.set_theta_zero_location('S')
        ax.set_theta_direction(-1)
        ax.set_ylim(0, y_max)
        if show_labels:
            ax.set_xticks(angles_radar)
            ax.set_xticklabels(radar_labels_list, fontsize=10, fontweight='bold')
        else:
            ax.set_xticks(angles_radar)
            ax.set_xticklabels(['', '', ''])
        ax.grid(True, alpha=0.25)
        thin_radial_grid(ax, n_ticks=3)
        if not show_rtick_numbers:
            ax.set_yticklabels([])

    for row_idx, (group, label) in enumerate(zip(row_groups, group_labels)):
        subjects_data = all_group_data.get(group, {})
        subject_ids = list(subjects_data.keys())
        colors = GROUP_MOUSE_COLORS.get(group, GROUP_MOUSE_COLORS['CTRL'])

        # Per-slice columns
        for col_idx in range(target_n_slices):
            ax_idx = row_idx * n_cols + col_idx + 1
            ax = fig.add_subplot(n_rows, n_cols, ax_idx, projection='polar')
            # Show radial tick numbers only on the first slice of the AP row
            show_rticks = (row_idx == 1 and col_idx == 0)
            _setup_radar_axes(ax, show_labels=(col_idx == 0),
                              show_rtick_numbers=show_rticks)

            for si, subj in enumerate(subject_ids):
                color = colors[si % len(colors)]
                marker = MOUSE_MARKERS[si % len(MOUSE_MARKERS)]
                data = subjects_data[subj]
                if col_idx in data:
                    dirs_key = 'scaled_directions' if 'scaled_directions' in data[col_idx] else 'directions'
                    _draw_radar(ax, data[col_idx][dirs_key], color, marker=marker,
                                zorder_base=len(subject_ids) - si + 2)

            if row_idx == 0:
                # Label slices starting from 1 (since we dropped original slice 0)
                ax.set_title(f'Slice {col_idx + 1}', fontsize=9, fontweight='bold', pad=8)

        # AVG/subj column
        ax_avg = fig.add_subplot(n_rows, n_cols,
                                 row_idx * n_cols + target_n_slices + 1, projection='polar')
        _setup_radar_axes(ax_avg)

        for si, subj in enumerate(subject_ids):
            color = colors[si % len(colors)]
            marker = MOUSE_MARKERS[si % len(MOUSE_MARKERS)]
            data = subjects_data[subj]
            sd1_list = [data[s]['scaled_directions'][0] for s in range(target_n_slices)
                        if s in data and 'scaled_directions' in data[s]]
            sd2_list = [data[s]['scaled_directions'][1] for s in range(target_n_slices)
                        if s in data and 'scaled_directions' in data[s]
                        and len(data[s]['scaled_directions']) > 1]
            avg_dirs = []
            if sd1_list:
                avg_dirs.append(np.mean(sd1_list, axis=0))
            if sd2_list:
                avg_dirs.append(np.mean(sd2_list, axis=0))
            if avg_dirs:
                _draw_radar(ax_avg, [avg_dirs[0]], color, marker=marker,
                            zorder_base=len(subject_ids) - si + 2)

        if row_idx == 0:
            ax_avg.set_title('AVG/subj', fontsize=8, fontweight='bold', pad=8, color='darkred')

        # GROUP AVG column
        ax_grp = fig.add_subplot(n_rows, n_cols,
                                 row_idx * n_cols + target_n_slices + 2, projection='polar')
        _setup_radar_axes(ax_grp)

        all_sd1, all_sd2 = [], []
        for subj in subject_ids:
            data = subjects_data[subj]
            for s in range(target_n_slices):
                if s not in data:
                    continue
                if 'scaled_directions' in data[s]:
                    all_sd1.append(data[s]['scaled_directions'][0])
                    if len(data[s]['scaled_directions']) > 1:
                        all_sd2.append(data[s]['scaled_directions'][1])

        grp_dirs = []
        grp_dirs_unit = []
        if all_sd1:
            grp_sd1 = np.mean(all_sd1, axis=0)
            grp_dirs.append(grp_sd1)
            norm1 = np.linalg.norm(grp_sd1)
            grp_dirs_unit.append(grp_sd1 / norm1 if norm1 > 1e-8 else grp_sd1)
        if all_sd2:
            grp_sd2 = np.mean(all_sd2, axis=0)
            grp_dirs.append(grp_sd2)
            norm2 = np.linalg.norm(grp_sd2)
            grp_dirs_unit.append(grp_sd2 / norm2 if norm2 > 1e-8 else grp_sd2)
        if grp_dirs:
            _draw_radar(ax_grp, [grp_dirs[0]], 'black',
                        linewidth_solid=3,
                        alpha_solid=0.85,
                        fill_solid=0.15,
                        show_marker=False)

        if row_idx == 0:
            ax_grp.set_title('GROUP\nAVG', fontsize=8, fontweight='bold', pad=8, color='darkred')

        # INFO panel
        ax_info = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + n_cols)
        ax_info.set_facecolor('white')
        ax_info.axis('off')

        n_subj = len(subject_ids)
        all_angles = []
        for subj in subject_ids:
            for s in range(target_n_slices):
                if s in subjects_data[subj] and 'angle_between' in subjects_data[subj][s]:
                    all_angles.append(subjects_data[subj][s]['angle_between'])
        mean_angle = np.mean(all_angles) if all_angles else 0

        all_norms = [np.linalg.norm(sd) for sd in all_sd1 + all_sd2]
        grp_mean_peak = float(np.mean(all_norms)) if all_norms else 0.0

        info_text = f'{label}\n\nn={n_subj} subjects\nmean Peak={grp_mean_peak:.3f}\n'
        if grp_dirs_unit:
            g1 = grp_dirs_unit[0]
            info_text += f'\nDir1: [{g1[0]:+.2f},{g1[1]:+.2f},{g1[2]:+.2f}]'
            dominant = ['AP', 'DV', 'ML'][int(np.argmax(np.abs(g1)))]
            info_text += f'\nDominant: {dominant}'
        if len(grp_dirs_unit) > 1:
            g2 = grp_dirs_unit[1]
            info_text += f'\n\nDir2: [{g2[0]:+.2f},{g2[1]:+.2f},{g2[2]:+.2f}]'
            dominant2 = ['AP', 'DV', 'ML'][int(np.argmax(np.abs(g2)))]
            info_text += f'\nDominant: {dominant2}'
        if all_angles:
            info_text += f'\n\n∠ Dir1-Dir2:\n{mean_angle:.1f}° (mean)'

        ax_info.text(0.5, 0.5, info_text, transform=ax_info.transAxes,
                     fontsize=8, ha='center', va='center',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Legend — per-group style
    legend_elements = []
    for group in row_groups:
        colors = GROUP_MOUSE_COLORS[group]
        subjects = list(all_group_data.get(group, {}).keys())
        for i, subj in enumerate(subjects):
            short = subj.replace('sub-', '')[:15]
            legend_elements.append(plt.Line2D(
                [0], [0], color=colors[i % len(colors)], linewidth=2,
                marker=MOUSE_MARKERS[i % len(MOUSE_MARKERS)],
                label=f'{group} {short}'))
    legend_elements.append(plt.Line2D([0], [0], color='black', linewidth=1.3,
                                       linestyle='-', label='Dir 1 (solid)'))
    fig.legend(handles=legend_elements, loc='upper right', fontsize=7,
               bbox_to_anchor=(0.99, 0.98), title='Subjects', ncol=1)

    fig.suptitle(
        f'Fixel K-means Direction × Peak Amplitude — Bilateral (L+R) v10 [newmask]\n'
        f'Radar scaled by FOD peak amplitude | Solid = Dir 1 only | '
        f'Outer ring = {y_max:.2f} | Slices start at 1 (slice 0 removed)',
        fontsize=13, fontweight='bold', y=0.99)

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])

    if save_png:
        fig.savefig(str(save_png), dpi=150, bbox_inches='tight', facecolor='white')
        print(f'  ✅ Saved PNG: {save_png}')
    if save_svg:
        with matplotlib.rc_context({'svg.fonttype': 'none'}):
            fig.savefig(str(save_svg), format='svg', facecolor='white')
        print(f'  ✅ Saved SVG: {save_svg}')

    plt.close(fig)
    return fig


# ============================================================
# Plot: 3-group component lines (same as v09 but slices start at 1)
# ============================================================
def plot_3group_components_v10(all_group_data, target_n_slices, y_max,
                                save_png=None, save_svg=None):
    """Plot per-slice directional component magnitudes for all 3 groups.

    Same as v09 but:
    - Slice labels start at 1 (slice 0 was dropped in data loading)
    """
    row_groups = ['CTRL', 'AP', 'APTR']
    group_titles = ['CTRL (isocortex)', 'AP (graft, no Tx)', 'XCX (graft + Tx)']
    group_bg = ['#e8f5e9', '#ffebee', '#e3f2fd']
    subject_markers = ['o', 's', '^', 'D']

    comp_colors = {'AP': '#e74c3c', 'DV': '#3498db', 'ML': '#27ae60'}
    comp_idx = {'AP': 0, 'DV': 1, 'ML': 2}

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False, sharey=True)

    for row_idx, (group, title, bg) in enumerate(zip(row_groups, group_titles, group_bg)):
        ax = axes[row_idx]
        ax.set_facecolor(bg)
        subjects_data = all_group_data.get(group, {})
        subject_ids = list(subjects_data.keys())

        max_slices = 0
        for si, subj in enumerate(subject_ids):
            data = subjects_data[subj]
            slices_sorted = sorted(data.keys())
            n_sl = len(slices_sorted)
            if n_sl > max_slices:
                max_slices = n_sl

            # x-axis labels start at 1
            x = np.arange(n_sl) + 1
            marker = subject_markers[si % len(subject_markers)]

            for comp_name, c_idx in comp_idx.items():
                color = comp_colors[comp_name]
                dirs_key = 'scaled_directions' if 'scaled_directions' in data.get(slices_sorted[0], {}) else 'directions'
                y_vals = []
                for s in slices_sorted:
                    sd = data[s].get(dirs_key, data[s].get('directions', [None]))[0]
                    if sd is not None:
                        y_vals.append(abs(sd[c_idx]))
                    else:
                        y_vals.append(0.0)

                y_vals = np.array(y_vals)
                label = f'{comp_name}' if si == 0 else None
                alpha = max(0.3, 1.0 - si * 0.15)
                ax.plot(x, y_vals, color=color, marker=marker, markersize=4,
                        linewidth=1.5, alpha=alpha, label=label)

        ax.set_ylabel('Peak × |component|', fontsize=10)
        ax.set_ylim(0, y_max)
        ax.set_title(f'{title}  (n={len(subject_ids)})', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        if max_slices > 0:
            ax.set_xlim(0.5, max(max_slices, target_n_slices) + 0.5)
            ax.set_xticks(range(1, max(max_slices, target_n_slices) + 1))

        if row_idx == 0:
            ax.legend(loc='upper right', fontsize=9, ncol=3,
                      title='Component (Dir1)', title_fontsize=8)

        legend_str = '  '.join([f'{subject_markers[i % len(subject_markers)]} {subj.replace("sub-", "")[:15]}'
                                 for i, subj in enumerate(subject_ids)])
        ax.text(0.01, 0.97, legend_str, transform=ax.transAxes,
                fontsize=6, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    axes[-1].set_xlabel('Coronal slice index (anterior → posterior, starting at 1)', fontsize=11)

    fig.suptitle(
        f'Directional Components per Slice — Bilateral (L+R) v10 [newmask]\n'
        f'Dir1 scaled by FOD peak amplitude | AP (red), DV (blue), ML (green) | Slice 0 removed',
        fontsize=13, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_png:
        fig.savefig(str(save_png), dpi=150, bbox_inches='tight', facecolor='white')
        print(f'  ✅ Saved PNG: {save_png}')
    if save_svg:
        with matplotlib.rc_context({'svg.fonttype': 'none'}):
            fig.savefig(str(save_svg), format='svg', facecolor='white')
        print(f'  ✅ Saved SVG: {save_svg}')

    plt.close(fig)
    return fig


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print('=' * 70)
    print('Fixel 3-group comparison v10 — bilateral ROIs, newmask pipeline')
    print('  Changes: 4x3-style radars, slice 0 removed')
    print('=' * 70)

    # Load all data (bilateral ROIs, slice 0 dropped)
    all_group_data, target_n_slices = load_all_group_data(verbose=True)

    if all_group_data is None:
        print('ERROR: No data loaded!')
        sys.exit(1)

    # Compute global max for y_max (outer ring = max value, no padding)
    y_max = compute_global_max(all_group_data)
    print(f'\n  Global y_max (outer ring): {y_max:.4f}')

    # Generate radar comparison figure
    print('\n--- Generating radar comparison ---')
    plot_3group_radar_v10(
        all_group_data, target_n_slices, y_max,
        save_png=REPORTS / 'fixel_3group_comparison_bilateral_v10.png',
        save_svg=REPORTS / 'fixel_3group_comparison_bilateral_v10.svg',
    )

    print('\n✅ All done!')
