"""Helper functions for figure generation in generate_figures_only.py.

This module collects all the plotting / data-loading helpers that were
previously inlined in the notebook cells, so the notebook stays focused on
the high-level figure composition while the reusable logic lives here.

Two figure-specific helper groups are provided:

* **Figure 2h** — screenshot rotation, radar drawing, scale bar, dashed margin.
* **Figure 3k** — ROI path resolution, fixel input validation, bilateral
  data loading, global-max computation, and the 3-group radar plot.

The figure-2h helpers accept the per-figure constants (color maps, markers,
group order) as arguments so they remain self-contained and testable.
"""

from __future__ import annotations

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
from pathlib import Path

from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans,
)


# ============================================================
# Figure 2h helpers
# ============================================================
def load_rotated(path):
    """Open PNG and rotate 90° CW (mrview captures are sideways)."""
    img = Image.open(path).rotate(-90, expand=True)
    return np.array(img)


def thin_radial_grid(ax, n_ticks=3):
    """Keep only n_ticks radial gridlines."""
    ylim = ax.get_ylim()
    ticks = np.linspace(0, ylim[1], n_ticks + 1)[1:]
    ax.set_rticks(ticks)
    ax.set_yticklabels([f"{t:.2f}" for t in ticks], fontsize=7, color="gray")


def draw_radar_mice(ax, mouse_data_list, title="", y_max=None, group="CTRL",
                    group_mouse_colors=None, mouse_markers=None):
    """Draw per-mouse radar: DV-dominant cluster (solid+marker) + other (dashed).

    Parameters
    ----------
    group_mouse_colors : dict[str, list[str]] | None
        Mapping ``{group: [color_0, color_1, ...]}``.  Falls back to a
        sensible default when *None*.
    mouse_markers : list[str] | None
        Marker styles cycled per mouse.  Falls back to
        ``["o", "s", "^", "D"]`` when *None*.
    """
    if group_mouse_colors is None:
        group_mouse_colors = {
            "CTRL": ["#1a1a1a", "#3d3d3d", "#5c5c5c", "#7a7a7a"],
            "AP":   ["#600030", "#800040", "#a0005c", "#c03080"],
            "APTR": ["#1a3fa0", "#3366cc", "#5588dd", "#88aaee"],
        }
    if mouse_markers is None:
        mouse_markers = ["o", "s", "^", "D"]

    colors = group_mouse_colors.get(group, group_mouse_colors["CTRL"])
    radar_labels = ["DV", "ML", "AP"]
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML → [DV, ML, AP]
    angles_closed = np.concatenate([angles, [angles[0]]])

    rmax = 0.0
    for k, md in enumerate(mouse_data_list):
        color = colors[k % len(colors)]
        # DV-dominant cluster (solid)
        sd = md["scaled_dir"]
        abs_d = np.abs(sd)[reorder]
        rmax = max(rmax, np.max(abs_d))
        vals = np.concatenate([abs_d, [abs_d[0]]])
        ax.plot(angles_closed, vals, color=color, linewidth=1.3,
                linestyle="-", label=md["label"], alpha=0.75,
                marker=mouse_markers[k % len(mouse_markers)],
                markersize=5, zorder=len(mouse_data_list) - k + 2)
        ax.fill(angles_closed, vals, color=color, alpha=0.18)
        # # Other cluster (dashed)
        # sd2 = md.get("scaled_dir_other")
        # if sd2 is not None:
        #     abs_d2 = np.abs(sd2)[reorder]
        #     rmax = max(rmax, np.max(abs_d2))
        #     vals2 = np.concatenate([abs_d2, [abs_d2[0]]])
        #     ax.plot(angles_closed, vals2, color=color, linewidth=1.0,
        #             linestyle="--", alpha=0.45,
        #             zorder=len(mouse_data_list) - k + 1)
        #     ax.fill(angles_closed, vals2, color=color, alpha=0.08)

    ax.set_theta_zero_location("S")
    ax.set_theta_direction(-1)
    ax.set_xticks(angles)
    ax.set_xticklabels(radar_labels, fontsize=10, fontweight="bold")
    if y_max is not None:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(0, max(rmax * 1.2, 0.01))
    thin_radial_grid(ax, n_ticks=3)
    ax.set_title(title, fontsize=10, pad=18)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.10), fontsize=6.5,
              framealpha=0.8, handlelength=1.5)


def add_scale_bar(ax, bar_mm=1.0, px_per_mm=None, position="br", color="white",
                  linewidth=2.5, pad_frac=0.05):
    """Add a scale bar to an image axis.

    px_per_mm: how many image pixels correspond to 1 mm. If None, auto-estimate
    from the image height assuming ~755 px ≈ 160 voxels × 0.1 mm/voxel = 16 mm
    field-of-view (i.e. px_per_mm ≈ 47). Adjust per actual capture.
    """
    if px_per_mm is None:
        # Heuristic: screenshots are ~755 px tall after rotation, covering
        # roughly 16 mm of tissue (160 slices × 0.1 mm). So ~47 px/mm.
        px_per_mm = 47.0
    bar_px = bar_mm * px_per_mm
    img_w = ax.get_xlim()[1]
    img_h = ax.get_ylim()[0]
    pad_x = img_w * pad_frac
    pad_y = img_h * pad_frac
    if position == "br":
        x0 = img_w - pad_x - bar_px
        y0 = img_h - pad_y
    elif position == "bl":
        x0 = pad_x
        y0 = img_h - pad_y
    elif position == "tr":
        x0 = img_w - pad_x - bar_px
        y0 = pad_y
    else:  # tl
        x0 = pad_x
        y0 = pad_y
    ax.plot([x0, x0 + bar_px], [y0, y0], color=color, linewidth=linewidth,
            solid_capstyle="butt", zorder=10)
    ax.text(x0 + bar_px / 2, y0 - img_h * 0.03, f"{bar_mm:.0f} mm",
            color=color, fontsize=8, ha="center", va="bottom",
            fontweight="bold", zorder=10)


def add_dashed_margin(ax, color="white", linewidth=1.5, linestyle="--", inset_frac=0.06):
    """Draw a dashed rectangle inset from the image border to indicate ROI margin."""
    img_w = ax.get_xlim()[1]
    img_h = ax.get_ylim()[0]
    dx = img_w * inset_frac
    dy = img_h * inset_frac
    rect = Rectangle((dx, dy), img_w - 2 * dx, img_h - 2 * dy,
                     linewidth=linewidth, edgecolor=color, facecolor="none",
                     linestyle=linestyle, zorder=9)
    ax.add_patch(rect)


# ============================================================
# Figure 3k helpers
# ============================================================
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


def missing_fixel_inputs(fixel_dir):
    """Return missing required files for an existing fod2fixel directory."""
    missing = []
    for name in ('index.mif', 'directions.mif'):
        if not (fixel_dir / name).exists():
            missing.append(name)
    if not (fixel_dir / 'afd.mif').exists() and not (fixel_dir / 'fd.mif').exists():
        missing.append('afd.mif or fd.mif')
    return missing


def validate_figure_3k_inputs(deriv, subject_map, max_examples=18):
    """Fail early when the installed data cannot support Figure 3k.

    Parameters
    ----------
    deriv : pathlib.Path
        Root derivatives directory (``mrtrix_newmask``).
    subject_map : dict[str, list[tuple[str, str]]]
        ``{group: [(subject_id, session_id), ...]}``.
    """
    missing = []

    for group, subjects in subject_map.items():
        for subj, ses in subjects:
            subj_dir = deriv / subj / ses
            if not subj_dir.exists():
                missing.append(f'{group} {subj}: missing subject/session directory')
                continue

            fixel_dir = subj_dir / 'dwi' / 'fixels_group'
            fixel_missing = missing_fixel_inputs(fixel_dir)
            if fixel_missing:
                missing.append(
                    f'{group} {subj}: missing dwi/fixels_group/'
                    f'{{{", ".join(fixel_missing)}}}'
                )

            roi_left, roi_right = get_roi_paths(subj_dir, group)
            for roi_path in (roi_left, roi_right):
                if not roi_path.exists():
                    missing.append(
                        f'{group} {subj}: missing '
                        f'{roi_path.relative_to(subj_dir)}'
                    )

    if not missing:
        return

    shown = '\n'.join(f'  - {item}' for item in missing[:max_examples])
    extra = len(missing) - max_examples
    if extra > 0:
        shown += f'\n  - ... {extra} more missing input(s)'

    raise RuntimeError(
        'Figure 3k inputs are incomplete under:\n'
        f'  {deriv}\n\n'
        'The 3k cell expects the post mrtrix_newmask outputs: per-subject '
        'dwi/fixels_group files plus FOD-space bilateral ROI masks.\n\n'
        f'Missing examples:\n{shown}\n\n'
        'This usually means the pre_preproc component is installed but the post '
        'component is missing or stale. Reinstall the post component with:\n'
        '  python3 utils/setup_data.py --components post --force\n\n'
        'If the post data lives outside this repo, set MBRAIN_DATA_BASE or '
        'MBRAIN_DERIVATIVES_DIR before running the Setup cell.'
    )


def load_all_group_data(deriv, subject_map, verbose=True):
    """Load fixels for all 12 subjects using bilateral (L|R) ROI masks.

    Parameters
    ----------
    deriv : pathlib.Path
        Root derivatives directory.
    subject_map : dict[str, list[tuple[str, str]]]
        ``{group: [(subject_id, session_id), ...]}``.

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

        for (subj, ses) in subject_map[group]:
            subj_dir = deriv / subj / ses
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

    for (subj, ses) in subject_map['CTRL']:
        subj_dir = deriv / subj / ses
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


def plot_3group_radar_v11(all_group_data, target_n_slices, y_max,
                          group_mouse_colors=None, mouse_markers=None,
                          save_png=None, save_svg=None):
    """Generate the 3-group K-means radar comparison figure (v11 style).

    Parameters
    ----------
    group_mouse_colors : dict[str, list[str]] | None
        Per-group color shades.  Falls back to a sensible default.
    mouse_markers : list[str] | None
        Marker styles cycled per mouse.
    """
    if group_mouse_colors is None:
        group_mouse_colors = {
            'CTRL': ['#1a1a1a', '#3d3d3d', '#5c5c5c', '#7a7a7a'],
            'AP':   ['#600030', '#800040', '#a0005c', '#c03080'],
            'APTR': ['#1a3fa0', '#3366cc', '#5588dd', '#88aaee'],
        }
    if mouse_markers is None:
        mouse_markers = ['o', 's', '^', 'D']

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
        colors = group_mouse_colors.get(group, group_mouse_colors['CTRL'])

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
                marker = mouse_markers[si % len(mouse_markers)]
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
            marker = mouse_markers[si % len(mouse_markers)]
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
        colors = group_mouse_colors[group]
        subjects = list(all_group_data.get(group, {}).keys())
        for i, subj in enumerate(subjects):
            short = subj.replace('sub-', '')[:15]
            legend_elements.append(plt.Line2D(
                [0], [0], color=colors[i % len(colors)], linewidth=2,
                marker=mouse_markers[i % len(mouse_markers)],
                label=f'{group} {short}'))
    legend_elements.append(plt.Line2D([0], [0], color='black', linewidth=1.3,
                                       linestyle='-', label='Dir 1 (solid)'))
    fig.legend(handles=legend_elements, loc='upper right', fontsize=7,
               bbox_to_anchor=(0.99, 0.98), title='Subjects', ncol=1)

    fig.suptitle(
        f'Fixel K-means Direction × Peak Amplitude — Bilateral (L+R) v11 [newmask]\n'
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
