#!/usr/bin/env python3
"""Test bar plot for APTR density — 4 subplots: combined + per-segment."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

csv_path = Path('reports/tract_metrics_v05_cleaned.csv')
df = pd.read_csv(csv_path)

# Filter APTR, left hemisphere only
aptr_left = df[(df['group'] == 'APTR') & (df['hemisphere'] == 'left')].copy()
targets_order = ['olfactory_cortex', 'striatum', 'thalamus', 'cerebellum']
target_labels = {'thalamus': 'Thalamus', 'striatum': 'Striatum',
                 'olfactory_cortex': 'Paleocortex', 'cerebellum': 'Cerebellum'}
segments_order = ['anterior', 'central', 'posterior']
seg_colors = {'anterior': '#E24A33', 'central': '#348ABD', 'posterior': '#988ED5'}

subjects = aptr_left['subject'].unique()
n_subj = len(subjects)
marker_list = ['o', 's', 'D', '^']
n_tgt = len(targets_order)


def _bar_panel(ax, data_df, seg_list, title, show_legend=False):
    n_seg = len(seg_list)
    bar_width = 0.22
    group_width = n_seg * bar_width + 0.15
    x_tgt = np.arange(n_tgt) * group_width

    for s_idx, seg in enumerate(seg_list):
        seg_data = data_df[data_df['segment'] == seg]
        means, sems = [], []
        for tgt in targets_order:
            vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
            means.append(np.mean(vals) if len(vals) > 0 else 0)
            sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        x_pos = x_tgt + s_idx * bar_width
        ax.bar(x_pos, means, bar_width * 0.85, yerr=sems, capsize=4,
               color=seg_colors[seg], alpha=0.55, edgecolor='black', linewidth=0.5,
               label=seg.capitalize(), zorder=2)
        for tgt_idx, tgt in enumerate(targets_order):
            tgt_rows = seg_data[seg_data['target'] == tgt]
            for _, row in tgt_rows.iterrows():
                subj_idx = list(subjects).index(row['subject'])
                jitter = (np.random.RandomState(subj_idx).rand() - 0.5) * bar_width * 0.5
                ax.scatter(x_pos[tgt_idx] + jitter, row['sift2_density'],
                           marker=marker_list[subj_idx % len(marker_list)],
                           color=seg_colors[seg], edgecolor='black', s=40, zorder=3,
                           linewidth=0.6)
    ax.set_xticks(x_tgt + bar_width * (n_seg - 1) / 2)
    ax.set_xticklabels([target_labels[t] for t in targets_order], fontsize=10)
    ax.set_title(title, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    if show_legend:
        ax.legend(title='Segment', fontsize=9)


def _single_seg_panel(ax, data_df, seg, color, title):
    seg_data = data_df[data_df['segment'] == seg]
    means, sems, all_vals = [], [], []
    for tgt in targets_order:
        vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
        means.append(np.mean(vals) if len(vals) > 0 else 0)
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        all_vals.append(vals)
    x_c = np.arange(n_tgt)
    ax.bar(x_c, means, 0.55, yerr=sems, capsize=4,
           color=color, alpha=0.55, edgecolor='black', linewidth=0.5, zorder=2)
    for tgt_idx in range(n_tgt):
        tgt_subjs = seg_data[seg_data['target'] == targets_order[tgt_idx]]['subject'].values
        for val_idx, val in enumerate(all_vals[tgt_idx]):
            subj_idx = list(subjects).index(tgt_subjs[val_idx])
            short = tgt_subjs[val_idx].replace('sub-', '')
            jitter = (np.random.RandomState(subj_idx).rand() - 0.5) * 0.3
            label = short if tgt_idx == 0 else None
            ax.scatter(x_c[tgt_idx] + jitter, val,
                       marker=marker_list[subj_idx % len(marker_list)],
                       color=color, edgecolor='black', s=50, zorder=3,
                       linewidth=0.6, label=label)
    ax.set_xticks(x_c)
    ax.set_xticklabels([target_labels[t] for t in targets_order], fontsize=9,
                       rotation=30, ha='right')
    ax.set_title(title, fontsize=12)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(title='Subject', fontsize=6, loc='upper right')


import matplotlib.gridspec as gridspec
from matplotlib.image import imread as mpl_imread

fig = plt.figure(figsize=(20, 14))
fig.suptitle('APTR — SIFT2 Density (Left Hemisphere, Cleaned Tracts)',
             fontsize=14, fontweight='bold')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3,
                       height_ratios=[1, 0.9])

# Row 0: 3 segment bar subplots
ax_ant = fig.add_subplot(gs[0, 0])
_single_seg_panel(ax_ant, aptr_left, 'anterior', seg_colors['anterior'],
                  'Anterior Graft → Destinations')
ax_ant.set_ylabel('SIFT2 Density\n(conn / seed_vox × 1000)', fontsize=10)

ax_cen = fig.add_subplot(gs[0, 1])
_single_seg_panel(ax_cen, aptr_left, 'central', seg_colors['central'],
                  'Central Graft → Destinations')
ax_cen.set_ylabel('SIFT2 Density', fontsize=10)

ax_pos = fig.add_subplot(gs[0, 2])
_single_seg_panel(ax_pos, aptr_left, 'posterior', seg_colors['posterior'],
                  'Posterior Graft → Destinations')
ax_pos.set_ylabel('SIFT2 Density', fontsize=10)

# Row 1: mrview screenshots by destination for sub-Apallial20R1F
# Save screenshots from mrview as:
#   reports/mrview_most_paleocortex.png  reports/mrview_most_striatum.png  reports/mrview_most_thalamus.png
destinations = ['paleocortex', 'striatum', 'thalamus']
dest_labels  = {'paleocortex': 'Paleocortex', 'striatum': 'Striatum', 'thalamus': 'Thalamus'}

for col_idx, dest in enumerate(destinations):
    ax = fig.add_subplot(gs[1, col_idx])
    img_path = Path(f'reports/mrview_most_{dest}.png')
    title = f'{dest_labels[dest]}\nsub-Apallial20R1F'
    if img_path.exists():
        img = mpl_imread(str(img_path))
        ax.imshow(img)
        ax.set_title(title, fontsize=10, fontweight='bold')
    else:
        ax.text(0.5, 0.5, f'Screenshot not yet saved:\n{img_path.name}',
                ha='center', va='center', fontsize=11, color='gray',
                transform=ax.transAxes)
        ax.set_title(title, fontsize=10, fontweight='bold', color='gray')
    ax.set_xticks([])
    ax.set_yticks([])

save_path = Path('reports/barplot_aptr_density_left.png')
fig.savefig(save_path, dpi=150, bbox_inches='tight')
print(f'Saved: {save_path}')
print(f'APTR subjects: {list(subjects)}')
print(f'N rows APTR left: {len(aptr_left)}')
