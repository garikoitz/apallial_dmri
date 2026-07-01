from __future__ import annotations

import numpy as np


def draw_single_seg_panel(
    ax,
    data_df,
    seg: str,
    color: str,
    title: str,
    *,
    targets_order: list[str],
    target_labels: dict[str, str],
    plot_subjects,
    marker_list: list[str],
):
    """Draw one APTR segment panel for the density bar plot."""
    seg_data = data_df[data_df['segment'] == seg]
    means, sems, all_vals = [], [], []
    for tgt in targets_order:
        vals = seg_data[seg_data['target'] == tgt]['sift2_density'].values
        means.append(np.mean(vals) if len(vals) > 0 else 0)
        sems.append(np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        all_vals.append(vals)

    n_tgt = len(targets_order)
    x_c = np.arange(n_tgt)
    ax.bar(x_c, means, 0.55, yerr=sems, capsize=4,
           color=color, alpha=0.55, edgecolor='black', linewidth=0.5, zorder=2)

    for tgt_idx in range(n_tgt):
        tgt_subjs = seg_data[seg_data['target'] == targets_order[tgt_idx]]['subject'].values
        for val_idx, val in enumerate(all_vals[tgt_idx]):
            subj_idx = list(plot_subjects).index(tgt_subjs[val_idx])
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