"""Regenerate shrink-only plots with unified y_max."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
import yaml, numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import (
    analyze_3group_kmeans, plot_3group_kmeans_comparison, plot_3group_component_lines
)

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)

out = Path('/Users/glerma/soft/mbrain/reports')

# Pass 1: data only, shrink only
print('Pass 1: computing data (shrink ROIs only)...')
all_results = {}
global_max = 0.0
for hemi in ['left', 'right']:
    tag = f'{hemi}_shrink'
    r = analyze_3group_kmeans(cfg, hemisphere=hemi, use_shrink=True,
                               output_dir=out, n_clusters=2, conda_env='mbrain',
                               force=False, show=False, verbose=True, plot=False)
    all_results[tag] = r
    if r:
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, data in r['all_group_data'].get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get('scaled_directions', []):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val

y_max = max(0.05, np.ceil(global_max * 20) / 20)
print(f'\n  Unified y_max = {y_max:.4f} (global max = {global_max:.4f})\n')

# Pass 2: plot with forced y_max
for hemi in ['left', 'right']:
    tag = f'{hemi}_shrink'
    r = all_results[tag]
    if not r:
        print(f'{tag}: no data'); continue
    tns = r['target_n_slices']

    fig_path = out / f'fixel_3group_comparison_{hemi}_shrink.png'
    plot_3group_kmeans_comparison(
        r['all_group_data'], hemi, tns,
        use_shrink=True, save_path=fig_path, show=False, y_max=y_max)

    lines_path = out / f'fixel_3group_components_{hemi}_shrink.png'
    plot_3group_component_lines(
        r['all_group_data'], hemi, tns,
        use_shrink=True, save_path=lines_path, show=False, y_max=y_max)

    print(f'{tag}: done')

print(f'\nAll plots use y_max = {y_max}')
