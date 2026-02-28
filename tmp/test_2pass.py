"""Test 2-pass approach with unified y_max."""
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
combos = [('left', False), ('left', True), ('right', False), ('right', True)]

# Pass 1: data only
print('Pass 1: computing data...')
all_results = {}
global_max = 0.0
for hemi, shrink in combos:
    tag = f'{hemi}{"_shrink" if shrink else ""}'
    r = analyze_3group_kmeans(cfg, hemisphere=hemi, use_shrink=shrink,
                               output_dir=out, n_clusters=2, conda_env='mbrain',
                               force=False, show=False, verbose=False, plot=False)
    all_results[tag] = (hemi, shrink, r)
    if r:
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, data in r['all_group_data'].get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get('scaled_directions', []):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val

y_max = max(0.05, np.ceil(global_max * 20) / 20)
print(f'  Unified y_max = {y_max:.4f} (global max = {global_max:.4f})\n')

# Pass 2: plot with forced y_max
for tag, (hemi, shrink, r) in all_results.items():
    if not r:
        print(f'{tag}: no data')
        continue
    suffix = '_shrink' if shrink else ''
    tns = r['target_n_slices']
    fig_path = out / f'fixel_3group_comparison_{hemi}{suffix}.png'
    lines_path = out / f'fixel_3group_components_{hemi}{suffix}.png'

    plot_3group_kmeans_comparison(
        r['all_group_data'], hemi, tns,
        use_shrink=shrink, save_path=fig_path, show=False, y_max=y_max)
    plot_3group_component_lines(
        r['all_group_data'], hemi, tns,
        use_shrink=shrink, save_path=lines_path, show=False, y_max=y_max)
    print(f'{tag}: radar={fig_path.name} lines={lines_path.name}')

print('\nDone! All 8 plots use y_max =', y_max)
