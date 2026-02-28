"""Check y_max across all 4 hemisphere/shrink combinations."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
import yaml, numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import analyze_3group_kmeans

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)

out = Path('/Users/glerma/soft/mbrain/reports')

for hemi in ['left', 'right']:
    for shrink in [False, True]:
        tag = f'{hemi}{"_shrink" if shrink else ""}'
        r = analyze_3group_kmeans(cfg, hemisphere=hemi, use_shrink=shrink,
                                  output_dir=out, n_clusters=2, conda_env='mbrain',
                                  force=False, show=False, verbose=False)
        if not r:
            print(f'{tag}: NO RESULT')
            continue
        gmax = 0
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, data in r['all_group_data'].get(group, {}).items():
                for s, res in data.items():
                    for sd in res.get('scaled_directions', []):
                        val = np.max(np.abs(sd))
                        if val > gmax:
                            gmax = val
        print(f'{tag}: global_max={gmax:.4f}  y_max_radar={max(0.05, np.ceil(gmax*20)/20):.4f}')
