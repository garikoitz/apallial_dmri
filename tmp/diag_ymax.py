"""Check the actual y_max values used in both plots."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
import yaml
import numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    slicewise_kmeans, analyze_3group_kmeans
)

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)

# Run the full pipeline for left, no shrink
result = analyze_3group_kmeans(
    cfg, hemisphere='left', use_shrink=False,
    output_dir=Path('/Users/glerma/soft/mbrain/reports'),
    n_clusters=2, conda_env='mbrain', force=False, show=False, verbose=False)

all_group_data = result['all_group_data']

# Replicate the global_max calculation
print("\n=== Per-group max |scaled_direction| component ===")
for group in ['CTRL', 'AP', 'APTR']:
    gmax = 0.0
    for subj, data in all_group_data.get(group, {}).items():
        for s_idx, res in data.items():
            for sd in res.get('scaled_directions', []):
                val = np.max(np.abs(sd))
                if val > gmax:
                    gmax = val
    print(f"  {group}: max |component| = {gmax:.4f}")

global_max = max(
    np.max(np.abs(sd))
    for group in ['CTRL', 'AP', 'APTR']
    for subj, data in all_group_data.get(group, {}).items()
    for s_idx, res in data.items()
    for sd in res.get('scaled_directions', [])
)
y_max_radar = max(0.05, np.ceil(global_max * 20) / 20)
y_max_lines = max(0.05, global_max * 1.15)
print(f"\n  Global max: {global_max:.4f}")
print(f"  Radar y_max: {y_max_radar:.4f}")
print(f"  Lines y_max: {y_max_lines:.4f}")

# Also check: what does a typical CTRL vs AP slice look like?
print("\n=== Sample scaled_directions (Dir1 only) ===")
for group in ['CTRL', 'AP', 'APTR']:
    subj = list(all_group_data[group].keys())[0]
    data = all_group_data[group][subj]
    first_s = sorted(data.keys())[0]
    sd = data[first_s]['scaled_directions'][0]
    print(f"  {group} {subj} slice {first_s}: |AP|={abs(sd[0]):.4f} |DV|={abs(sd[1]):.4f} |ML|={abs(sd[2]):.4f}  |total|={np.linalg.norm(sd):.4f}")
