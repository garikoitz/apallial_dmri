"""Quick test: run analyze_3group_kmeans for one combo and check outputs."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
import yaml
from pathlib import Path
from functions.fixel_roi_analysis import analyze_3group_kmeans

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)

out_dir = Path('/Users/glerma/soft/mbrain/reports')

result = analyze_3group_kmeans(
    cfg,
    hemisphere='left',
    use_shrink=False,
    output_dir=out_dir,
    n_clusters=2,
    conda_env='mbrain',
    force=False,
    show=False,
    verbose=True,
)

if result:
    print(f"\n=== RESULT ===")
    print(f"Radar: {result.get('figure_path', 'MISSING')}")
    print(f"Lines: {result.get('lines_path', 'MISSING')}")
    radar_exists = Path(result['figure_path']).exists() if result.get('figure_path') else False
    lines_exists = Path(result['lines_path']).exists() if result.get('lines_path') else False
    print(f"Radar file exists: {radar_exists}")
    print(f"Lines file exists: {lines_exists}")
else:
    print("No result!")
