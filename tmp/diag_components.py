"""Compare per-component peak-weighted averages across groups."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors, slicewise_kmeans
)
from pathlib import Path
import numpy as np
import yaml

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)

deriv = Path(cfg['paths']['derivatives_dir'])

subjects_check = {
    'CTRL': ('sub-129-921-R1F', 'ses-260115121231', 'atlas/rois/fod_space/roi_isocortex_left.nii.gz'),
    'AP':   ('sub-596-C-OO-NmM', 'ses-240911103349', 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'),
    'APTR': ('sub-AP72-862-901-NmF', 'ses-251210093342', 'manual_roi/fod_space/roi_Graft_left_flip_fod.mif'),
}

for grp, (subj, ses, roi_rel) in subjects_check.items():
    base = deriv / subj / ses
    fd = load_fixel_dir(base / 'dwi' / 'fixels_group', 'mbrain')
    roi = load_roi_mask(base / roi_rel, 'mbrain')
    sv, sp = extract_slicewise_vectors(fd, roi, slice_axis=0)
    skm = slicewise_kmeans(sv, n_clusters=2, slice_afds=sp)

    print(f"\n{'='*60}")
    print(f"{grp} {subj}")
    print(f"{'='*60}")
    for s_idx in sorted(skm.keys())[:3]:
        r = skm[s_idx]
        d = r['directions'][0]
        sd = r['scaled_directions'][0]
        mp = r['mean_peaks'][0]
        print(f"  Slice {s_idx}: unit_dir=[{d[0]:+.3f},{d[1]:+.3f},{d[2]:+.3f}]")
        print(f"            scaled  =[{sd[0]:+.3f},{sd[1]:+.3f},{sd[2]:+.3f}]  |scaled|={np.linalg.norm(sd):.4f}")
        print(f"            mean_peak={mp:.4f}, n_fixels={r['n_fixels']}, cluster_sizes={r['cluster_sizes']}")
        # Component comparison
        print(f"            AP_comp={abs(sd[0]):.4f}  DV_comp={abs(sd[1]):.4f}  ML_comp={abs(sd[2]):.4f}")
