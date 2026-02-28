"""Test voxel-normalized scaling."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors, slicewise_kmeans
)
from pathlib import Path
import numpy as np, yaml

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)
deriv = Path(cfg['paths']['derivatives_dir'])

checks = {
    'CTRL': ('sub-129-921-R1F', 'ses-260115121231',
             'atlas/rois/fod_space/roi_isocortex_left.nii.gz'),
    'AP':   ('sub-596-C-OO-NmM', 'ses-240911103349',
             'manual_roi/fod_space/shrink/roi_Graft_left_flip_fod_shrink.mif'),
    'APTR': ('sub-AP72-862-901-NmF', 'ses-251210093342',
             'manual_roi/fod_space/shrink/roi_Graft_left_flip_fod_shrink.mif'),
}

for grp, (subj, ses, roi_rel) in checks.items():
    base = deriv / subj / ses
    fd = load_fixel_dir(base / 'dwi' / 'fixels_group', 'mbrain')
    roi = load_roi_mask(base / roi_rel, 'mbrain')
    sv, sp, svox = extract_slicewise_vectors(fd, roi, slice_axis=0)
    skm = slicewise_kmeans(sv, n_clusters=2, slice_afds=sp, slice_roi_voxels=svox)
    print(f'\n{grp:5s} {subj}:')
    for s_idx in sorted(skm.keys())[:4]:
        r = skm[s_idx]
        sd = r['scaled_directions'][0]
        print(f'  Slice {s_idx}: voxels={r["n_roi_voxels"]:4d} fixels={r["n_fixels"]:4d}  '
              f'AP={abs(sd[0]):.5f} DV={abs(sd[1]):.5f} ML={abs(sd[2]):.5f}  '
              f'|sd|={np.linalg.norm(sd):.5f}  '
              f'fix/vox={r["n_fixels"]/r["n_roi_voxels"]:.2f}')
