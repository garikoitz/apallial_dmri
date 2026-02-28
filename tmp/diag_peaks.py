"""Diagnostic: check actual peak amplitude values across groups."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from functions.fixel_roi_analysis import load_fixel_dir, load_roi_mask, extract_slicewise_vectors
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

    fixel_dir = base / 'dwi' / 'fixels_group'
    fd = load_fixel_dir(fixel_dir, 'mbrain')
    print(f"\n{grp} {subj}")
    print(f"  Total fixels: {len(fd['fd'])}")
    print(f"  peak_amp range: {fd['peak_amp'].min():.6f} - {fd['peak_amp'].max():.6f}, mean={fd['peak_amp'].mean():.6f}")
    print(f"  fd range: {fd['fd'].min():.6f} - {fd['fd'].max():.6f}, mean={fd['fd'].mean():.6f}")

    roi_path = base / roi_rel
    if not roi_path.exists():
        print(f"  ROI not found: {roi_path}")
        continue

    roi_mask = load_roi_mask(roi_path, 'mbrain')
    sv, sp = extract_slicewise_vectors(fd, roi_mask, slice_axis=0)
    if sp:
        all_peaks = np.concatenate(list(sp.values()))
        all_vecs = np.concatenate(list(sv.values()))
        print(f"  Slices with fixels: {len(sv)}")
        print(f"  peak_amp in ROI: range={all_peaks.min():.6f}-{all_peaks.max():.6f}, mean={all_peaks.mean():.6f}")
        print(f"  Total fixels in ROI: {len(all_peaks)}")

        # Show what scaled_directions look like
        from functions.fixel_roi_analysis import slicewise_kmeans
        skm = slicewise_kmeans(sv, n_clusters=2, slice_afds=sp)
        if skm:
            first_s = sorted(skm.keys())[0]
            print(f"  First slice ({first_s}): dir1={skm[first_s]['directions'][0]}, "
                  f"mean_peak={skm[first_s]['mean_peaks'][0]:.6f}, "
                  f"scaled={skm[first_s]['scaled_directions'][0]}")
    else:
        print("  No slices with fixels in ROI!")
