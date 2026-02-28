"""Deep diagnostic: trace AP vs CTRL fixel data step by step."""
import sys
sys.path.insert(0, '/Users/glerma/soft/mbrain')
import yaml, numpy as np
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors, slicewise_kmeans,
    _load_mif_as_array
)

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
    fixel_dir = base / 'dwi' / 'fixels_group'

    print(f"\n{'='*70}")
    print(f"{grp} {subj}")
    print(f"{'='*70}")

    # 1. Load raw fixel data
    fd = load_fixel_dir(fixel_dir, 'mbrain')
    print(f"  Total fixels (whole brain): {len(fd['fd'])}")
    print(f"  peak_amp whole brain: min={fd['peak_amp'].min():.6f} max={fd['peak_amp'].max():.6f} mean={fd['peak_amp'].mean():.6f}")

    # 2. Load ROI
    roi_path = base / roi_rel
    if not roi_path.exists():
        print(f"  ROI NOT FOUND: {roi_path}")
        continue
    roi_mask = load_roi_mask(roi_path, 'mbrain')
    n_voxels = np.sum(roi_mask)
    print(f"  ROI voxels: {int(n_voxels)}")

    # 3. Extract fixels in ROI manually to see raw numbers
    index = fd['index']
    directions = fd['directions']
    peak_amp = fd['peak_amp']
    afd = fd['fd']

    # Count fixels per voxel in ROI
    fixels_in_roi = []
    peaks_in_roi = []
    afds_in_roi = []
    for x in range(index.shape[0]):
        for y in range(index.shape[1]):
            for z in range(index.shape[2]):
                if not roi_mask[x, y, z]:
                    continue
                n_fix = int(index[x, y, z, 0])
                offset = int(index[x, y, z, 1])
                for fi in range(n_fix):
                    idx = offset + fi
                    if idx < len(peak_amp):
                        fixels_in_roi.append(directions[idx])
                        peaks_in_roi.append(peak_amp[idx])
                        afds_in_roi.append(afd[idx])

    peaks_in_roi = np.array(peaks_in_roi)
    afds_in_roi = np.array(afds_in_roi)
    fixels_in_roi = np.array(fixels_in_roi)

    print(f"  Fixels in ROI: {len(peaks_in_roi)}")
    if len(peaks_in_roi) == 0:
        print("  NO FIXELS IN ROI!")
        continue
    print(f"  Fixels per ROI voxel: {len(peaks_in_roi)/n_voxels:.1f}")
    print(f"  peak_amp in ROI: min={peaks_in_roi.min():.6f} max={peaks_in_roi.max():.6f} mean={peaks_in_roi.mean():.6f} median={np.median(peaks_in_roi):.6f}")
    print(f"  afd in ROI:      min={afds_in_roi.min():.6f} max={afds_in_roi.max():.6f} mean={afds_in_roi.mean():.6f} median={np.median(afds_in_roi):.6f}")

    # Distribution of peak amplitudes
    thresholds = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
    print(f"  Peak amplitude distribution:")
    for t in thresholds:
        pct = 100 * np.sum(peaks_in_roi > t) / len(peaks_in_roi)
        print(f"    > {t:.2f}: {pct:.1f}%  ({np.sum(peaks_in_roi > t)} fixels)")

    # 4. Check direction norms (should be ~1 for real fixels)
    norms = np.linalg.norm(fixels_in_roi, axis=1)
    print(f"  Direction norms: min={norms.min():.4f} max={norms.max():.4f} mean={norms.mean():.4f}")

    # 5. Run slicewise extraction and kmeans
    sv, sp = extract_slicewise_vectors(fd, roi_mask, slice_axis=0)
    print(f"  Slices with fixels: {len(sv)}")
    skm = slicewise_kmeans(sv, n_clusters=2, slice_afds=sp)
    print(f"  Slices after kmeans: {len(skm)}")

    for s_idx in sorted(skm.keys())[:3]:
        r = skm[s_idx]
        sd = r['scaled_directions'][0]
        d = r['directions'][0]
        mp = r['mean_peaks'][0]
        print(f"    Slice {s_idx}: n_fixels={r['n_fixels']:4d}  mean_peak={mp:.4f}  "
              f"|scaled|={np.linalg.norm(sd):.4f}  "
              f"AP={abs(sd[0]):.4f} DV={abs(sd[1]):.4f} ML={abs(sd[2]):.4f}")
