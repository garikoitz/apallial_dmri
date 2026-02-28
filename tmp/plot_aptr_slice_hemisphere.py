"""Pick a random APTR graft slice and generate the 3D hemisphere + radar figure."""
import sys, yaml, random
import numpy as np
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    cluster_fixel_directions, plot_fixel_directions,
)

# Load config
with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)
deriv = Path(cfg['paths']['derivatives_dir'])

# Pick a random APTR subject
aptr_subjects = [(s['id'], s['session']) for s in cfg['groups']['APTR']]
random.seed(42)
subj_id, ses = random.choice(aptr_subjects)
print(f'Selected APTR subject: {subj_id} / {ses}')

base = deriv / subj_id / ses

# Load fixels and graft ROI (left hemisphere, shrink)
fixel_dir = base / 'dwi' / 'fixels_group'
roi_path = base / 'manual_roi' / 'fod_space' / 'shrink' / 'roi_Graft_left_flip_fod_shrink.mif'

print(f'Fixel dir: {fixel_dir}')
print(f'ROI: {roi_path}')

fd = load_fixel_dir(fixel_dir, 'mbrain')
roi = load_roi_mask(roi_path, 'mbrain')

# Extract per-slice vectors
sv, sp, svox = extract_slicewise_vectors(fd, roi, slice_axis=0)
print(f'\nAvailable slices: {sorted(sv.keys())}')
for s_idx in sorted(sv.keys()):
    print(f'  Slice {s_idx}: {len(sv[s_idx])} fixels, {svox[s_idx]} ROI voxels')

# Pick a slice with a good number of fixels (middle-ish)
slice_indices = sorted(sv.keys())
mid_idx = len(slice_indices) // 2
chosen_slice = slice_indices[mid_idx]
print(f'\nChosen slice: {chosen_slice} ({len(sv[chosen_slice])} fixels)')

# Get vectors and peaks for the chosen slice
all_vectors = sv[chosen_slice]
slice_peaks = sp[chosen_slice]
# Create dummy peak indices (all from "peak 1" since we loaded from directions.mif)
all_peak_indices = np.ones(len(all_vectors), dtype=int)
n_peaks = 1

# Cluster
result = cluster_fixel_directions(
    all_vectors, all_peak_indices, n_peaks,
    n_clusters=2, verbose=True,
)

# Plot
out_path = Path('/Users/glerma/soft/mbrain/reports') / f'fixel_hemisphere_{subj_id}_slice{chosen_slice}.png'
fig = plot_fixel_directions(
    result,
    subject_id=f'{subj_id} (graft left, slice {chosen_slice})',
    save_path=out_path,
    show=False,
    peaks=slice_peaks,
)
print(f'\nSaved: {out_path}')
