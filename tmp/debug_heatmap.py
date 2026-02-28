"""Debug: show whiskers AND heatmap together to compare."""
import sys, yaml, random
import numpy as np
sys.path.insert(0, '/Users/glerma/soft/mbrain')
from pathlib import Path
from functions.fixel_roi_analysis import (
    load_fixel_dir, load_roi_mask, extract_slicewise_vectors,
    cluster_fixel_directions,
)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import to_rgba

with open('/Users/glerma/soft/mbrain/config1.yaml') as f:
    cfg = yaml.safe_load(f)
deriv = Path(cfg['paths']['derivatives_dir'])

# Same APTR subject/slice as before
subj_id, ses = 'sub-AP72-862-901-NmF', 'ses-251210093342'
base = deriv / subj_id / ses
fd = load_fixel_dir(base / 'dwi' / 'fixels_group', 'mbrain')
roi = load_roi_mask(base / 'manual_roi' / 'fod_space' / 'shrink' / 'roi_Graft_left_flip_fod_shrink.mif', 'mbrain')
sv, sp, svox = extract_slicewise_vectors(fd, roi, slice_axis=0)
chosen = sorted(sv.keys())[len(sv) // 2]
vecs = sv[chosen].copy()
peaks = sp[chosen].copy()

# Flip to DV >= 0
flip = vecs[:, 1] < 0
vecs[flip] *= -1

print(f'Slice {chosen}: {len(vecs)} fixels')
print(f'  AP range: [{vecs[:,0].min():.3f}, {vecs[:,0].max():.3f}]')
print(f'  DV range: [{vecs[:,1].min():.3f}, {vecs[:,1].max():.3f}]')
print(f'  ML range: [{vecs[:,2].min():.3f}, {vecs[:,2].max():.3f}]')
print(f'  Peak range: [{peaks.min():.3f}, {peaks.max():.3f}]')

# Spherical coords
elev = np.arcsin(np.clip(vecs[:, 1], -1, 1))
az = np.arctan2(vecs[:, 2], vecs[:, 0]) % (2 * np.pi)
print(f'  Azimuth range: [{np.degrees(az.min()):.1f}, {np.degrees(az.max()):.1f}] deg')
print(f'  Elevation range: [{np.degrees(elev.min()):.1f}, {np.degrees(elev.max()):.1f}] deg')

# Bin into grid
n_u, n_v = 40, 20
u_edges = np.linspace(0, 2 * np.pi, n_u + 1)
v_edges = np.linspace(0, np.pi / 2, n_v + 1)
peak_grid = np.zeros((n_u, n_v))
count_grid = np.zeros((n_u, n_v), dtype=int)
for fi in range(len(vecs)):
    ui = np.searchsorted(u_edges, az[fi], side='right') - 1
    vi = np.searchsorted(v_edges, elev[fi], side='right') - 1
    ui = np.clip(ui, 0, n_u - 1)
    vi = np.clip(vi, 0, n_v - 1)
    peak_grid[ui, vi] += peaks[fi]
    count_grid[ui, vi] += 1

n_nonempty = np.sum(peak_grid > 0)
print(f'\n  Non-empty quads: {n_nonempty} / {n_u * n_v}')
print(f'  Max peak sum in a quad: {peak_grid.max():.3f}')

# Show whiskers + heatmap side by side
fig = plt.figure(figsize=(20, 8))

# Left: whiskers
ax1 = fig.add_subplot(131, projection='3d')
ax1.set_title('Whiskers (reference)')
pk_max = peaks.max()
scale = 0.4 / pk_max
rng = np.random.RandomState(42)
idx = rng.choice(len(vecs), min(400, len(vecs)), replace=False)
for i in idx:
    v = vecs[i]
    l = peaks[i] * scale
    ax1.plot([v[0], v[0]*(1+l)], [v[2], v[2]*(1+l)], [v[1], v[1]*(1+l)],
             color='#E24A33', alpha=0.4, linewidth=0.8)
ax1.scatter(vecs[idx, 0], vecs[idx, 2], vecs[idx, 1], c='#348ABD', alpha=0.3, s=4)
ax1.set_xlabel('AP'); ax1.set_ylabel('ML'); ax1.set_zlabel('DV')
ax1.set_xlim(-1.1, 1.1); ax1.set_ylim(-1.1, 1.1); ax1.set_zlim(-0.1, 1.1)
ax1.view_init(elev=25, azim=-60)

# Middle: heatmap (only non-empty quads)
ax2 = fig.add_subplot(132, projection='3d')
ax2.set_title('Heatmap (non-empty only)')
peak_max = peak_grid.max()
quads = []
fcolors = []
for i in range(n_u):
    for j in range(n_v):
        if peak_grid[i, j] <= 0:
            continue
        corners = []
        for a, e in [(u_edges[i], v_edges[j]), (u_edges[i+1], v_edges[j]),
                      (u_edges[i+1], v_edges[j+1]), (u_edges[i], v_edges[j+1])]:
            cx = np.cos(a) * np.cos(e)
            cy = np.sin(e)
            cz = np.sin(a) * np.cos(e)
            corners.append((cx, cz, cy))  # X=AP, Y=ML, Z=DV
        quads.append(corners)
        rel = peak_grid[i, j] / peak_max
        fcolors.append(to_rgba('#E24A33', alpha=max(0.1, rel)))

poly = Poly3DCollection(quads, linewidths=0.3, edgecolors=(0.3, 0.3, 0.3, 0.3))
poly.set_facecolors(fcolors)
ax2.add_collection3d(poly)
ax2.set_xlabel('AP'); ax2.set_ylabel('ML'); ax2.set_zlabel('DV')
ax2.set_xlim(-1.1, 1.1); ax2.set_ylim(-1.1, 1.1); ax2.set_zlim(-0.1, 1.1)
ax2.view_init(elev=25, azim=-60)

# Right: both overlaid
ax3 = fig.add_subplot(133, projection='3d')
ax3.set_title('Both overlaid')
poly2 = Poly3DCollection(quads, linewidths=0.3, edgecolors=(0.3, 0.3, 0.3, 0.2))
poly2.set_facecolors(fcolors)
ax3.add_collection3d(poly2)
for i in idx:
    v = vecs[i]
    l = peaks[i] * scale
    ax3.plot([v[0], v[0]*(1+l)], [v[2], v[2]*(1+l)], [v[1], v[1]*(1+l)],
             color='blue', alpha=0.4, linewidth=0.8)
ax3.set_xlabel('AP'); ax3.set_ylabel('ML'); ax3.set_zlabel('DV')
ax3.set_xlim(-1.1, 1.1); ax3.set_ylim(-1.1, 1.1); ax3.set_zlim(-0.1, 1.1)
ax3.view_init(elev=25, azim=-60)

plt.tight_layout()
out = '/Users/glerma/soft/mbrain/reports/debug_heatmap_vs_whiskers.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f'\nSaved: {out}')
