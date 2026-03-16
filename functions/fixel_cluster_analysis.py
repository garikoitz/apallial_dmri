"""
Fixel Cluster Analysis — Spatial + Directional 3D Clustering
=============================================================

Identify "structures" within a graft (or any ROI) by clustering voxels that
share similar dominant fixel directions **and** are spatially close.

Algorithm
---------
1. Load bilateral graft ROI mask (left | right)
2. Load fixel data (index, directions, peak amplitude)
3. For each ROI voxel: extract the **dominant fixel** (highest peak amplitude)
4. Build feature matrix: [x_norm, y_norm, z_norm, |dx|, |dy|, |dz|]
   - Spatial coords normalised to [0, 1]
   - Absolute direction values (antipodal symmetry: v ≡ −v)
5. K-means silhouette scan k=2..8, pick best k
6. Label each cluster by its principal direction (DV / ML / AP)
7. Write outputs to ``fixel_clusters/`` inside subject dwi dir

Outputs
-------
- ``cluster_labels.mif``      — 3D integer volume for mrview overlay
- ``cluster_info.json``        — per-cluster: id, n_voxels, main_dir, mean_vec
- ``silhouette_scores.png``    — elbow/silhouette plot

Usage
-----
    python functions/fixel_cluster_analysis.py          # runs on AP72-NmF
    python functions/fixel_cluster_analysis.py --all     # (future) all subjects
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.preprocessing import MinMaxScaler
from scipy.ndimage import label as ndimage_label
from scipy.sparse import lil_matrix

# ---- project paths --------------------------------------------------------
DERIV = Path('/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask')

# Subject catalogue (can grow; for now just APTR AP72-NmF)
SUBJECT_MAP = {
    'sub-AP72-862-901-NmF': {
        'ses': 'ses-251210093342',
        'group': 'APTR',
        'roi_type': 'graft',   # AP/APTR → graft;  CTRL → isocortex
    },
    'sub-AP72-862-901-R1M': {
        'ses': 'ses-251210100043',
        'group': 'APTR',
        'roi_type': 'graft',
    },
    'sub-Apallial17NmF': {
        'ses': 'ses-240816091633',
        'group': 'APTR',
        'roi_type': 'graft',
    },
    'sub-Apallial20R1F': {
        'ses': 'ses-240909094232',
        'group': 'APTR',
        'roi_type': 'graft',
    },
}

# K-range for silhouette scan
K_MIN, K_MAX = 2, 20

# Direction labels (axis-index convention after flip: 0=AP, 1=DV, 2=ML)
DIR_LABELS = {0: 'AP', 1: 'DV', 2: 'ML'}

# Balancing weight between spatial and direction features.
# alpha > 1 favours spatial proximity → more spatially contiguous clusters.
SPATIAL_WEIGHT = 2.5

# Minimum cluster size (voxels).  Tiny fragments get merged into nearest cluster.
MIN_CLUSTER_SIZE = 20

# Desired k range — the automatic selector is clamped to this window.
K_CLAMP_MIN, K_CLAMP_MAX = 6, 15

# Use spatially-constrained agglomerative clustering by default
USE_AGGLOMERATIVE = True

CONDA_ENV = 'mbrain'   # for micromamba fallback


# ===========================================================================
#  Low-level helpers (mirror fixel_roi_analysis.py for consistency)
# ===========================================================================

def _load_mif_as_array(mif_path: str | Path) -> np.ndarray:
    """Convert .mif → temp NIfTI and return numpy array."""
    mif_path = Path(mif_path)
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        cmd_base = f'mrconvert "{mif_path}" "{tmp_path}" -force -quiet'
        if shutil.which('mrconvert'):
            subprocess.run(cmd_base, shell=True, check=True, capture_output=True)
        else:
            subprocess.run(
                f'micromamba run -n {CONDA_ENV} {cmd_base}',
                shell=True, check=True, capture_output=True,
            )
        data = nib.load(tmp_path).get_fdata()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return data


def _save_array_as_mif(array: np.ndarray, ref_mif: str | Path,
                        out_mif: str | Path) -> None:
    """Save *array* as .mif, stealing the header/transform from *ref_mif*.

    Workflow: array → temp NIfTI → mrconvert with -template ref → out .mif.
    We use the *index.mif* spatial grid (first 3 dims) as reference so that
    the cluster volume is in the same voxel space.
    """
    ref_mif = Path(ref_mif)
    out_mif = Path(out_mif)
    out_mif.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_nii = tmp.name
    try:
        # Load reference NIfTI header (we need its affine)
        with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp2:
            ref_nii = tmp2.name
        try:
            cmd_ref = f'mrconvert "{ref_mif}" "{ref_nii}" -force -quiet'
            if shutil.which('mrconvert'):
                subprocess.run(cmd_ref, shell=True, check=True, capture_output=True)
            else:
                subprocess.run(
                    f'micromamba run -n {CONDA_ENV} {cmd_ref}',
                    shell=True, check=True, capture_output=True,
                )
            ref_img = nib.load(ref_nii)
            affine = ref_img.affine
        finally:
            if os.path.exists(ref_nii):
                os.remove(ref_nii)

        # Write array with same affine
        img = nib.Nifti1Image(array.astype(np.float32), affine)
        nib.save(img, tmp_nii)

        # Convert back to .mif so mrview can overlay
        cmd_out = f'mrconvert "{tmp_nii}" "{out_mif}" -force -quiet'
        if shutil.which('mrconvert'):
            subprocess.run(cmd_out, shell=True, check=True, capture_output=True)
        else:
            subprocess.run(
                f'micromamba run -n {CONDA_ENV} {cmd_out}',
                shell=True, check=True, capture_output=True,
            )
    finally:
        if os.path.exists(tmp_nii):
            os.remove(tmp_nii)


def load_fixel_dir(fixel_dir: str | Path) -> dict:
    """Load fixel data from fod2fixel output directory.

    Returns dict: index (X,Y,Z,2), directions (N,3), fd (N,), peak_amp (N,).
    """
    fixel_dir = Path(fixel_dir)
    index = _load_mif_as_array(fixel_dir / 'index.mif')
    directions_raw = _load_mif_as_array(fixel_dir / 'directions.mif')

    afd_path = fixel_dir / 'afd.mif'
    fd_path = fixel_dir / 'fd.mif'
    density_path = afd_path if afd_path.exists() else fd_path
    fd = _load_mif_as_array(density_path).flatten()

    peak_path = fixel_dir / 'peak.mif'
    if not peak_path.exists():
        peak_path = fixel_dir / 'peak_amp.mif'
    peak_amp = _load_mif_as_array(peak_path).flatten() if peak_path.exists() else fd.copy()

    dirs = directions_raw[:, :, 0] if directions_raw.ndim == 3 else directions_raw
    return dict(index=index, directions=dirs, fd=fd, peak_amp=peak_amp)


def load_roi_mask(roi_path: str | Path) -> np.ndarray:
    """Load ROI mask (.mif or .nii.gz) as boolean 3D array."""
    roi_path = Path(roi_path)
    if roi_path.suffix == '.mif':
        return _load_mif_as_array(roi_path) > 0.5
    else:
        return nib.load(str(roi_path)).get_fdata() > 0.5


# ===========================================================================
#  Core analysis
# ===========================================================================

def extract_dominant_fixels(fixel_data: dict, roi_mask: np.ndarray,
                             min_fd: float = 0.001):
    """For every ROI voxel, pick the fixel with the highest peak amplitude.

    Returns
    -------
    coords : np.ndarray (M, 3)   — voxel coordinates (ijk)
    directions : np.ndarray (M, 3) — dominant fixel unit-direction
    peaks : np.ndarray (M,)       — peak amplitude of that fixel
    """
    index = fixel_data['index']
    dirs = fixel_data['directions']
    fd = fixel_data['fd']
    peak_amp = fixel_data['peak_amp']

    xs, ys, zs = np.where(roi_mask)
    coord_list = []
    dir_list = []
    peak_list = []

    for x, y, z in zip(xs, ys, zs):
        count = int(index[x, y, z, 0])
        offset = int(index[x, y, z, 1])
        if count == 0:
            continue

        best_peak = -1.0
        best_dir = None
        for i in range(count):
            fi = offset + i
            if fi >= len(fd):
                break
            if fd[fi] < min_fd:
                continue
            if peak_amp[fi] > best_peak:
                best_peak = peak_amp[fi]
                d = dirs[fi]
                norm = np.linalg.norm(d)
                if norm > 0.01:
                    best_dir = d / norm
                    best_peak = peak_amp[fi]

        if best_dir is not None:
            coord_list.append([x, y, z])
            dir_list.append(best_dir)
            peak_list.append(best_peak)

    return (np.array(coord_list),
            np.array(dir_list),
            np.array(peak_list))


def build_feature_matrix(coords: np.ndarray, directions: np.ndarray,
                          spatial_weight: float = SPATIAL_WEIGHT) -> np.ndarray:
    """Combine normalised spatial coords + |direction| into a feature matrix.

    Parameters
    ----------
    coords : (M, 3)  — raw voxel coordinates
    directions : (M, 3) — unit direction vectors
    spatial_weight : float — scales spatial columns relative to direction

    Returns
    -------
    features : (M, 6)
    """
    # Normalise spatial coords to [0, 1]
    scaler = MinMaxScaler()
    coords_norm = scaler.fit_transform(coords.astype(float))

    # Absolute direction (antipodal symmetry)
    dirs_abs = np.abs(directions)

    # Weight spatial block
    features = np.hstack([coords_norm * spatial_weight, dirs_abs])
    return features


def _build_connectivity_graph(coords: np.ndarray) -> lil_matrix:
    """Build a sparse adjacency matrix connecting 26-neighbour voxels.

    Only voxels present in *coords* are nodes.  Two nodes are connected
    iff their Chebyshev distance ≤ 1  (26-connectivity).
    """
    n = len(coords)
    # Map (x,y,z) → index for fast lookup
    coord_to_idx = {}
    for i, (x, y, z) in enumerate(coords):
        coord_to_idx[(int(x), int(y), int(z))] = i

    graph = lil_matrix((n, n), dtype=np.int8)
    for i, (x, y, z) in enumerate(coords):
        x, y, z = int(x), int(y), int(z)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    nb = (x + dx, y + dy, z + dz)
                    j = coord_to_idx.get(nb)
                    if j is not None:
                        graph[i, j] = 1
    return graph.tocsr()


def _merge_tiny_clusters(labels: np.ndarray, features: np.ndarray,
                          min_size: int = MIN_CLUSTER_SIZE) -> np.ndarray:
    """Merge clusters smaller than *min_size* into their nearest large cluster."""
    labels = labels.copy()
    unique, counts = np.unique(labels, return_counts=True)
    large = set(unique[counts >= min_size])
    small = set(unique[counts < min_size])
    if not small:
        return labels

    # Compute centroids of large clusters
    large_centroids = {}
    for cid in large:
        large_centroids[cid] = features[labels == cid].mean(axis=0)

    for cid in small:
        centroid = features[labels == cid].mean(axis=0)
        # Find nearest large cluster
        best_dist = np.inf
        best_target = None
        for lid, lc in large_centroids.items():
            d = np.linalg.norm(centroid - lc)
            if d < best_dist:
                best_dist = d
                best_target = lid
        if best_target is not None:
            labels[labels == cid] = best_target

    # Re-number labels contiguously from 0
    mapping = {old: new for new, old in enumerate(sorted(set(labels)))}
    return np.array([mapping[l] for l in labels])


def silhouette_scan(features: np.ndarray,
                     k_min: int = K_MIN, k_max: int = K_MAX,
                     random_state: int = 42,
                     connectivity=None,
                     use_agglomerative: bool = USE_AGGLOMERATIVE):
    """Try k=k_min..k_max with agglomerative (+ connectivity) or KMeans.

    Returns dict {k: (silhouette_score, calinski_score, labels_array)}.
    """
    results = {}
    for k in range(k_min, k_max + 1):
        if use_agglomerative:
            model = AgglomerativeClustering(
                n_clusters=k,
                connectivity=connectivity,
                linkage='ward',
            )
            labels = model.fit_predict(features)
        else:
            km = KMeans(n_clusters=k, n_init=10, random_state=random_state)
            labels = km.fit_predict(features)

        sil = silhouette_score(features, labels)
        cal = calinski_harabasz_score(features, labels)
        results[k] = (sil, cal, labels)
        print(f'  k={k}  silhouette={sil:.4f}  calinski={cal:.1f}')
    return results


def label_direction(mean_dir: np.ndarray) -> str:
    """Return 'AP', 'DV', or 'ML' based on the axis with the largest |component|."""
    idx = int(np.argmax(np.abs(mean_dir)))
    return DIR_LABELS.get(idx, f'axis{idx}')


def cluster_principal_direction(directions: np.ndarray) -> tuple[np.ndarray, float]:
    """Scatter-matrix eigenvector method for proper axial mean direction.

    Returns (unit_dir, coherence).
    """
    scatter = np.zeros((3, 3))
    for v in directions:
        scatter += np.outer(v, v)
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    mean_dir = eigenvectors[:, -1]
    coherence = eigenvalues[-1] / eigenvalues.sum() if eigenvalues.sum() > 0 else 0.0
    return mean_dir, coherence


# ===========================================================================
#  Main pipeline
# ===========================================================================

def run_fixel_clustering(subject_id: str,
                          session: str,
                          group: str,
                          roi_type: str = 'graft',
                          deriv: Path = DERIV,
                          spatial_weight: float = SPATIAL_WEIGHT,
                          k_min: int = K_MIN,
                          k_max: int = K_MAX,
                          k_force: int | None = None,
                          min_cluster_size: int = MIN_CLUSTER_SIZE,
                          use_agglomerative: bool = USE_AGGLOMERATIVE,
                          random_state: int = 42):
    """Full pipeline: load → cluster → save outputs."""

    subj_dir = deriv / subject_id / session / 'dwi'
    fixel_dir = subj_dir / 'fixels_group'
    out_dir = subj_dir / 'fixel_clusters'
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- 1.  ROI paths (bilateral) ----------------------------------------
    if roi_type == 'graft':
        roi_base = subj_dir.parent / 'manual_roi' / 'fod_space'
        roi_left = roi_base / 'roi_Graft_left_flip_fod.mif'
        roi_right = roi_base / 'roi_Graft_right_flip_fod.mif'
    else:  # isocortex (CTRL)
        roi_base = subj_dir.parent / 'atlas' / 'rois' / 'fod_space'
        roi_left = roi_base / 'roi_isocortex_left.nii.gz'
        roi_right = roi_base / 'roi_isocortex_right.nii.gz'

    print(f'[{subject_id}]  Loading fixels from {fixel_dir.name}/ …')
    fixel_data = load_fixel_dir(fixel_dir)
    spatial_shape = fixel_data['index'].shape[:3]
    print(f'  Fixel grid: {spatial_shape}   total fixels: {len(fixel_data["fd"])}')

    # --- 2.  Bilateral ROI mask -------------------------------------------
    mask_L = load_roi_mask(roi_left)
    mask_R = load_roi_mask(roi_right)
    # Squeeze extra dims (some masks are 4-D with trailing 1s)
    if mask_L.ndim > 3:
        mask_L = mask_L.squeeze()
    if mask_R.ndim > 3:
        mask_R = mask_R.squeeze()
    roi_mask = mask_L | mask_R
    n_roi_vox = int(roi_mask.sum())
    print(f'  Bilateral ROI voxels: {n_roi_vox}  (L={int(mask_L.sum())}, R={int(mask_R.sum())})')

    # --- 3.  Dominant fixel per voxel -------------------------------------
    coords, directions, peaks = extract_dominant_fixels(fixel_data, roi_mask)
    print(f'  Voxels with valid dominant fixel: {len(coords)}')

    if len(coords) < k_max + 1:
        print(f'  WARNING: too few voxels ({len(coords)}) for k_max={k_max}. Aborting.')
        return

    # --- 4.  Feature matrix -----------------------------------------------
    features = build_feature_matrix(coords, directions, spatial_weight)
    print(f'  Feature matrix: {features.shape}')

    # --- 4b. Spatial connectivity for agglomerative clustering -------------
    connectivity = None
    if use_agglomerative:
        print('  Building spatial connectivity graph (26-neighbourhood) …')
        connectivity = _build_connectivity_graph(coords)
        print(f'  Connectivity edges: {connectivity.nnz // 2}')

    # --- 5.  Silhouette scan -----------------------------------------------
    print('  Cluster scan:')
    scan = silhouette_scan(features, k_min, k_max, random_state,
                           connectivity=connectivity,
                           use_agglomerative=use_agglomerative)

    # --- 5a. Select best k -------------------------------------------------
    if k_force is not None:
        # User forced a specific k
        best_k = k_force
        if best_k not in scan:
            print(f'  k_force={best_k} not in scan range; running single fit …')
            if use_agglomerative:
                model = AgglomerativeClustering(
                    n_clusters=best_k, connectivity=connectivity, linkage='ward')
                forced_labels = model.fit_predict(features)
            else:
                km = KMeans(n_clusters=best_k, n_init=10, random_state=random_state)
                forced_labels = km.fit_predict(features)
            best_sil = silhouette_score(features, forced_labels)
            best_cal = calinski_harabasz_score(features, forced_labels)
            scan[best_k] = (best_sil, best_cal, forced_labels)
        best_sil = scan[best_k][0]
        best_cal = scan[best_k][1]
        labels = scan[best_k][2]
        print(f'  → Forced k = {best_k}  (silhouette={best_sil:.4f}  calinski={best_cal:.1f})')
    else:
        # Pick best k using elbow detection on the Calinski-Harabasz curve,
        # clamped to the [K_CLAMP_MIN, K_CLAMP_MAX] window.
        ks = sorted(scan.keys())
        cal_arr = np.array([scan[k][1] for k in ks])
        clamp_lo = K_CLAMP_MIN
        clamp_hi = K_CLAMP_MAX
        window_idx = [i for i, k in enumerate(ks) if clamp_lo <= k <= clamp_hi]

        if len(window_idx) >= 3:
            cal_window = cal_arr[window_idx]
            d2 = np.diff(cal_window, n=2)
            elbow_local = int(np.argmax(np.abs(d2)))
            best_idx = window_idx[elbow_local]
        elif window_idx:
            best_idx = window_idx[int(np.argmax(cal_arr[window_idx]))]
        else:
            mid = (clamp_lo + clamp_hi) // 2
            best_idx = int(np.argmin([abs(k - mid) for k in ks]))

        best_k = ks[best_idx]
        best_sil = scan[best_k][0]
        best_cal = scan[best_k][1]
        labels = scan[best_k][2]
        print(f'  → Best k = {best_k}  (silhouette={best_sil:.4f}  calinski={best_cal:.1f})')

    # --- 5b. Merge tiny clusters -------------------------------------------
    labels = _merge_tiny_clusters(labels, features, min_size=min_cluster_size)
    final_k = len(set(labels))
    if final_k != best_k:
        print(f'  After merging tiny clusters (<{min_cluster_size} vox): {best_k} → {final_k} clusters')
        best_k = final_k

    # --- 6.  Cluster info --------------------------------------------------
    cluster_info = []
    for cid in range(best_k):
        mask_c = labels == cid
        c_dirs = directions[mask_c]
        c_peaks = peaks[mask_c]
        mean_dir, coherence = cluster_principal_direction(c_dirs)
        main_label = label_direction(mean_dir)
        avg_peak = float(c_peaks.mean())

        cluster_info.append({
            'cluster_id': int(cid + 1),         # 1-based for display
            'n_voxels': int(mask_c.sum()),
            'main_direction': main_label,
            'mean_direction_vector': [round(float(v), 4) for v in mean_dir],
            'coherence': round(float(coherence), 4),
            'avg_peak_amplitude': round(avg_peak, 4),
        })
        print(f'    Cluster {cid+1}: {mask_c.sum()} vox  dir={main_label}  '
              f'coh={coherence:.3f}  peak={avg_peak:.4f}  vec={np.round(mean_dir, 3)}')

    # --- 7.  Build label volume -------------------------------------------
    label_vol = np.zeros(spatial_shape, dtype=np.float32)
    for i, (x, y, z) in enumerate(coords):
        label_vol[x, y, z] = labels[i] + 1   # 1-based labels; 0 = background

    # Save as .mif (uses index.mif header for voxel grid)
    ref_mif = fixel_dir / 'index.mif'
    out_mif = out_dir / 'cluster_labels.mif'
    _save_array_as_mif(label_vol, ref_mif, out_mif)
    print(f'  Saved label volume → {out_mif}')

    # Also save NIfTI for convenience
    out_nii = out_dir / 'cluster_labels.nii.gz'
    _save_nii_from_ref(label_vol, ref_mif, out_nii)
    print(f'  Saved NIfTI        → {out_nii}')

    # --- 8.  Save cluster info JSON ----------------------------------------
    info_path = out_dir / 'cluster_info.json'
    meta = {
        'subject': subject_id,
        'session': session,
        'group': group,
        'roi_type': roi_type,
        'n_roi_voxels': n_roi_vox,
        'n_clustered_voxels': len(coords),
        'best_k': best_k,
        'silhouette_score': round(float(best_sil), 4),
        'spatial_weight': spatial_weight,
        'clusters': cluster_info,
        'use_agglomerative': use_agglomerative,
        'min_cluster_size': min_cluster_size,
        'silhouette_all_k': {str(k): round(float(v[0]), 4) for k, v in scan.items()},
        'calinski_all_k': {str(k): round(float(v[1]), 1) for k, v in scan.items()},
    }
    with open(info_path, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f'  Saved cluster info → {info_path}')

    # --- 9.  Silhouette plot -----------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4))
    ks = sorted(scan.keys())
    sil_scores = [scan[k][0] for k in ks]
    cal_scores = [scan[k][1] for k in ks]
    ax.plot(ks, sil_scores, 'o-', color='steelblue', label='silhouette')
    ax2 = ax.twinx()
    ax2.plot(ks, cal_scores, 's--', color='seagreen', label='calinski-harabasz')
    ax.axvline(best_k, ls='--', color='tomato', label=f'best k={best_k}')
    ax.set_xlabel('k')
    ax.set_ylabel('Silhouette score', color='steelblue')
    ax2.set_ylabel('Calinski-Harabasz', color='seagreen')
    ax.set_title(f'{subject_id}  fixel clustering (agglom+connectivity)')
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    fig.tight_layout()
    sil_path = out_dir / 'silhouette_scores.png'
    fig.savefig(sil_path, dpi=150)
    plt.close(fig)
    print(f'  Saved silhouette plot → {sil_path}')

    # --- 10. Per-cluster label volumes (individual masks) -------------------
    for ci in cluster_info:
        cid = ci['cluster_id']
        single = np.zeros(spatial_shape, dtype=np.float32)
        mask_c = labels == (cid - 1)
        for i, (x, y, z) in enumerate(coords):
            if mask_c[i]:
                single[x, y, z] = 1.0
        fname = f'cluster_{cid}_{ci["main_direction"]}.mif'
        _save_array_as_mif(single, ref_mif, out_dir / fname)
        print(f'  Saved individual mask → {fname}')

    print(f'\n  Done — {best_k} clusters in {out_dir}')
    return meta


# ===========================================================================
#  Convenience: save NIfTI with ref header
# ===========================================================================

def _save_nii_from_ref(array: np.ndarray, ref_mif: str | Path,
                        out_nii: str | Path) -> None:
    """Save array as NIfTI using the affine from ref_mif."""
    ref_mif = Path(ref_mif)
    out_nii = Path(out_nii)
    out_nii.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        ref_nii = tmp.name
    try:
        cmd = f'mrconvert "{ref_mif}" "{ref_nii}" -force -quiet'
        if shutil.which('mrconvert'):
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
        else:
            subprocess.run(
                f'micromamba run -n {CONDA_ENV} {cmd}',
                shell=True, check=True, capture_output=True,
            )
        affine = nib.load(ref_nii).affine
    finally:
        if os.path.exists(ref_nii):
            os.remove(ref_nii)

    img = nib.Nifti1Image(array.astype(np.float32), affine)
    nib.save(img, str(out_nii))


# ===========================================================================
#  CLI entry point
# ===========================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fixel cluster analysis')
    parser.add_argument('--subject', '-s', default='sub-AP72-862-901-NmF',
                        help='Subject ID (default: AP72-NmF)')
    parser.add_argument('--all', action='store_true',
                        help='Run on ALL subjects in SUBJECT_MAP')
    parser.add_argument('--group', default=None,
                        help='Run on all subjects in this group (e.g. APTR)')
    parser.add_argument('--spatial-weight', type=float, default=SPATIAL_WEIGHT,
                        help='α weight for spatial vs direction features (default: 2.5)')
    parser.add_argument('--k-min', type=int, default=K_MIN)
    parser.add_argument('--k-max', type=int, default=K_MAX)
    parser.add_argument('--k-force', type=int, default=None,
                        help='Force exactly this many clusters (skip auto selection)')
    parser.add_argument('--min-cluster-size', type=int, default=MIN_CLUSTER_SIZE,
                        help='Merge clusters smaller than this into nearest neighbour')
    parser.add_argument('--no-agglomerative', action='store_true',
                        help='Use KMeans instead of agglomerative clustering')
    args = parser.parse_args()

    # Build subject list
    if args.all:
        subjects = list(SUBJECT_MAP.keys())
    elif args.group:
        subjects = [s for s, v in SUBJECT_MAP.items() if v['group'] == args.group]
        if not subjects:
            raise SystemExit(f'No subjects in group {args.group}')
    else:
        subjects = [args.subject]

    all_results = []
    for sid in subjects:
        if sid not in SUBJECT_MAP:
            print(f'WARNING: {sid} not in SUBJECT_MAP. Skipping.')
            continue
        info = SUBJECT_MAP[sid]
        meta = run_fixel_clustering(
            subject_id=sid,
            session=info['ses'],
            group=info['group'],
            roi_type=info['roi_type'],
            spatial_weight=args.spatial_weight,
            k_min=args.k_min,
            k_max=args.k_max,
            k_force=args.k_force,
            min_cluster_size=args.min_cluster_size,
            use_agglomerative=not args.no_agglomerative,
        )
        if meta:
            all_results.append(meta)

    # Print summary table if multiple subjects
    if len(all_results) > 1:
        print('\n' + '=' * 100)
        print('SUMMARY TABLE  (all subjects)')
        print('=' * 100)
        print(f'{"Subject":<28} {"k":>3} {"Sil":>6} {"Cluster":>8} {"Voxels":>7} '
              f'{"Dir":>4} {"Coh":>6} {"Peak":>6} {"MeanVec":>24}')
        print('-' * 100)
        for meta in all_results:
            sid = meta['subject']
            for ci in meta['clusters']:
                vec_str = ','.join(f'{v:6.3f}' for v in ci['mean_direction_vector'])
                print(f'{sid:<28} {meta["best_k"]:>3} {meta["silhouette_score"]:>6.4f} '
                      f'{ci["cluster_id"]:>8} {ci["n_voxels"]:>7} '
                      f'{ci["main_direction"]:>4} {ci["coherence"]:>6.3f} '
                      f'{ci["avg_peak_amplitude"]:>6.4f} [{vec_str}]')
        print('=' * 100)
