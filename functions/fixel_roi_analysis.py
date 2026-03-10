"""
Fixel ROI Direction Analysis
=============================

Analyze main fiber directions within an arbitrary ROI mask using fixel decomposition.

Pipeline:
  1. Regrid ROI mask to FOD space (if grids differ)
  2. fod2fixel within the ROI → fixel directory
  3. fixel2peaks → peak direction vectors
  4. K-means clustering (k=n_clusters) on |v| for antipodal-safe grouping
  5. Scatter matrix eigenvector method for proper axial averaging
  6. 3D hemisphere scatter + radar plot visualization

Usage:
  from functions.fixel_roi_analysis import run_fixel_roi_analysis

  results = run_fixel_roi_analysis(
      roi_path='/path/to/roi.mif',
      fod_path='/path/to/wm_fod_norm_group.mif',
      output_dir='/path/to/dwi/',
      subject_id='sub-XXX',
      report_dir='/path/to/reports/',
      n_clusters=2,
  )
"""

import subprocess
import numpy as np
import nibabel as nib
from pathlib import Path
from sklearn.cluster import KMeans


# ---------------------------------------------------------------------------
# Core math
# ---------------------------------------------------------------------------

def get_directional_mean(vecs):
    """Compute the true directional mean via the scatter matrix eigenvector method.

    For axial data (v ≡ −v), the arithmetic mean can cancel out.
    The principal eigenvector of S = Σ v_i v_i^T is the correct average.

    Returns
    -------
    mean_dir : np.ndarray, shape (3,)
        Unit vector of the principal direction.
    coherence : float
        Largest eigenvalue / sum of eigenvalues (1.0 = perfectly aligned).
    """
    scatter = np.zeros((3, 3))
    for v in vecs:
        scatter += np.outer(v, v)
    eigenvalues, eigenvectors = np.linalg.eigh(scatter)
    mean_dir = eigenvectors[:, -1]
    coherence = eigenvalues[-1] / eigenvalues.sum() if eigenvalues.sum() > 0 else 0.0
    return mean_dir, coherence


def extract_peak_vectors(peaks_nii_path):
    """Load a fixel2peaks NIfTI and return normalised direction vectors + peak indices.

    Parameters
    ----------
    peaks_nii_path : str or Path
        Path to the peaks NIfTI image (4-D, volumes = 3 * n_peaks_per_voxel).

    Returns
    -------
    all_vectors : np.ndarray, shape (N, 3)
    all_peak_indices : np.ndarray, shape (N,)  — 1-based peak index
    n_peaks_per_voxel : int
    """
    peaks_img = nib.load(str(peaks_nii_path))
    peaks_data = peaks_img.get_fdata()

    n_vols = peaks_data.shape[3] if peaks_data.ndim == 4 else 0
    n_peaks_per_voxel = n_vols // 3 if n_vols >= 3 else 0

    all_vectors = []
    all_peak_indices = []

    for p in range(n_peaks_per_voxel):
        vx = peaks_data[..., p * 3 + 0]
        vy = peaks_data[..., p * 3 + 1]
        vz = peaks_data[..., p * 3 + 2]
        magnitude = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        valid = magnitude > 0.01
        dirs = np.column_stack([vx[valid], vy[valid], vz[valid]])
        norms = np.linalg.norm(dirs, axis=1, keepdims=True)
        dirs = dirs / np.clip(norms, 1e-8, None)
        all_vectors.append(dirs)
        all_peak_indices.extend([p + 1] * len(dirs))

    if all_vectors:
        all_vectors = np.vstack(all_vectors)
    else:
        all_vectors = np.empty((0, 3))
    all_peak_indices = np.array(all_peak_indices)

    return all_vectors, all_peak_indices, n_peaks_per_voxel


def cluster_fixel_directions(all_vectors, all_peak_indices, n_peaks_per_voxel,
                             n_clusters=2, random_state=42, verbose=True):
    """Cluster fixel direction vectors and compute scatter-matrix eigenvector means.

    Parameters
    ----------
    all_vectors : np.ndarray (N, 3)
    all_peak_indices : np.ndarray (N,)
    n_peaks_per_voxel : int
    n_clusters : int
    random_state : int
    verbose : bool

    Returns
    -------
    dict with keys:
        directions : list of np.ndarray  — unit vectors, sorted by cluster size (largest first)
        coherences : list of float
        labels : np.ndarray (N,)
        angle_between : float (degrees, only if n_clusters == 2)
        cluster_sizes : list of int
        all_vectors : np.ndarray (N, 3)  — original vectors (unmodified)
        all_peak_indices : np.ndarray (N,)
        n_peaks_per_voxel : int
        per_cluster : list of dict with per-cluster stats
    """
    axis_names = ['AP (axis 0)', 'DV (axis 1)', 'ML (axis 2)']

    # Cluster on |v| for antipodal safety
    abs_vectors = np.abs(all_vectors)
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
    labels = kmeans.fit_predict(abs_vectors)

    # Scatter matrix eigenvector for each cluster
    per_cluster = []
    for k in range(n_clusters):
        cluster_vecs = all_vectors[labels == k]
        cluster_peaks = all_peak_indices[labels == k]
        mean_dir, coherence = get_directional_mean(cluster_vecs)
        # Ensure consistent sign (flip to z >= 0)
        if mean_dir[2] < 0:
            mean_dir *= -1

        # Angular dispersion
        dots = np.clip(np.abs(np.dot(cluster_vecs, mean_dir)), 0, 1)
        angles = np.degrees(np.arccos(dots))
        mean_angle = float(np.mean(angles))
        dominant_axis = int(np.argmax(np.abs(mean_dir)))

        info = dict(
            direction=mean_dir,
            coherence=coherence,
            n_fixels=len(cluster_vecs),
            mean_angle=mean_angle,
            dominant_axis=dominant_axis,
            dominant_axis_name=axis_names[dominant_axis],
            n_from_peak1=int(np.sum(cluster_peaks == 1)),
            n_from_peak2=int(np.sum(cluster_peaks == 2)) if n_peaks_per_voxel > 1 else 0,
        )

        if verbose:
            d = mean_dir
            print(f'\n  ===== Direction {k+1} ({info["n_fixels"]} fixels) =====')
            print(f'  Mean vector (eigenvector): [{d[0]:+.4f}, {d[1]:+.4f}, {d[2]:+.4f}]')
            print(f'  Dominant axis:  {info["dominant_axis_name"]}')
            print(f'  Coherence:      {coherence:.3f} (1.0 = perfectly aligned)')
            print(f'  Angular spread: {mean_angle:.1f}° (mean deviation from center)')
            print(f'  From peak 1: {info["n_from_peak1"]}, peak 2: {info["n_from_peak2"]}')
            theta = np.degrees(np.arccos(np.clip(np.abs(d[1]), 0, 1)))
            phi = np.degrees(np.arctan2(np.abs(d[2]), np.abs(d[0])))
            print(f'  Spherical:      θ={theta:.1f}° from DV, φ={phi:.1f}° in AP-ML plane')

        per_cluster.append(info)

    # Sort by cluster size (largest first)
    order = sorted(range(n_clusters), key=lambda i: per_cluster[i]['n_fixels'], reverse=True)
    per_cluster = [per_cluster[i] for i in order]
    # Remap labels
    remap = {old: new for new, old in enumerate(order)}
    labels = np.array([remap[lb] for lb in labels])

    directions = [pc['direction'] for pc in per_cluster]
    coherences = [pc['coherence'] for pc in per_cluster]
    cluster_sizes = [pc['n_fixels'] for pc in per_cluster]

    result = dict(
        directions=directions,
        coherences=coherences,
        labels=labels,
        cluster_sizes=cluster_sizes,
        all_vectors=all_vectors,
        all_peak_indices=all_peak_indices,
        n_peaks_per_voxel=n_peaks_per_voxel,
        per_cluster=per_cluster,
    )

    # Angle between first two directions
    if n_clusters >= 2:
        d1, d2 = directions[0], directions[1]
        result['angle_between'] = float(
            np.degrees(np.arccos(np.clip(np.abs(np.dot(d1, d2)), 0, 1)))
        )
        if verbose:
            print(f'\n  ===== Summary =====')
            print(f'  Direction 1: [{d1[0]:+.4f}, {d1[1]:+.4f}, {d1[2]:+.4f}]')
            print(f'  Direction 2: [{d2[0]:+.4f}, {d2[1]:+.4f}, {d2[2]:+.4f}]')
            print(f'  Angle between directions: {result["angle_between"]:.1f}°')
            print(f'  Total fixels analyzed: {len(all_vectors)}')
            print(f'  Method: scatter matrix eigenvector (proper axial averaging)')

    return result


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def draw_radial_whiskers(ax3d, vectors, peaks, color='#E24A33', alpha=0.3,
                         linewidth=0.6, max_length=0.4):
    """Draw radial whiskers from unit-sphere points, length ∝ peak amplitude.

    Parameters
    ----------
    ax3d : mpl_toolkits.mplot3d.Axes3D
    vectors : np.ndarray (N, 3)  — unit directions (original axis order: 0=AP, 1=DV, 2=ML)
    peaks : np.ndarray (N,)  — scalar amplitudes
    color, alpha, linewidth : styling
    max_length : float  — length of whisker for the maximum peak value
    """
    pk_max = peaks.max() if len(peaks) > 0 and peaks.max() > 0 else 1.0
    scale = max_length / pk_max
    for i in range(len(vectors)):
        v = vectors[i]
        length = peaks[i] * scale
        # Swapped axes: X=AP(0), Y=ML(2), Z=DV(1)
        ax3d.plot([v[0], v[0] * (1 + length)],
                  [v[2], v[2] * (1 + length)],
                  [v[1], v[1] * (1 + length)],
                  color=color, alpha=alpha, linewidth=linewidth)


def plot_fixel_directions(result, subject_id='', save_path=None, show=True,
                          peaks=None):
    """Generate 3D hemisphere scatter + radar plot for fixel direction results.

    Parameters
    ----------
    result : dict
        Output from cluster_fixel_directions().
    subject_id : str
        For the plot title.
    save_path : str or Path or None
        If provided, save the figure here.
    show : bool
        Whether to call plt.show().
    peaks : np.ndarray or None
        Per-fixel peak FOD amplitudes (same length as all_vectors).
        If provided, radial whiskers are drawn from each point on the
        sphere with length proportional to the peak value.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    directions = result['directions']
    labels = result['labels']
    all_vectors = result['all_vectors'].copy()
    n_clusters = len(directions)
    angle_between = result.get('angle_between', None)

    # Flip vectors to DV >= 0 for viz (DV = axis 1)
    flip = all_vectors[:, 1] < 0
    all_vectors[flip] *= -1

    # Prepare peaks array (flip in sync with vectors)
    if peaks is not None:
        peaks = np.asarray(peaks, dtype=float).copy()
        # peaks don't need sign flip, they are scalar magnitudes

    colors_list = ['#E24A33', '#348ABD', '#988ED5', '#8EBA42', '#FFB347', '#77BEDB']
    colors = colors_list[:n_clusters]

    fig = plt.figure(figsize=(16, 7))
    fig.suptitle(f'Fixel ROI — Main Fiber Directions ({subject_id})',
                 fontsize=14, fontweight='bold')

    # ---- Left: 3D binned sticks (origin → bin direction, length = sum of peaks) ----
    ax3d = fig.add_subplot(121, projection='3d')

    from matplotlib import cm
    pk_arr = np.asarray(peaks, dtype=float) if peaks is not None else np.ones(len(all_vectors))
    cmap = cm.get_cmap('hot_r')

    # Bin fixels on a spherical grid and sum peaks per bin
    n_az, n_el = 16, 8
    az_edges = np.linspace(0, 2 * np.pi, n_az + 1)
    el_edges = np.linspace(0, np.pi / 2, n_el + 1)

    # Spherical coords: elevation = arcsin(DV), azimuth = atan2(ML, AP)
    elev_vals = np.arcsin(np.clip(all_vectors[:, 1], -1, 1))
    az_vals = np.arctan2(all_vectors[:, 2], all_vectors[:, 0]) % (2 * np.pi)

    peak_grid = np.zeros((n_az, n_el))
    for fi in range(len(all_vectors)):
        ui = np.clip(np.searchsorted(az_edges, az_vals[fi], side='right') - 1, 0, n_az - 1)
        vi = np.clip(np.searchsorted(el_edges, elev_vals[fi], side='right') - 1, 0, n_el - 1)
        peak_grid[ui, vi] += pk_arr[fi]

    pk_max = peak_grid.max() if peak_grid.max() > 0 else 1.0

    # Draw one stick per non-empty bin: direction = bin centre, length = sum(peaks)
    for i in range(n_az):
        for j in range(n_el):
            if peak_grid[i, j] <= 0:
                continue
            # Bin centre in spherical coords
            az_c = (az_edges[i] + az_edges[i + 1]) / 2
            el_c = (el_edges[j] + el_edges[j + 1]) / 2
            # Convert to cartesian (original axes: AP, DV, ML)
            dx = np.cos(az_c) * np.cos(el_c)  # AP
            dy = np.sin(el_c)                   # DV
            dz = np.sin(az_c) * np.cos(el_c)   # ML
            length = peak_grid[i, j] / pk_max   # normalised to [0, 1]
            # Plot axes: X=AP, Y=ML, Z=DV
            ax3d.plot([0, dx * length],
                      [0, dz * length],
                      [0, dy * length],
                      color=cmap(length), alpha=0.8, linewidth=1.5)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=0, vmax=pk_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax3d, shrink=0.5, pad=0.08, label='Sum peak amplitude')

    ax3d.set_xlabel('AP (axis 0)')
    ax3d.set_ylabel('ML (axis 2)')
    ax3d.set_zlabel('DV (axis 1)')
    ax3d.set_xlim(-1.05, 1.05)
    ax3d.set_ylim(-1.05, 1.05)
    ax3d.set_zlim(-0.05, 1.05)
    ax3d.set_title('Binned Fixel Sticks (length = sum peak amp)')
    ax3d.view_init(elev=25, azim=-60)

    # Remove background panes and grid, keep only axes
    ax3d.xaxis.pane.fill = False
    ax3d.yaxis.pane.fill = False
    ax3d.zaxis.pane.fill = False
    ax3d.xaxis.pane.set_edgecolor('none')
    ax3d.yaxis.pane.set_edgecolor('none')
    ax3d.zaxis.pane.set_edgecolor('none')
    ax3d.grid(False)

    # ---- Right: Radar plot (DV south, ML upper-right, AP upper-left) ----
    # Triangles are scaled by mean peak amplitude per cluster, not unitless.
    ax_radar = fig.add_subplot(122, projection='polar')
    ax_radar.set_theta_zero_location('S')
    ax_radar.set_theta_direction(-1)

    radar_labels_list = ['DV', 'ML', 'AP']
    angles_radar = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML → [DV, ML, AP]

    linestyles = ['-', '--', '-.', ':']

    # Compute mean peak amplitude per cluster for scaling
    cluster_mean_peaks = []
    for k_idx in range(n_clusters):
        mask_k = labels == k_idx
        if peaks is not None and np.any(mask_k):
            cluster_mean_peaks.append(float(np.mean(peaks[mask_k])))
        else:
            cluster_mean_peaks.append(1.0)  # fallback: unit vectors

    radar_max = 0.0
    for k_idx, d in enumerate(directions):
        # Scale direction components by mean peak amplitude of this cluster
        abs_d = np.abs(d)[reorder] * cluster_mean_peaks[k_idx]
        radar_max = max(radar_max, np.max(abs_d))
        vals = np.concatenate([abs_d, [abs_d[0]]])
        angles_closed = np.concatenate([angles_radar, [angles_radar[0]]])
        ls = linestyles[k_idx % len(linestyles)]
        n_k = int(np.sum(labels == k_idx))
        mean_pk = cluster_mean_peaks[k_idx]
        ax_radar.plot(angles_closed, vals, color=colors[k_idx], linewidth=2.5,
                      linestyle=ls,
                      label=f'Dir {k_idx + 1} (n={n_k}, peak={mean_pk:.3f})')
        ax_radar.fill(angles_closed, vals, color=colors[k_idx],
                      alpha=0.15 if k_idx == 0 else 0.08)

    # Isotropic reference scaled by overall mean peak
    if peaks is not None and len(peaks) > 0:
        overall_mean_peak = float(np.mean(peaks))
    else:
        overall_mean_peak = 1.0
    iso_val = overall_mean_peak / np.sqrt(3)
    iso_angles = np.linspace(0, 2 * np.pi, 100)
    ax_radar.plot(iso_angles, [iso_val] * len(iso_angles), color='gray',
                  linewidth=0.8, linestyle=':', alpha=0.5, label='isotropic')

    ax_radar.set_xticks(angles_radar)
    ax_radar.set_xticklabels(radar_labels_list, fontsize=13, fontweight='bold')
    # Dynamic y-axis based on actual peak-scaled values
    y_max = max(radar_max * 1.1, iso_val * 1.2)
    ax_radar.set_ylim(0, y_max)
    # Auto-compute nice tick spacing
    tick_step = round(y_max / 4, 3)
    if tick_step > 0:
        yticks = [tick_step * i for i in range(1, 5) if tick_step * i < y_max]
        ax_radar.set_yticks(yticks)
        ax_radar.set_yticklabels([f'{t:.3f}' for t in yticks], fontsize=8, color='gray')
    angle_str = f'(∠ = {angle_between:.1f}°)' if angle_between is not None else ''
    ax_radar.set_title(f'Peak-weighted Direction\n{angle_str}', pad=20, fontsize=12)
    ax_radar.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), fontsize=9, framealpha=0.9)

    plt.tight_layout()

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')

    if show:
        plt.show()

    return fig


# ---------------------------------------------------------------------------
# MRtrix shell helpers
# ---------------------------------------------------------------------------

def _run_cmd(cmd, timeout=300):
    """Run a shell command, return (returncode, stdout, stderr)."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        print(f'  ⚠ Command failed (exit {result.returncode}):')
        if result.stderr:
            print(f'  STDERR: {result.stderr.strip()}')
    return result.returncode, result.stdout, result.stderr


def generate_fixel_peaks(roi_path, fod_path, output_dir, conda_env='mbrain', force=True):
    """Run the MRtrix pipeline: regrid ROI → fod2fixel → fixel2peaks → NIfTI export.

    Parameters
    ----------
    roi_path : str or Path
        Path to the ROI mask (any grid — will be regridded to FOD space if needed).
    fod_path : str or Path
        Path to the FOD image (e.g. wm_fod_norm_group.mif).
    output_dir : str or Path
        Directory for intermediate files (fixels, peaks, etc.).
    conda_env : str
        Name of the conda/micromamba environment with MRtrix.
    force : bool
        Whether to overwrite existing outputs.

    Returns
    -------
    peaks_nii : Path
        Path to the exported peaks NIfTI.
    """
    roi_path = Path(roi_path)
    fod_path = Path(fod_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = roi_path.stem.replace('.mif', '').replace('.nii', '')
    roi_fod = output_dir / f'{prefix}_fod.mif'
    fixel_dir = output_dir / f'fixels_{prefix}'
    peaks_mif = output_dir / f'{prefix}_peaks.mif'
    peaks_nii = output_dir / f'{prefix}_peaks.nii'
    force_flag = '-force' if force else ''
    env = f'micromamba run -n {conda_env}'

    # Phase 0: Regrid ROI to FOD space
    print(f'\n--- Regrid ROI to FOD space ---')
    rc, _, _ = _run_cmd(
        f'{env} mrgrid {roi_path} regrid {roi_fod} '
        f'-template {fod_path} -interp nearest -datatype bit {force_flag} -quiet',
        timeout=120)
    if rc != 0:
        raise RuntimeError(f'mrgrid regrid failed for {roi_path}')

    # Verify non-zero voxels
    _, stdout, _ = _run_cmd(f'{env} mrstats {roi_fod} -output mean -output count', timeout=30)
    print(f'  ROI (FOD space): {stdout.strip()}')

    # Phase 1: fod2fixel
    print(f'\n--- fod2fixel within ROI ---')
    rc, _, _ = _run_cmd(
        f'{env} fod2fixel {fod_path} {fixel_dir} '
        f'-mask {roi_fod} -afd afd.mif -peak peak_amp.mif '
        f'-disp dispersion.mif {force_flag} -quiet',
        timeout=300)
    if rc != 0:
        raise RuntimeError(f'fod2fixel failed')
    for f in ['index.mif', 'directions.mif', 'afd.mif', 'peak_amp.mif']:
        p = fixel_dir / f
        status = '✅' if p.exists() else '❌ MISSING'
        print(f'  {status}: {f}')

    # Phase 2: fixel2peaks
    print(f'\n--- fixel2peaks ---')
    rc, _, _ = _run_cmd(
        f'{env} fixel2peaks {fixel_dir}/directions.mif {peaks_mif} {force_flag} -quiet',
        timeout=120)
    if rc != 0:
        raise RuntimeError(f'fixel2peaks failed')

    # Export to NIfTI
    _run_cmd(f'{env} mrconvert {peaks_mif} {peaks_nii} {force_flag} -quiet', timeout=60)
    print(f'  ✅ Peaks exported: {peaks_nii}')

    return peaks_nii


# ---------------------------------------------------------------------------
# High-level API
# ---------------------------------------------------------------------------

def run_fixel_roi_analysis(roi_path, fod_path, output_dir, subject_id='',
                           report_dir=None, n_clusters=2, conda_env='mbrain',
                           force=True, show=True, verbose=True):
    """Full pipeline: regrid → fixels → cluster → plot.

    Parameters
    ----------
    roi_path : str or Path
        ROI mask (any grid).
    fod_path : str or Path
        FOD image.
    output_dir : str or Path
        Working directory for fixel intermediates.
    subject_id : str
        For titles and filenames.
    report_dir : str or Path or None
        Where to save the figure. If None, uses output_dir.
    n_clusters : int
        Number of fiber populations to identify (default 2).
    conda_env : str
        Micromamba environment name.
    force : bool
        Whether to overwrite existing files.
    show : bool
        Whether to display the plot inline.
    verbose : bool
        Print progress and statistics.

    Returns
    -------
    dict with keys: directions, coherences, labels, angle_between, figure, ...
    """
    roi_path = Path(roi_path)
    output_dir = Path(output_dir)

    if not roi_path.exists():
        raise FileNotFoundError(f'ROI not found: {roi_path}')
    if not Path(fod_path).exists():
        raise FileNotFoundError(f'FOD not found: {fod_path}')

    # Step 1: Generate fixel peaks
    peaks_nii = generate_fixel_peaks(roi_path, fod_path, output_dir,
                                     conda_env=conda_env, force=force)

    # Step 2: Extract and cluster
    all_vectors, all_peak_indices, n_peaks = extract_peak_vectors(peaks_nii)
    if len(all_vectors) == 0:
        print('  ⚠ No valid fixel vectors found in the ROI.')
        return None

    if verbose:
        print(f'\n  Total fixel vectors: {len(all_vectors)}')
        print(f'    Peak 1: {np.sum(all_peak_indices == 1)}')
        if n_peaks > 1:
            print(f'    Peak 2: {np.sum(all_peak_indices == 2)}')

    result = cluster_fixel_directions(all_vectors, all_peak_indices, n_peaks,
                                      n_clusters=n_clusters, verbose=verbose)

    # Step 3: Plot
    if report_dir is None:
        report_dir = output_dir
    report_dir = Path(report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    roi_stem = roi_path.stem.replace('.mif', '').replace('.nii', '')
    fig_name = f'fixel_roi_{roi_stem}_{subject_id}.png' if subject_id else f'fixel_roi_{roi_stem}.png'
    save_path = report_dir / fig_name

    fig = plot_fixel_directions(result, subject_id=subject_id,
                                save_path=save_path, show=show)
    result['figure'] = fig
    result['figure_path'] = save_path

    return result


# ---------------------------------------------------------------------------
# Fixel loading helpers (for pre-existing fixel directories)
# ---------------------------------------------------------------------------

def _load_mif_as_array(mif_path, conda_env='mbrain'):
    """Convert .mif to temporary NIfTI and load as numpy array."""
    import tempfile, os, shutil
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        # Try bare mrconvert first (works when mrtrix is loaded via modules)
        if shutil.which('mrconvert'):
            subprocess.run(
                f'mrconvert "{mif_path}" "{tmp_path}" -force -quiet',
                shell=True, check=True, capture_output=True)
        else:
            # Fall back to micromamba wrapper
            subprocess.run(
                f'micromamba run -n {conda_env} mrconvert "{mif_path}" "{tmp_path}" -force -quiet',
                shell=True, check=True, capture_output=True)
        data = nib.load(tmp_path).get_fdata()
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return data


def load_fixel_dir(fixel_dir, conda_env='mbrain'):
    """Load fixel data from an fod2fixel output directory.

    Returns dict with keys: index (X,Y,Z,2), directions (N,3), fd (N,), peak_amp (N,).
    """
    fixel_dir = Path(fixel_dir)
    index = _load_mif_as_array(fixel_dir / 'index.mif', conda_env)
    directions_raw = _load_mif_as_array(fixel_dir / 'directions.mif', conda_env)
    afd_path = fixel_dir / 'afd.mif'
    fd_path = fixel_dir / 'fd.mif'
    density_path = afd_path if afd_path.exists() else fd_path
    fd = _load_mif_as_array(density_path, conda_env).flatten()

    # Load peak FOD amplitude (more sensitive to fiber presence vs noise)
    peak_path = fixel_dir / 'peak.mif'
    if not peak_path.exists():
        peak_path = fixel_dir / 'peak_amp.mif'
    if peak_path.exists():
        peak_amp = _load_mif_as_array(peak_path, conda_env).flatten()
    else:
        peak_amp = fd.copy()  # fallback to AFD

    dirs = directions_raw[:, :, 0] if directions_raw.ndim == 3 else directions_raw
    return dict(index=index, directions=dirs, fd=fd, peak_amp=peak_amp)


def load_roi_mask(roi_path, conda_env='mbrain'):
    """Load an ROI mask, handling both .mif and .nii/.nii.gz."""
    roi_path = Path(roi_path)
    if roi_path.suffix == '.mif':
        return _load_mif_as_array(roi_path, conda_env) > 0.5
    else:
        return nib.load(str(roi_path)).get_fdata() > 0.5


def extract_slicewise_vectors(fixel_data, roi_mask, slice_axis=0, min_fixels=5):
    """Extract fixel direction vectors, AFD, and peak amplitude per slice within an ROI.

    Parameters
    ----------
    fixel_data : dict  (from load_fixel_dir)
    roi_mask : np.ndarray  (bool, same spatial shape as fixel index)
    slice_axis : int  (0=AP, 1=DV, 2=ML in flip1 convention)
    min_fixels : int  (skip slices with fewer fixels)

    Returns
    -------
    slice_vectors : dict mapping slice_index -> np.ndarray of shape (N, 3)
    slice_peaks : dict mapping slice_index -> np.ndarray of shape (N,)  [peak FOD amplitude]
    slice_roi_voxels : dict mapping slice_index -> int  (total ROI voxels in that slice)
    """
    index = fixel_data['index']
    dirs = fixel_data['directions']
    fd = fixel_data['fd']
    peak_amp = fixel_data.get('peak_amp', fd)  # fallback to fd if no peak

    n_slices = roi_mask.shape[slice_axis]
    slice_vectors = {}
    slice_peaks = {}
    slice_roi_voxels = {}

    for s in range(n_slices):
        slicer = [slice(None)] * 3
        slicer[slice_axis] = s
        roi_slice = roi_mask[tuple(slicer)]
        n_roi_vox = int(roi_slice.sum())
        if n_roi_vox == 0:
            continue

        mask_3d = np.zeros_like(roi_mask, dtype=bool)
        mask_3d[tuple(slicer)] = roi_slice

        coords = np.where(mask_3d)
        vecs = []
        peaks = []
        for x, y, z in zip(coords[0], coords[1], coords[2]):
            count = int(index[x, y, z, 0])
            offset = int(index[x, y, z, 1])
            if count == 0:
                continue
            for i in range(count):
                fi = offset + i
                if fi >= len(fd):
                    break
                if fd[fi] < 0.001:
                    continue
                d = dirs[fi]
                norm = np.linalg.norm(d)
                if norm > 0.01:
                    vecs.append(d / norm)
                    peaks.append(peak_amp[fi] if fi < len(peak_amp) else fd[fi])

        if len(vecs) >= min_fixels:
            slice_vectors[s] = np.array(vecs)
            slice_peaks[s] = np.array(peaks)
            slice_roi_voxels[s] = n_roi_vox

    return slice_vectors, slice_peaks, slice_roi_voxels


def slicewise_kmeans(slice_vectors, n_clusters=2, random_state=42,
                     slice_afds=None, slice_roi_voxels=None):
    """Run K-means + peak-weighted per-component averaging on each slice.

    For each cluster, compute the reference direction via scatter-matrix
    eigenvector, flip all fixel vectors to align with it, then compute
    ``sum(peak_i * aligned_direction_i) / n_roi_voxels`` per component.
    Dividing by **ROI voxels** (not fixels) means empty voxels contribute
    zero, correctly reflecting regions with sparse/no real fibers.

    Parameters
    ----------
    slice_vectors : dict  (slice_idx -> (N, 3) array of unit directions)
    n_clusters : int
    random_state : int
    slice_afds : dict or None  (slice_idx -> (N,) array of peak amplitude values)
    slice_roi_voxels : dict or None  (slice_idx -> int, total ROI voxels per slice)
        If provided, scaled_directions are sum(peak*dir)/n_voxels.
        If None, falls back to mean over fixels (old behavior).

    Returns
    -------
    dict mapping slice_idx -> dict with keys:
        directions : list of np.ndarray  (unit eigenvectors, sorted by cluster size)
        scaled_directions : list of np.ndarray  (per-component peak-weighted, voxel-normalized)
        mean_peaks : list of float  (mean peak amplitude per cluster)
        coherences : list of float
        angle_between : float (degrees, n_clusters==2 only)
        n_fixels : int
        n_roi_voxels : int  (total ROI voxels in this slice)
        cluster_sizes : list of int
    """
    results = {}
    for s_idx, vecs in slice_vectors.items():
        if len(vecs) < n_clusters:
            continue
        abs_v = np.abs(vecs)
        km = KMeans(n_clusters=n_clusters, n_init=10, random_state=random_state)
        labels = km.fit_predict(abs_v)

        # Get per-fixel peak amplitude for this slice (if available)
        peaks = slice_afds[s_idx] if (slice_afds is not None and s_idx in slice_afds) else None

        cluster_info = []
        for k in range(n_clusters):
            mask_k = labels == k
            cv = vecs[mask_k]
            if len(cv) == 0:
                cluster_info.append(dict(direction=np.array([0, 0, 0.]),
                                         scaled_dir=np.array([0, 0, 0.]),
                                         coherence=0.0, n=0, mean_peak=0.0))
                continue

            # 1) Reference direction via scatter-matrix eigenvector
            mean_dir, coh = get_directional_mean(cv)
            if mean_dir[2] < 0:
                mean_dir *= -1

            # 2) Align all fixel vectors to reference direction
            signs = np.sign(cv @ mean_dir)
            signs[signs == 0] = 1.0
            aligned = cv * signs[:, None]  # (N, 3)

            # 3) Per-component peak-weighted sum, divided by ROI voxels
            n_vox = (slice_roi_voxels[s_idx]
                     if slice_roi_voxels is not None and s_idx in slice_roi_voxels
                     else len(cv))  # fallback: divide by fixel count
            if peaks is not None:
                pk = peaks[mask_k]
                weighted_sum = np.sum(aligned * pk[:, None], axis=0)  # (3,)
                scaled_dir = weighted_sum / n_vox
                mean_peak = float(np.sum(pk)) / n_vox
            else:
                scaled_dir = np.sum(aligned, axis=0) / n_vox
                mean_peak = 1.0

            cluster_info.append(dict(direction=mean_dir, scaled_dir=scaled_dir,
                                     coherence=coh, n=len(cv), mean_peak=mean_peak))

        # Sort by cluster size (largest first)
        order = sorted(range(n_clusters), key=lambda i: cluster_info[i]['n'], reverse=True)
        cluster_info = [cluster_info[i] for i in order]

        n_vox = (slice_roi_voxels[s_idx]
                 if slice_roi_voxels is not None and s_idx in slice_roi_voxels
                 else len(vecs))
        res = dict(
            directions=[ci['direction'] for ci in cluster_info],
            scaled_directions=[ci['scaled_dir'] for ci in cluster_info],
            mean_peaks=[ci['mean_peak'] for ci in cluster_info],
            coherences=[ci['coherence'] for ci in cluster_info],
            n_fixels=len(vecs),
            n_roi_voxels=n_vox,
            cluster_sizes=[ci['n'] for ci in cluster_info],
        )
        if n_clusters >= 2:
            d1, d2 = res['directions'][0], res['directions'][1]
            res['angle_between'] = float(
                np.degrees(np.arccos(np.clip(np.abs(np.dot(d1, d2)), 0, 1))))
        results[s_idx] = res

    return results


# ---------------------------------------------------------------------------
# 3-group comparison figure (K-means per slice)
# ---------------------------------------------------------------------------

def plot_3group_kmeans_comparison(all_group_data, hemisphere, target_n_slices,
                                  use_shrink=False, save_path=None, show=False,
                                  y_max=None):
    """Generate the 3-group fixel K-means comparison figure.

    Radar amplitude is per-component peak-weighted, so regions with no
    real fibers (e.g. AP graft) collapse to near-zero.

    Parameters
    ----------
    all_group_data : dict
        {group_name: {subject_id: {slice_idx: kmeans_result_dict, ...}, ...}, ...}
    hemisphere : str  ('left' or 'right')
    target_n_slices : int
    use_shrink : bool
    save_path : str or Path or None
    show : bool
    y_max : float or None
        Fixed y-axis maximum for consistent cross-plot comparison.
        If None, computed from data.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    row_groups = ['CTRL', 'AP', 'APTR']
    group_labels = ['CTRL\n(isocortex)', 'AP\n(graft, no Tx)', 'APTR\n(graft + Tx)']
    group_bg = ['#e8f5e9', '#ffebee', '#e3f2fd']
    subject_colors = ['#e74c3c', '#3498db', '#27ae60', '#9b59b6', '#f39c12']

    # Radar setup: DV at south, ML upper-right, AP upper-left
    radar_labels_list = ['DV', 'ML', 'AP']
    angles_radar = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    reorder = [1, 2, 0]  # d[0]=AP, d[1]=DV, d[2]=ML -> [DV, ML, AP]
    angles_closed = np.concatenate([angles_radar, [angles_radar[0]]])

    n_cols = target_n_slices + 3   # slices + AVG/subj + GROUP AVG + INFO
    n_rows = len(row_groups)

    # ---- Compute y-axis max (use forced value or derive from data) ----
    if y_max is None:
        global_max = 0.0
        for group in row_groups:
            for subj, data in all_group_data.get(group, {}).items():
                for s_idx, res in data.items():
                    for sd in res.get('scaled_directions', res.get('directions', [])):
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val
        # Round up with some headroom; ensure minimum
        y_max = max(0.05, np.ceil(global_max * 20) / 20)  # round up to nearest 0.05

    fig = plt.figure(figsize=(2.4 * n_cols, 3.5 * n_rows))

    def _draw_radar(ax, directions_list, color, linestyle='-', linewidth=1.5, alpha_fill=0.08):
        """Draw a single 2-direction radar on a polar axes."""
        for k_idx, d in enumerate(directions_list):
            abs_d = np.abs(d)[reorder]
            vals = np.concatenate([abs_d, [abs_d[0]]])
            ls = '-' if k_idx == 0 else '--'
            ax.plot(angles_closed, vals, color=color, linewidth=linewidth,
                    linestyle=ls, alpha=0.85)
            ax.fill(angles_closed, vals, color=color,
                    alpha=alpha_fill if k_idx == 0 else alpha_fill * 0.5)

    def _setup_radar_axes(ax, bg, show_labels=False):
        """Common radar axes setup."""
        ax.set_facecolor(bg)
        ax.set_theta_zero_location('S')
        ax.set_theta_direction(-1)
        ax.set_ylim(0, y_max * 1.05)
        # Nice tick spacing
        tick_step = max(0.01, round(y_max / 4, 2))
        ticks = np.arange(tick_step, y_max + 0.001, tick_step)
        ax.set_yticks(ticks)
        if show_labels:
            ax.set_xticks(angles_radar)
            ax.set_xticklabels(radar_labels_list, fontsize=7, fontweight='bold')
            tick_labels = [f'{t:.2f}' if t < 0.1 else f'{t:.1f}' for t in ticks]
            # Show only half, first empty
            labels = ['' if i % 2 == 0 else tick_labels[i] for i in range(len(ticks))]
            ax.set_yticklabels(labels, fontsize=5, color='gray')
        else:
            ax.set_xticks([])
            ax.set_yticklabels([])
        ax.grid(True, alpha=0.25)

    for row_idx, (group, label, bg) in enumerate(zip(row_groups, group_labels, group_bg)):
        subjects_data = all_group_data.get(group, {})
        subject_ids = list(subjects_data.keys())

        # ---- Per-slice columns ----
        for col_idx in range(target_n_slices):
            ax_idx = row_idx * n_cols + col_idx + 1
            ax = fig.add_subplot(n_rows, n_cols, ax_idx, projection='polar')
            _setup_radar_axes(ax, bg, show_labels=(col_idx == 0))

            for si, subj in enumerate(subject_ids):
                color = subject_colors[si % len(subject_colors)]
                data = subjects_data[subj]
                if col_idx in data:
                    dirs_key = 'scaled_directions' if 'scaled_directions' in data[col_idx] else 'directions'
                    _draw_radar(ax, data[col_idx][dirs_key], color)

            if row_idx == 0:
                ax.set_title(f'Slice {col_idx}', fontsize=9, fontweight='bold', pad=8)

        # ---- AVG/subj column: average scaled_directions per-component across slices ----
        ax_avg = fig.add_subplot(n_rows, n_cols,
                                 row_idx * n_cols + target_n_slices + 1, projection='polar')
        _setup_radar_axes(ax_avg, '#fffde7')

        for si, subj in enumerate(subject_ids):
            color = subject_colors[si % len(subject_colors)]
            data = subjects_data[subj]
            # Collect scaled_directions (already per-component peak-weighted)
            sd1_list = [data[s]['scaled_directions'][0] for s in range(target_n_slices)
                        if s in data and 'scaled_directions' in data[s]]
            sd2_list = [data[s]['scaled_directions'][1] for s in range(target_n_slices)
                        if s in data and 'scaled_directions' in data[s]
                        and len(data[s]['scaled_directions']) > 1]

            avg_dirs = []
            if sd1_list:
                avg_dirs.append(np.mean(sd1_list, axis=0))
            if sd2_list:
                avg_dirs.append(np.mean(sd2_list, axis=0))
            if avg_dirs:
                _draw_radar(ax_avg, avg_dirs, color)

        if row_idx == 0:
            ax_avg.set_title('AVG/subj', fontsize=8, fontweight='bold', pad=8, color='darkred')

        # ---- GROUP AVG column: average scaled_directions across all subjects × slices ----
        ax_grp = fig.add_subplot(n_rows, n_cols,
                                 row_idx * n_cols + target_n_slices + 2, projection='polar')
        _setup_radar_axes(ax_grp, '#ffecb3')

        all_sd1, all_sd2 = [], []
        for subj in subject_ids:
            data = subjects_data[subj]
            for s in range(target_n_slices):
                if s not in data:
                    continue
                if 'scaled_directions' in data[s]:
                    all_sd1.append(data[s]['scaled_directions'][0])
                    if len(data[s]['scaled_directions']) > 1:
                        all_sd2.append(data[s]['scaled_directions'][1])

        grp_dirs = []
        grp_dirs_unit = []  # for INFO panel
        if all_sd1:
            grp_sd1 = np.mean(all_sd1, axis=0)
            grp_dirs.append(grp_sd1)
            norm1 = np.linalg.norm(grp_sd1)
            grp_dirs_unit.append(grp_sd1 / norm1 if norm1 > 1e-8 else grp_sd1)
        if all_sd2:
            grp_sd2 = np.mean(all_sd2, axis=0)
            grp_dirs.append(grp_sd2)
            norm2 = np.linalg.norm(grp_sd2)
            grp_dirs_unit.append(grp_sd2 / norm2 if norm2 > 1e-8 else grp_sd2)
        if grp_dirs:
            _draw_radar(ax_grp, grp_dirs, 'black', linewidth=3, alpha_fill=0.15)

        if row_idx == 0:
            ax_grp.set_title('GROUP\nAVG', fontsize=8, fontweight='bold', pad=8, color='darkred')

        # ---- INFO panel ----
        ax_info = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + n_cols)
        ax_info.set_facecolor(bg)
        ax_info.axis('off')

        n_subj = len(subject_ids)
        all_angles = []
        for subj in subject_ids:
            for s in range(target_n_slices):
                if s in subjects_data[subj] and 'angle_between' in subjects_data[subj][s]:
                    all_angles.append(subjects_data[subj][s]['angle_between'])
        mean_angle = np.mean(all_angles) if all_angles else 0

        # Compute group mean peak amplitude from scaled directions magnitudes
        all_norms = [np.linalg.norm(sd) for sd in all_sd1 + all_sd2]
        grp_mean_peak = float(np.mean(all_norms)) if all_norms else 0.0

        info_text = f'{label}\n\nn={n_subj} subjects\nmean Peak={grp_mean_peak:.3f}\n'
        if grp_dirs_unit:
            g1 = grp_dirs_unit[0]
            info_text += f'\nDir1: [{g1[0]:+.2f},{g1[1]:+.2f},{g1[2]:+.2f}]'
            dominant = ['AP', 'DV', 'ML'][int(np.argmax(np.abs(g1)))]
            info_text += f'\nDominant: {dominant}'
        if len(grp_dirs_unit) > 1:
            g2 = grp_dirs_unit[1]
            info_text += f'\n\nDir2: [{g2[0]:+.2f},{g2[1]:+.2f},{g2[2]:+.2f}]'
            dominant2 = ['AP', 'DV', 'ML'][int(np.argmax(np.abs(g2)))]
            info_text += f'\nDominant: {dominant2}'
        if all_angles:
            info_text += f'\n\n∠ Dir1-Dir2:\n{mean_angle:.1f}° (mean)'

        ax_info.text(0.5, 0.5, info_text, transform=ax_info.transAxes,
                     fontsize=8, ha='center', va='center',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Legend
    max_subj = max(len(all_group_data.get(g, {})) for g in row_groups) if all_group_data else 1
    legend_elements = [plt.Line2D([0], [0], color=subject_colors[i], linewidth=2,
                                   marker='o', label=f'Mouse {i+1}')
                       for i in range(min(5, max_subj))]
    legend_elements.append(plt.Line2D([0], [0], color='black', linewidth=1.5,
                                       linestyle='-', label='Dir 1 (solid)'))
    legend_elements.append(plt.Line2D([0], [0], color='black', linewidth=1.5,
                                       linestyle='--', label='Dir 2 (dashed)'))
    fig.legend(handles=legend_elements, loc='upper right', fontsize=8,
               bbox_to_anchor=(0.99, 0.98), title='Subjects')

    shrink_str = ' (Shrunk ROI -1 voxel)' if use_shrink else ''
    fig.suptitle(
        f'Fixel K-means Direction × Peak Amplitude - {hemisphere.upper()} Hemisphere{shrink_str}\n'
        f'Radar scaled by mean FOD peak amplitude (amplitude ∝ fiber strength) | '
        f'Solid = Dir 1, Dashed = Dir 2',
        fontsize=13, fontweight='bold', y=0.99)

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')

    if show:
        plt.show()

    return fig


def plot_3group_component_lines(all_group_data, hemisphere, target_n_slices,
                                 use_shrink=False, save_path=None, show=False,
                                 y_max=None):
    """Plot per-slice directional component magnitudes for all 3 groups.

    For each group (CTRL, AP, APTR), shows Dir1 component magnitudes (AP, DV, ML)
    scaled by peak amplitude across coronal slices. Each subject is a separate line;
    lines may differ in length (number of slices).

    Parameters
    ----------
    all_group_data : dict
        {group_name: {subject_id: {slice_idx: kmeans_result_dict, ...}, ...}, ...}
    hemisphere : str
    target_n_slices : int  (used only for title; actual slices come from data)
    use_shrink : bool
    save_path : str or Path or None
    show : bool
    y_max : float or None
        Fixed y-axis maximum for consistent cross-plot comparison.
        If None, computed from data.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    row_groups = ['CTRL', 'AP', 'APTR']
    group_titles = ['CTRL (isocortex)', 'AP (graft, no Tx)', 'APTR (graft + Tx)']
    group_bg = ['#e8f5e9', '#ffebee', '#e3f2fd']
    subject_markers = ['o', 's', '^', 'D', 'v']

    # Component colors:  AP=red, DV=blue, ML=green
    comp_colors = {'AP': '#e74c3c', 'DV': '#3498db', 'ML': '#27ae60'}
    comp_idx = {'AP': 0, 'DV': 1, 'ML': 2}  # index into direction vector [AP, DV, ML]

    # Compute y-axis max (use forced value or derive from data)
    if y_max is None:
        global_max = 0.0
        for group in row_groups:
            for subj, data in all_group_data.get(group, {}).items():
                for s_idx, res in data.items():
                    sd = res.get('scaled_directions', res.get('directions', [None]))[0]
                    if sd is not None:
                        val = np.max(np.abs(sd))
                        if val > global_max:
                            global_max = val
        y_max = max(0.05, global_max * 1.15)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=False, sharey=True)

    for row_idx, (group, title, bg) in enumerate(zip(row_groups, group_titles, group_bg)):
        ax = axes[row_idx]
        ax.set_facecolor(bg)
        subjects_data = all_group_data.get(group, {})
        subject_ids = list(subjects_data.keys())

        # Track max slice count for x-axis
        max_slices = 0

        for si, subj in enumerate(subject_ids):
            data = subjects_data[subj]
            slices_sorted = sorted(data.keys())
            n_sl = len(slices_sorted)
            if n_sl > max_slices:
                max_slices = n_sl

            x = np.arange(n_sl)
            marker = subject_markers[si % len(subject_markers)]
            # Short label: last part of subject ID
            short_id = subj.replace('sub-', '')
            if len(short_id) > 12:
                short_id = short_id[:12]

            for comp_name, c_idx in comp_idx.items():
                color = comp_colors[comp_name]
                dirs_key = 'scaled_directions' if 'scaled_directions' in data.get(slices_sorted[0], {}) else 'directions'
                y_vals = []
                for s in slices_sorted:
                    sd = data[s].get(dirs_key, data[s].get('directions', [None]))[0]
                    if sd is not None:
                        y_vals.append(abs(sd[c_idx]))
                    else:
                        y_vals.append(0.0)

                y_vals = np.array(y_vals)
                # Only label first subject's components (for legend clarity)
                label = f'{comp_name}' if si == 0 else None
                alpha = max(0.3, 1.0 - si * 0.15)
                ax.plot(x, y_vals, color=color, marker=marker, markersize=4,
                        linewidth=1.5, alpha=alpha, label=label)

        ax.set_ylabel('Peak × |component|', fontsize=10)
        ax.set_ylim(0, y_max)
        ax.set_title(f'{title}  (n={len(subject_ids)})', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        if max_slices > 0:
            ax.set_xlim(-0.5, max(max_slices - 0.5, target_n_slices - 0.5))
            ax.set_xticks(range(max(max_slices, target_n_slices)))

        if row_idx == 0:
            ax.legend(loc='upper right', fontsize=9, ncol=3,
                      title='Component (Dir1)', title_fontsize=8)

        # Add subject legend as text
        legend_str = '  '.join([f'{subject_markers[i % len(subject_markers)]} {subj.replace("sub-", "")[:15]}'
                                 for i, subj in enumerate(subject_ids)])
        ax.text(0.01, 0.97, legend_str, transform=ax.transAxes,
                fontsize=6, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    axes[-1].set_xlabel('Coronal slice index (anterior → posterior)', fontsize=11)

    shrink_str = ' (Shrunk ROI)' if use_shrink else ''
    fig.suptitle(
        f'Directional Components per Slice - {hemisphere.upper()} Hemisphere{shrink_str}\n'
        f'Dir1 scaled by FOD peak amplitude | AP (red), DV (blue), ML (green)',
        fontsize=13, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')

    if show:
        plt.show()

    return fig


def analyze_3group_kmeans(cfg, hemisphere='left', use_shrink=False,
                           output_dir=None, n_clusters=2, conda_env='mbrain',
                           force=False, show=False, verbose=True,
                           y_max=None, plot=True):
    """Full 3-group K-means comparison pipeline.

    Loads fixels for all subjects in CTRL/AP/APTR, runs per-slice K-means,
    and generates the comparison figure.

    Parameters
    ----------
    cfg : dict
        Loaded config1.yaml
    hemisphere : str  ('left' or 'right')
    use_shrink : bool
        Use shrunk graft ROIs for AP/APTR
    output_dir : str or Path or None
        Where to save the figure. Defaults to REPO_DIR/reports/.
    n_clusters : int
    conda_env : str
    force : bool
    show : bool
    verbose : bool
    y_max : float or None
        Fixed y-axis maximum for consistent cross-plot comparison.
        If None, auto-computed per plot.
    plot : bool
        Whether to generate plots. Set False for data-only pass.

    Returns
    -------
    dict with all_group_data + figure path
    """
    deriv = Path(cfg['paths']['derivatives_flip1'])

    if output_dir is None:
        output_dir = Path(cfg['paths'].get('repo_dir', '.')) / 'reports'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Subject definitions from config (skip test subjects)
    groups_cfg = cfg.get('groups', {})
    subject_map = {}
    for grp_name in ['CTRL', 'AP', 'APTR']:
        subject_map[grp_name] = []
        for entry in groups_cfg.get(grp_name, []):
            sid = entry['id']
            ses = entry['session']
            if 'test' in sid.lower():
                continue
            subject_map[grp_name].append((sid, ses))

    all_group_data = {}
    graft_slice_counts = []

    # Process AP and APTR first (to determine target_n_slices from graft)
    for group in ['AP', 'APTR']:
        if verbose:
            print(f'\n--- {group} ({hemisphere}, shrink={use_shrink}) ---')
        all_group_data[group] = {}

        for (subj, ses) in subject_map[group]:
            subj_dir = deriv / subj / ses
            fixel_dir = subj_dir / 'dwi' / 'fixels_group'
            fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'

            if not fod_file.exists():
                if verbose:
                    print(f'  ⚠ Skipping {subj}: FOD not found')
                continue

            # Ensure fixels exist — reuse if present, never auto-delete
            if (fixel_dir / 'index.mif').exists():
                if verbose:
                    print(f'  ✓ Reusing existing fixels for {subj}')
            else:
                if verbose:
                    print(f'  Generating fixels for {subj}...')
                if fixel_dir.exists() and any(fixel_dir.iterdir()):
                    print(f'  ⚠ Fixel dir not empty but index.mif missing: {fixel_dir}')
                    print(f'    Delete manually to regenerate: rm -rf "{fixel_dir}"')
                    continue
                _run_cmd(
                    f'micromamba run -n {conda_env} fod2fixel "{fod_file}" "{fixel_dir}" '
                    f'-mask "{subj_dir}/dwi/mask_anat.mif" '
                    f'-afd afd.mif -peak peak_amp.mif -disp dispersion.mif '
                    f'-force -quiet', timeout=600)
                if not (fixel_dir / 'index.mif').exists():
                    if verbose:
                        print(f'  ⚠ Skipping {subj}: fixels not generated')
                    continue

            if verbose:
                print(f'  Loading {subj}...')
            fixel_data = load_fixel_dir(fixel_dir, conda_env)

            if use_shrink:
                roi_path = (subj_dir / 'manual_roi' / 'fod_space' / 'shrink' /
                            f'roi_Graft_{hemisphere}_flip_fod_shrink.mif')
            else:
                roi_path = (subj_dir / 'manual_roi' / 'fod_space' /
                            f'roi_Graft_{hemisphere}_flip_fod.mif')

            if not roi_path.exists():
                if verbose:
                    print(f'  ⚠ ROI not found: {roi_path.name}')
                continue

            roi_mask = load_roi_mask(roi_path, conda_env)
            slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(fixel_data, roi_mask, slice_axis=0)
            slice_km = slicewise_kmeans(slice_vecs, n_clusters=n_clusters,
                                        slice_afds=slice_peaks, slice_roi_voxels=slice_voxels)

            if slice_km:
                all_group_data[group][subj] = slice_km
                graft_slice_counts.append(len(slice_km))
                if verbose:
                    print(f'    {len(slice_km)} slices with fixel data')

    if not graft_slice_counts:
        print('  ⚠ No graft fixel data found for AP/APTR!')
        return None

    target_n_slices = min(graft_slice_counts)
    if verbose:
        print(f'\n  Target slice count (smallest graft): {target_n_slices}')

    # Re-index from anterior (take first target_n_slices sorted slices)
    for group in ['AP', 'APTR']:
        for subj, subj_data in all_group_data[group].items():
            sorted_slices = sorted(subj_data.keys())[:target_n_slices]
            all_group_data[group][subj] = {i: subj_data[s]
                                            for i, s in enumerate(sorted_slices)}

    # Process CTRL
    if verbose:
        print(f'\n--- CTRL ({hemisphere}) ---')
    all_group_data['CTRL'] = {}

    for (subj, ses) in subject_map['CTRL']:
        subj_dir = deriv / subj / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'

        if not fod_file.exists():
            if verbose:
                print(f'  ⚠ Skipping {subj}: FOD not found')
            continue

        # Ensure fixels exist — reuse if present, never auto-delete
        if (fixel_dir / 'index.mif').exists():
            if verbose:
                print(f'  ✓ Reusing existing fixels for {subj}')
        else:
            if verbose:
                print(f'  Generating fixels for {subj}...')
            if fixel_dir.exists() and any(fixel_dir.iterdir()):
                print(f'  ⚠ Fixel dir not empty but index.mif missing: {fixel_dir}')
                print(f'    Delete manually to regenerate: rm -rf "{fixel_dir}"')
                continue
            _run_cmd(
                f'micromamba run -n {conda_env} fod2fixel "{fod_file}" "{fixel_dir}" '
                f'-mask "{subj_dir}/dwi/mask_anat.mif" '
                f'-afd afd.mif -peak peak_amp.mif -disp dispersion.mif '
                f'-force -quiet', timeout=600)
            if not (fixel_dir / 'index.mif').exists():
                if verbose:
                    print(f'  ⚠ Skipping {subj}: fixels not generated')
                continue

        if verbose:
            print(f'  Loading {subj}...')
        fixel_data = load_fixel_dir(fixel_dir, conda_env)

        roi_path = (subj_dir / 'atlas' / 'rois' / 'fod_space' /
                    f'roi_isocortex_{hemisphere}.nii.gz')
        if not roi_path.exists():
            if verbose:
                print(f'  ⚠ ROI not found: {roi_path.name}')
            continue

        roi_mask = load_roi_mask(roi_path, conda_env)
        slice_vecs, slice_peaks, slice_voxels = extract_slicewise_vectors(fixel_data, roi_mask, slice_axis=0)
        slice_km = slicewise_kmeans(slice_vecs, n_clusters=n_clusters,
                                    slice_afds=slice_peaks, slice_roi_voxels=slice_voxels)

        if slice_km:
            sorted_slices = sorted(slice_km.keys())[:target_n_slices]
            all_group_data['CTRL'][subj] = {i: slice_km[s]
                                             for i, s in enumerate(sorted_slices)}
            if verbose:
                print(f'    {len(sorted_slices)} slices with fixel data')

    suffix = '_shrink' if use_shrink else ''
    result = dict(all_group_data=all_group_data,
                  target_n_slices=target_n_slices)

    if not plot:
        return result

    # Generate radar figure
    fig_path = output_dir / f'fixel_3group_comparison_{hemisphere}{suffix}.png'
    fig = plot_3group_kmeans_comparison(
        all_group_data, hemisphere, target_n_slices,
        use_shrink=use_shrink, save_path=fig_path, show=show,
        y_max=y_max)

    # Generate per-slice component line plot
    lines_path = output_dir / f'fixel_3group_components_{hemisphere}{suffix}.png'
    fig_lines = plot_3group_component_lines(
        all_group_data, hemisphere, target_n_slices,
        use_shrink=use_shrink, save_path=lines_path, show=show,
        y_max=y_max)

    result.update(figure=fig, figure_path=fig_path,
                  figure_lines=fig_lines, lines_path=lines_path)
    return result
