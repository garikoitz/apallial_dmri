#!/usr/bin/env python3
"""
FOD Direction Analysis - Graft ROI
Analyzes FOD directionality within graft regions for APTR (and AP) mice.

Creates 3D orientation plots showing mean FOD orientations per A-P slice,
with proper handling of the full 3D orientation space.

Axes convention (in FOD/voxel space after flip1):
  - Axis 0 (X): Anterior-Posterior (A→P)
  - Axis 1 (Y): Dorsal-Ventral (D→V)
  - Axis 2 (Z): Left-Right (L→R)

Usage:
    python functions/fod_direction_analysis.py [--group APTR|AP]
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
import os
import sys
import subprocess
import tempfile
import argparse

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions.config_loader import load_config


def get_sphere_vertices(n_points=100):
    """Generate evenly distributed points on a unit sphere using Fibonacci spiral."""
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + np.sqrt(5)) * indices
    
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    
    return np.column_stack([x, y, z])


def compute_sh_basis(directions, lmax=8):
    """
    Compute real spherical harmonic basis for given directions.
    For lmax=8, we have 45 coefficients.
    """
    from scipy.special import sph_harm
    
    x, y, z = directions[:, 0], directions[:, 1], directions[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(np.clip(z / (r + 1e-10), -1, 1))
    phi = np.arctan2(y, x)
    
    n_dirs = len(directions)
    n_coeffs = int((lmax + 1) * (lmax + 2) / 2)
    
    basis = np.zeros((n_dirs, n_coeffs))
    idx = 0
    for l in range(0, lmax + 1, 2):
        for m in range(-l, l + 1):
            if m < 0:
                Y = np.sqrt(2) * (-1)**m * np.imag(sph_harm(abs(m), l, phi, theta))
            elif m == 0:
                Y = np.real(sph_harm(0, l, phi, theta))
            else:
                Y = np.sqrt(2) * (-1)**m * np.real(sph_harm(m, l, phi, theta))
            
            basis[:, idx] = np.real(Y)
            idx += 1
            if idx >= n_coeffs:
                break
        if idx >= n_coeffs:
            break
    
    return basis


def load_mif(filepath):
    """Load a .mif file and return data."""
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    
    subprocess.run(['mrconvert', filepath, tmp_path, '-force', '-quiet'], check=True)
    img = nib.load(tmp_path)
    data = img.get_fdata()
    os.remove(tmp_path)
    
    return data


def create_shrunk_graft_rois(config):
    """
    Create shrunk (eroded by 1 voxel) versions of graft ROIs for AP and APTR groups.
    Saves to manual_roi/fod_space/shrink/ folder.
    
    Uses scipy.ndimage.binary_erosion with a 3D ball structuring element.
    """
    from scipy import ndimage
    
    deriv = config['paths']['derivatives_flip1']
    
    # All subjects with graft ROIs
    all_subjects = {
        'AP': [
            ('sub-596-C-OO-NmM', 'ses-240911103349'),
            ('sub-769-C-O5-R1M', 'ses-240911115733'),
            ('sub-769-C-OO-NmM', 'ses-240911111853'),
            ('sub-875-C-OO-NM', 'ses-240911100253'),
        ],
        'APTR': [
            ('sub-AP72-862-901-NmF', 'ses-251210093342'),
            ('sub-AP72-862-901-R1M', 'ses-251210100043'),
            ('sub-Apallial17NmF', 'ses-240816091633'),
            ('sub-Apallial20R1F', 'ses-240909094232'),
        ],
    }
    
    # 3D ball structuring element (radius 1 = 6-connectivity erosion)
    struct = ndimage.generate_binary_structure(3, 1)
    
    for group, subjects in all_subjects.items():
        print(f"\n--- Creating shrunk ROIs for {group} ---")
        
        for subj, ses in subjects:
            base = os.path.join(deriv, subj, ses, 'manual_roi', 'fod_space')
            shrink_dir = os.path.join(base, 'shrink')
            os.makedirs(shrink_dir, exist_ok=True)
            
            for hemi in ['left', 'right']:
                roi_path = os.path.join(base, f'roi_Graft_{hemi}_flip_fod.mif')
                
                if not os.path.exists(roi_path):
                    print(f"  WARNING: {roi_path} not found, skipping")
                    continue
                
                # Load ROI
                roi_data = load_mif(roi_path) > 0.5
                original_voxels = roi_data.sum()
                
                # Erode by 1 voxel
                eroded = ndimage.binary_erosion(roi_data, structure=struct)
                eroded_voxels = eroded.sum()
                
                # Save as .mif
                out_path = os.path.join(shrink_dir, f'roi_Graft_{hemi}_flip_fod_shrink.mif')
                
                # Create temporary nifti then convert to mif
                with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
                    tmp_nii = tmp.name
                
                # Need to get affine from original - load with nibabel via mrconvert
                with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp2:
                    tmp_orig = tmp2.name
                subprocess.run(['mrconvert', roi_path, tmp_orig, '-force', '-quiet'], check=True)
                orig_img = nib.load(tmp_orig)
                
                # Save eroded mask
                eroded_img = nib.Nifti1Image(eroded.astype(np.float32), orig_img.affine, orig_img.header)
                nib.save(eroded_img, tmp_nii)
                
                # Convert to mif
                subprocess.run(['mrconvert', tmp_nii, out_path, '-force', '-quiet'], check=True)
                
                # Cleanup
                os.remove(tmp_nii)
                os.remove(tmp_orig)
                
                print(f"  {subj} {hemi}: {original_voxels} -> {eroded_voxels} voxels "
                      f"({100*eroded_voxels/original_voxels:.1f}% retained)")
    
    print("\nShrunk ROIs created successfully!")


def analyze_graft_fod(subject, session, hemisphere, config, roi_type='graft', use_shrink=False):
    """
    Analyze FOD directions within graft/ROI for one subject/hemisphere.
    
    Parameters:
    -----------
    subject : str
    session : str
    hemisphere : str ('left' or 'right')
    config : dict
    roi_type : str ('graft' for AP/APTR, 'isocortex' for CTRL - though for now using graft)
    use_shrink : bool - if True, use the eroded (shrunk) version of the graft ROI
    
    Returns:
    --------
    slice_data : dict with slice_idx -> {n_voxels, ap_sum, dv_sum, lr_sum, ap_fraction, dv_fraction, lr_fraction, total_amplitude}
    """
    deriv = config['paths']['derivatives_flip1']
    base = os.path.join(deriv, subject, session)
    
    # Load FOD
    fod_path = os.path.join(base, 'dwi', 'wm_fod_norm.mif')
    fod_data = load_mif(fod_path)
    print(f"  FOD shape: {fod_data.shape}")
    
    # Load graft mask (shrunk or original)
    if use_shrink:
        mask_path = os.path.join(base, 'manual_roi', 'fod_space', 'shrink',
                                 f'roi_Graft_{hemisphere}_flip_fod_shrink.mif')
        print(f"  Using SHRUNK mask")
    else:
        mask_path = os.path.join(base, 'manual_roi', 'fod_space', 
                                 f'roi_Graft_{hemisphere}_flip_fod.mif')
    mask_data = load_mif(mask_path) > 0.5
    print(f"  Mask voxels: {mask_data.sum()}")
    
    # Find A-P slice range (axis 0)
    mask_slices = np.any(mask_data, axis=(1, 2))
    slice_indices = np.where(mask_slices)[0]
    print(f"  Graft spans slices: {slice_indices.min()} to {slice_indices.max()} ({len(slice_indices)} slices)")
    
    # Precompute SH basis for sampling
    n_sample = 100
    sample_dirs = get_sphere_vertices(n_sample)
    sh_basis = compute_sh_basis(sample_dirs, lmax=8)
    
    # Analyze each slice
    slice_data = {}
    
    for slice_idx in slice_indices:
        slice_mask = mask_data[slice_idx, :, :]
        slice_fod = fod_data[slice_idx, :, :, :]
        
        y_coords, z_coords = np.where(slice_mask)
        
        # Accumulate weighted sums for each direction (like your example code)
        total_ap = 0.0  # Anterior-Posterior
        total_dv = 0.0  # Dorsal-Ventral
        total_lr = 0.0  # Left-Right
        total_wm = 0.0  # WM content from c[0]
        n_valid_voxels = 0
        
        for y, z in zip(y_coords, z_coords):
            sh_coeffs = slice_fod[y, z, :]
            
            if np.max(np.abs(sh_coeffs)) < 0.01:
                continue
            
            # WM content from c[0] (isotropic component)
            total_wm += sh_coeffs[0]
            
            # Sample FOD at many directions
            amplitudes = sh_basis @ sh_coeffs
            amplitudes = np.maximum(amplitudes, 0)  # FOD should be non-negative
            
            # For each sampled direction, weight by amplitude
            # This mimics: sum(|component| * magnitude) from the peaks approach
            # Here, amplitude IS the magnitude in each direction
            for d, amp in zip(sample_dirs, amplitudes):
                if amp > 0.01:  # Threshold small values
                    # |component| × amplitude (like your formula)
                    total_ap += np.abs(d[0]) * amp  # A-P component
                    total_dv += np.abs(d[1]) * amp  # D-V component
                    total_lr += np.abs(d[2]) * amp  # L-R component
            
            n_valid_voxels += 1
        
        if n_valid_voxels > 0:
            # Store raw totals (like your approach)
            total_all = total_ap + total_dv + total_lr
            
            slice_data[slice_idx] = {
                'n_voxels': n_valid_voxels,
                # Raw weighted sums (for comparison across slices)
                'ap_sum': total_ap,
                'dv_sum': total_dv,
                'lr_sum': total_lr,
                # Fractions (for radar plots)
                'ap_fraction': total_ap / total_all if total_all > 0 else 0,
                'dv_fraction': total_dv / total_all if total_all > 0 else 0,
                'lr_fraction': total_lr / total_all if total_all > 0 else 0,
                'total_amplitude': total_all,
                # WM content
                'wm_content': total_wm
            }
    
    return slice_data


def analyze_isocortex_fod(subject, session, hemisphere, config):
    """
    Analyze FOD directions within isocortex ROI for CTRL subjects.
    Uses the Allen atlas isocortex ROI (already in fod_space).
    
    Parameters:
    -----------
    subject : str
    session : str
    hemisphere : str ('left' or 'right')
    config : dict
    
    Returns:
    --------
    slice_data : dict with slice_idx -> same structure as analyze_graft_fod
    """
    deriv = config['paths']['derivatives_flip1']
    base = os.path.join(deriv, subject, session)
    
    # Load FOD
    fod_path = os.path.join(base, 'dwi', 'wm_fod_norm.mif')
    fod_data = load_mif(fod_path)
    print(f"  FOD shape: {fod_data.shape}")
    
    # Load isocortex mask from atlas (already in fod_space)
    mask_path = os.path.join(base, 'atlas', 'rois', 'fod_space', 
                             f'roi_isocortex_{hemisphere}.nii.gz')
    mask_img = nib.load(mask_path)
    mask_data = mask_img.get_fdata() > 0.5
    print(f"  Isocortex voxels: {mask_data.sum()}")
    
    # Find A-P slice range (axis 0)
    mask_slices = np.any(mask_data, axis=(1, 2))
    slice_indices = np.where(mask_slices)[0]
    print(f"  Isocortex spans slices: {slice_indices.min()} to {slice_indices.max()} ({len(slice_indices)} slices)")
    
    # Precompute SH basis for sampling
    n_sample = 100
    sample_dirs = get_sphere_vertices(n_sample)
    sh_basis = compute_sh_basis(sample_dirs, lmax=8)
    
    # Analyze each slice
    slice_data = {}
    
    for slice_idx in slice_indices:
        slice_mask = mask_data[slice_idx, :, :]
        slice_fod = fod_data[slice_idx, :, :, :]
        
        y_coords, z_coords = np.where(slice_mask)
        
        total_ap = 0.0
        total_dv = 0.0
        total_lr = 0.0
        total_wm = 0.0
        n_valid_voxels = 0
        
        for y, z in zip(y_coords, z_coords):
            sh_coeffs = slice_fod[y, z, :]
            
            if np.max(np.abs(sh_coeffs)) < 0.01:
                continue
            
            total_wm += sh_coeffs[0]
            
            amplitudes = sh_basis @ sh_coeffs
            amplitudes = np.maximum(amplitudes, 0)
            
            for d, amp in zip(sample_dirs, amplitudes):
                if amp > 0.01:
                    total_ap += np.abs(d[0]) * amp
                    total_dv += np.abs(d[1]) * amp
                    total_lr += np.abs(d[2]) * amp
            
            n_valid_voxels += 1
        
        if n_valid_voxels > 0:
            total_all = total_ap + total_dv + total_lr
            
            slice_data[slice_idx] = {
                'n_voxels': n_valid_voxels,
                'ap_sum': total_ap,
                'dv_sum': total_dv,
                'lr_sum': total_lr,
                'ap_fraction': total_ap / total_all if total_all > 0 else 0,
                'dv_fraction': total_dv / total_all if total_all > 0 else 0,
                'lr_fraction': total_lr / total_all if total_all > 0 else 0,
                'total_amplitude': total_all,
                'wm_content': total_wm
            }
    
    return slice_data


def analyze_wholebrain_fod(subject, session, config, target_n_slices=None):
    """
    Analyze FOD directions across whole brain mask for CTRL subjects.
    If target_n_slices is specified, average slices to match that count.
    
    Parameters:
    -----------
    subject : str
    session : str
    config : dict
    target_n_slices : int, optional - target number of slices (for matching graft slice count)
    
    Returns:
    --------
    slice_data : dict with slice_idx -> same structure as analyze_graft_fod
    """
    deriv = config['paths']['derivatives_flip1']
    base = os.path.join(deriv, subject, session)
    
    # Load FOD
    fod_path = os.path.join(base, 'dwi', 'wm_fod_norm.mif')
    fod_data = load_mif(fod_path)
    print(f"  FOD shape: {fod_data.shape}")
    
    # Load whole brain mask (dwi/mask_anat.mif)
    mask_path = os.path.join(base, 'dwi', 'mask_anat.mif')
    mask_data = load_mif(mask_path) > 0.5
    print(f"  Mask voxels: {mask_data.sum()}")
    
    # Find A-P slice range (axis 0)
    mask_slices = np.any(mask_data, axis=(1, 2))
    slice_indices = np.where(mask_slices)[0]
    print(f"  Brain spans slices: {slice_indices.min()} to {slice_indices.max()} ({len(slice_indices)} slices)")
    
    # Precompute SH basis for sampling
    n_sample = 100
    sample_dirs = get_sphere_vertices(n_sample)
    sh_basis = compute_sh_basis(sample_dirs, lmax=8)
    
    # Analyze each slice
    raw_slice_data = {}
    
    for slice_idx in slice_indices:
        slice_mask = mask_data[slice_idx, :, :]
        slice_fod = fod_data[slice_idx, :, :, :]
        
        y_coords, z_coords = np.where(slice_mask)
        
        # Accumulate weighted sums
        total_ap = 0.0
        total_dv = 0.0
        total_lr = 0.0
        total_wm = 0.0  # c[0] coefficient - WM content
        n_valid_voxels = 0
        
        for y, z in zip(y_coords, z_coords):
            sh_coeffs = slice_fod[y, z, :]
            
            if np.max(np.abs(sh_coeffs)) < 0.01:
                continue
            
            # WM content from c[0]
            total_wm += sh_coeffs[0]
            
            amplitudes = sh_basis @ sh_coeffs
            amplitudes = np.maximum(amplitudes, 0)
            
            for d, amp in zip(sample_dirs, amplitudes):
                if amp > 0.01:
                    total_ap += np.abs(d[0]) * amp
                    total_dv += np.abs(d[1]) * amp
                    total_lr += np.abs(d[2]) * amp
            
            n_valid_voxels += 1
        
        if n_valid_voxels > 0:
            total_all = total_ap + total_dv + total_lr
            
            raw_slice_data[slice_idx] = {
                'n_voxels': n_valid_voxels,
                'ap_sum': total_ap,
                'dv_sum': total_dv,
                'lr_sum': total_lr,
                'ap_fraction': total_ap / total_all if total_all > 0 else 0,
                'dv_fraction': total_dv / total_all if total_all > 0 else 0,
                'lr_fraction': total_lr / total_all if total_all > 0 else 0,
                'total_amplitude': total_all,
                'wm_content': total_wm
            }
    
    # If target_n_slices specified, bin/average to match
    if target_n_slices is not None and len(raw_slice_data) > target_n_slices:
        n_raw = len(raw_slice_data)
        bin_size = n_raw / target_n_slices
        sorted_slices = sorted(raw_slice_data.keys())
        
        slice_data = {}
        for i in range(target_n_slices):
            start_idx = int(i * bin_size)
            end_idx = int((i + 1) * bin_size)
            bin_slices = sorted_slices[start_idx:end_idx]
            
            # Average the bin
            n_vox = sum(raw_slice_data[s]['n_voxels'] for s in bin_slices)
            ap_s = sum(raw_slice_data[s]['ap_sum'] for s in bin_slices)
            dv_s = sum(raw_slice_data[s]['dv_sum'] for s in bin_slices)
            lr_s = sum(raw_slice_data[s]['lr_sum'] for s in bin_slices)
            wm_c = sum(raw_slice_data[s]['wm_content'] for s in bin_slices)
            total_all = ap_s + dv_s + lr_s
            
            slice_data[i] = {
                'n_voxels': n_vox,
                'ap_sum': ap_s,
                'dv_sum': dv_s,
                'lr_sum': lr_s,
                'ap_fraction': ap_s / total_all if total_all > 0 else 0,
                'dv_fraction': dv_s / total_all if total_all > 0 else 0,
                'lr_fraction': lr_s / total_all if total_all > 0 else 0,
                'total_amplitude': total_all,
                'wm_content': wm_c
            }
        
        return slice_data
    
    return raw_slice_data


def plot_polar_histogram(ax, directions, magnitudes, color, alpha=0.4, n_bins=18):
    """
    Plot a semicircular histogram of ORIENTATIONS in the D-V vs L-R plane.
    Only shows [0, π) since orientations are antipodally symmetric.
    """
    # Project to Y (D-V) vs Z (L-R) plane
    angles = np.arctan2(directions[:, 1], directions[:, 2])
    
    # Fold to [0, π) - orientations, not directions
    angles_folded = np.where(angles < 0, angles + np.pi, angles)
    
    # Histogram over semicircle [0, π)
    bars, bin_edges = np.histogram(angles_folded, bins=n_bins, weights=magnitudes, range=(0, np.pi))
    
    if bars.max() > 0:
        bars = bars / bars.max()
    
    width = np.pi / n_bins
    ax.bar(bin_edges[:-1], bars, width=width, color=color, alpha=alpha, edgecolor='none')
    
    return bars


def plot_mean_arrow(ax, mean_direction, mean_magnitude, color, max_mag=None):
    """Plot arrow showing the mean ORIENTATION."""
    angle = np.arctan2(mean_direction[1], mean_direction[2])
    # Fold to [0, π)
    if angle < 0:
        angle += np.pi
    
    if max_mag is not None and max_mag > 0:
        r = 0.8 * (mean_magnitude / max_mag)
    else:
        r = 0.5
    r = min(r, 0.9)
    
    ax.annotate('', xy=(angle, r), xytext=(angle, 0),
                arrowprops=dict(arrowstyle='->', color=color, lw=2))


def orientation_to_spherical(direction):
    """
    Convert 3D orientation vector to spherical coordinates.
    
    Returns:
        elevation: angle from horizontal plane (D-V), -90 to +90 degrees
        azimuth: angle in horizontal plane from A-P axis, 0 to 180 degrees (orientation)
    
    Convention:
        - elevation = 0: horizontal (in A-P / L-R plane)
        - elevation = +90: dorsal
        - azimuth = 0: anterior-posterior axis
        - azimuth = 90: left-right axis
    """
    # Ensure we use consistent orientation (fold antipodal)
    if direction[1] < 0:  # If pointing ventral, flip
        direction = -direction
    
    x, y, z = direction[0], direction[1], direction[2]  # A-P, D-V, L-R
    
    # Horizontal distance
    r_horiz = np.sqrt(x**2 + z**2)
    
    # Elevation (angle from horizontal)
    elevation = np.degrees(np.arctan2(y, r_horiz))
    
    # Azimuth (angle in horizontal plane from A-P axis)
    azimuth = np.degrees(np.arctan2(abs(z), abs(x)))  # 0-90 range
    
    return elevation, azimuth


def plot_3d_orientations(subject_data, hemisphere, group, config, output_dir):
    """
    Create intuitive visualizations of mean orientations.
    Main figure: Radar plots per slice showing FOD amplitude in each direction.
    """
    all_slices = set()
    for subj_info in subject_data.values():
        all_slices.update(subj_info['data'].keys())
    slices_sorted = sorted(all_slices)
    n_slices = len(slices_sorted)
    
    # ===== FIGURE 1: Radar plot per slice =====
    n_cols = min(6, n_slices)
    n_rows = int(np.ceil(n_slices / n_cols))
    
    fig = plt.figure(figsize=(4 * n_cols, 4 * n_rows))
    
    # Order like your example: DV -> AP -> LR (with DV at South/bottom)
    categories = ['Dor-Ven', 'Ant-Post', 'Left-Right']
    n_cats = len(categories)
    angles = [n / float(n_cats) * 2 * np.pi for n in range(n_cats)]
    angles += angles[:1]  # Close the polygon
    
    # Find max value for consistent scaling
    max_val = 0
    for subj_info in subject_data.values():
        for sd in subj_info['data'].values():
            max_val = max(max_val, sd['dv_fraction'], sd['ap_fraction'], sd['lr_fraction'])
    
    for i, slice_idx in enumerate(slices_sorted):
        ax = fig.add_subplot(n_rows, n_cols, i + 1, projection='polar')
        
        for subj, subj_info in subject_data.items():
            data = subj_info['data']
            color = subj_info['color']
            
            if slice_idx in data:
                sd = data[slice_idx]
                # Order: DV, AP, LR (matching categories)
                values = [sd['dv_fraction'], sd['ap_fraction'], sd['lr_fraction']]
                values += values[:1]  # Close polygon
                
                ax.plot(angles, values, 'o-', linewidth=2, color=color, 
                       markersize=6, alpha=0.8)
                ax.fill(angles, values, alpha=0.15, color=color)
        
        # DV at South (bottom), like your example
        ax.set_theta_zero_location("S")
        ax.set_theta_direction(1)  # Counter-clockwise
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=9, fontweight='bold')
        ax.set_ylim(0, max_val * 1.1)
        ax.set_yticks([])  # Hide radial ticks like your example
        ax.set_title(f'Slice {slice_idx}', fontsize=10, fontweight='bold', pad=15)
    
    # Legend
    legend_elements = [plt.Line2D([0], [0], color=info['color'], linewidth=2, 
                                   label=subj.replace('sub-', '')[:15])
                      for subj, info in subject_data.items()]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=9,
              bbox_to_anchor=(0.99, 0.99))
    
    fig.suptitle(f'{group} {hemisphere.upper()} Graft - FOD Orientation per Slice\n'
                 f'(Radar: fraction of total FOD amplitude in each direction)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])
    
    output_path = os.path.join(output_dir, f'fod_radar_{group.lower()}_{hemisphere}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()
    
    # ===== FIGURE 2: Heatmap + summary =====
    fig2 = plt.figure(figsize=(16, 10))
    
    subjects = list(subject_data.keys())
    n_subj = len(subjects)
    
    # Panel 1: RGB heatmap
    ax1 = fig2.add_subplot(2, 2, 1)
    rgb_matrix = np.zeros((n_slices, n_subj, 3))
    
    for j, (subj, subj_info) in enumerate(subject_data.items()):
        data = subj_info['data']
        for i, slice_idx in enumerate(slices_sorted):
            if slice_idx in data:
                sd = data[slice_idx]
                rgb_matrix[i, j, 0] = sd['lr_fraction']  # L-R -> Red
                rgb_matrix[i, j, 1] = sd['ap_fraction']  # A-P -> Green
                rgb_matrix[i, j, 2] = sd['dv_fraction']  # D-V -> Blue
            else:
                rgb_matrix[i, j, :] = 0.3
    
    ax1.imshow(rgb_matrix, aspect='auto', interpolation='nearest')
    ax1.set_yticks(range(n_slices))
    ax1.set_yticklabels([str(s) for s in slices_sorted], fontsize=8)
    ax1.set_xticks(range(n_subj))
    ax1.set_xticklabels([s.replace('sub-', '')[:12] for s in subjects], 
                        rotation=45, ha='right', fontsize=8)
    ax1.set_ylabel('A-P Slice')
    ax1.set_xlabel('Subject')
    ax1.set_title('FOD Orientation RGB Map\n(Red=L-R, Green=A-P, Blue=D-V)', fontsize=11)
    
    # Panel 2: Line plot per component
    ax2 = fig2.add_subplot(2, 2, 2)
    
    for subj, subj_info in subject_data.items():
        data = subj_info['data']
        color = subj_info['color']
        
        slices = []
        ap_vals = []
        
        for slice_idx in slices_sorted:
            if slice_idx in data:
                slices.append(slice_idx)
                ap_vals.append(data[slice_idx]['ap_fraction'])
        
        ax2.plot(slices, ap_vals, 'o-', color=color, linewidth=2, alpha=0.8,
                label=subj.replace('sub-', '')[:12], markersize=4)
    
    ax2.set_xlabel('A-P Slice', fontsize=10)
    ax2.set_ylabel('A-P Fraction (FOD amplitude weighted)', fontsize=10)
    ax2.set_title('Anterior-Posterior Component per Slice', fontsize=11)
    ax2.set_ylim(0, 0.7)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(1/3, color='gray', linestyle='--', alpha=0.5, label='Isotropic')
    ax2.legend(loc='upper right', fontsize=7)
    
    # Panel 3: Stacked area (average across subjects)
    ax3 = fig2.add_subplot(2, 2, 3)
    
    avg_ap, avg_dv, avg_lr = [], [], []
    for slice_idx in slices_sorted:
        ap_sum, dv_sum, lr_sum, count = 0, 0, 0, 0
        for subj_info in subject_data.values():
            if slice_idx in subj_info['data']:
                sd = subj_info['data'][slice_idx]
                ap_sum += sd['ap_fraction']
                dv_sum += sd['dv_fraction']
                lr_sum += sd['lr_fraction']
                count += 1
        
        if count > 0:
            avg_ap.append(ap_sum / count)
            avg_dv.append(dv_sum / count)
            avg_lr.append(lr_sum / count)
        else:
            avg_ap.append(0)
            avg_dv.append(0)
            avg_lr.append(0)
    
    ax3.stackplot(slices_sorted, avg_ap, avg_dv, avg_lr,
                  labels=['A-P (along cord)', 'D-V (dorsal-ventral)', 'L-R (crossing)'],
                  colors=['#2ecc71', '#3498db', '#e74c3c'],
                  alpha=0.8)
    ax3.set_xlabel('A-P Slice', fontsize=10)
    ax3.set_ylabel('Fraction', fontsize=10)
    ax3.set_title('Mean FOD Composition per Slice\n(averaged across subjects)', fontsize=11)
    ax3.legend(loc='upper right', fontsize=8)
    ax3.set_ylim(0, 1)
    
    # Panel 4: Overall radar per subject
    ax4 = fig2.add_subplot(2, 2, 4, projection='polar')
    
    # Same order as per-slice plots: DV, AP, LR
    categories_radar = ['Dor-Ven', 'Ant-Post', 'Left-Right']
    angles_radar = [n / 3.0 * 2 * np.pi for n in range(3)]
    angles_radar += angles_radar[:1]
    
    for subj, subj_info in subject_data.items():
        data = subj_info['data']
        color = subj_info['color']
        
        # Sum across slices (like your approach - raw sums, not fractions)
        total_ap, total_dv, total_lr = 0, 0, 0
        for sd in data.values():
            total_ap += sd['ap_sum']
            total_dv += sd['dv_sum']
            total_lr += sd['lr_sum']
        
        total_all = total_ap + total_dv + total_lr
        if total_all > 0:
            # Order: DV, AP, LR
            values = [total_dv/total_all, total_ap/total_all, total_lr/total_all]
        else:
            values = [0, 0, 0]
        values += values[:1]
        
        ax4.plot(angles_radar, values, 'o-', linewidth=2, color=color, 
                label=subj.replace('sub-', '')[:12], markersize=8)
        ax4.fill(angles_radar, values, alpha=0.15, color=color)
    
    ax4.set_theta_zero_location("S")
    ax4.set_theta_direction(1)
    ax4.set_xticks(angles_radar[:-1])
    ax4.set_xticklabels(categories_radar, fontsize=10, fontweight='bold')
    ax4.set_ylim(0, 0.6)
    ax4.set_yticks([])
    ax4.set_title('Overall FOD Profile\n(sum across all slices)', fontsize=11)
    ax4.legend(loc='upper right', bbox_to_anchor=(1.35, 1), fontsize=8)
    
    fig2.suptitle(f'{group} {hemisphere.upper()} Graft - FOD Orientation Summary', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path2 = os.path.join(output_dir, f'fod_3d_{group.lower()}_{hemisphere}.png')
    plt.savefig(output_path2, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path2}")
    plt.close()
    
    return output_path


def plot_composite_summary(subject_data, hemisphere, group, config, output_dir):
    """
    Create a clean summary figure with simple bar charts.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    subjects = list(subject_data.keys())
    n_subj = len(subjects)
    x_pos = np.arange(n_subj)
    labels = [s.replace('sub-', '')[:12] for s in subjects]
    colors_subj = [subject_data[s]['color'] for s in subjects]
    
    # Collect overall sums per subject (using the new data structure)
    mean_ap, mean_dv, mean_lr = [], [], []
    
    for subj, subj_info in subject_data.items():
        data = subj_info['data']
        
        # Sum across all slices
        total_ap = sum(sd['ap_sum'] for sd in data.values())
        total_dv = sum(sd['dv_sum'] for sd in data.values())
        total_lr = sum(sd['lr_sum'] for sd in data.values())
        
        total_all = total_ap + total_dv + total_lr
        if total_all > 0:
            mean_ap.append(total_ap / total_all)
            mean_dv.append(total_dv / total_all)
            mean_lr.append(total_lr / total_all)
        else:
            mean_ap.append(0)
            mean_dv.append(0)
            mean_lr.append(0)
    
    # Plot 1: Grouped bars
    ax1 = axes[0]
    width = 0.25
    ax1.bar(x_pos - width, mean_ap, width, label='A-P (along cord)', color='#2ecc71', edgecolor='black')
    ax1.bar(x_pos, mean_dv, width, label='D-V (dorsal-ventral)', color='#3498db', edgecolor='black')
    ax1.bar(x_pos + width, mean_lr, width, label='L-R (crossing)', color='#e74c3c', edgecolor='black')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax1.set_ylabel('Orientation Component', fontsize=10)
    ax1.set_title('Directional Components per Subject', fontsize=11)
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_ylim(0, 1)
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Horizontal bars showing dominant direction
    ax2 = axes[1]
    
    # Calculate dominance ratio: how much the top component dominates
    dominant_names = []
    dominant_ratios = []
    
    for i, (ap, dv, lr) in enumerate(zip(mean_ap, mean_dv, mean_lr)):
        components = {'A-P': ap, 'D-V': dv, 'L-R': lr}
        sorted_comp = sorted(components.items(), key=lambda x: x[1], reverse=True)
        dominant_names.append(sorted_comp[0][0])
        # Ratio: top component / second component
        if sorted_comp[1][1] > 0:
            dominant_ratios.append(sorted_comp[0][1] / sorted_comp[1][1])
        else:
            dominant_ratios.append(3.0)  # Cap at 3
    
    bar_colors = {'A-P': '#2ecc71', 'D-V': '#3498db', 'L-R': '#e74c3c'}
    for i, (name, ratio, label) in enumerate(zip(dominant_names, dominant_ratios, labels)):
        ax2.barh(i, min(ratio, 3), color=bar_colors[name], edgecolor='black', alpha=0.8)
        ax2.text(min(ratio, 3) + 0.1, i, f'{name}', va='center', fontsize=9)
    
    ax2.set_yticks(range(n_subj))
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.set_xlabel('Dominance Ratio (1st / 2nd component)', fontsize=10)
    ax2.set_title('Dominant Orientation', fontsize=11)
    ax2.axvline(1, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlim(0, 3.5)
    
    # Plot 3: Summary text box
    ax3 = axes[2]
    ax3.axis('off')
    
    summary_text = f"SUMMARY: {group} {hemisphere.upper()} Graft\n"
    summary_text += "=" * 35 + "\n\n"
    
    for i, (subj, ap, dv, lr) in enumerate(zip(labels, mean_ap, mean_dv, mean_lr)):
        components = {'A-P': ap, 'D-V': dv, 'L-R': lr}
        dominant = max(components, key=components.get)
        pct = components[dominant] / (ap + dv + lr) * 100
        summary_text += f"{subj}:\n"
        summary_text += f"  Dominant: {dominant} ({pct:.0f}%)\n"
        summary_text += f"  A-P={ap:.2f}, D-V={dv:.2f}, L-R={lr:.2f}\n\n"
    
    # Group average
    avg_ap = np.mean(mean_ap)
    avg_dv = np.mean(mean_dv)
    avg_lr = np.mean(mean_lr)
    group_dom = max({'A-P': avg_ap, 'D-V': avg_dv, 'L-R': avg_lr}, 
                    key=lambda k: {'A-P': avg_ap, 'D-V': avg_dv, 'L-R': avg_lr}[k])
    
    summary_text += "-" * 35 + "\n"
    summary_text += f"GROUP MEAN:\n"
    summary_text += f"  A-P={avg_ap:.2f}, D-V={avg_dv:.2f}, L-R={avg_lr:.2f}\n"
    summary_text += f"  Dominant: {group_dom}\n"
    
    ax3.text(0.1, 0.95, summary_text, transform=ax3.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle(f'{group} {hemisphere.upper()} - Overall Orientation Summary', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, f'fod_summary_{group.lower()}_{hemisphere}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()
    
    return output_path


def analyze_all_groups_comparison(hemisphere, config, output_dir, use_shrink=False):
    """
    Create a 3-row comparison figure: CTRL, AP, APTR
    Each row shows radar plots per slice with WM content circles in center.
    Slices aligned from anterior, ends when smallest graft ends.
    
    Parameters:
    -----------
    use_shrink : bool - if True, use eroded (shrunk by 1 voxel) graft ROIs for AP/APTR
    """
    shrink_label = " (SHRUNK ROIs)" if use_shrink else ""
    print(f"\n{'='*60}")
    print(f"3-GROUP COMPARISON - {hemisphere.upper()} hemisphere{shrink_label}")
    print(f"{'='*60}")
    
    # Subject definitions for all groups
    all_subjects = {
        'CTRL': [
            ('sub-129-921-R1F', 'ses-260115121231'),
            ('sub-129-921-R1L1F', 'ses-260116085433'),
            ('sub-129-1014-L1M', 'ses-260115113814'),
            ('sub-138-1010-L1F', 'ses-260115101441'),
        ],
        'AP': [
            ('sub-596-C-OO-NmM', 'ses-240911103349'),
            ('sub-769-C-O5-R1M', 'ses-240911115733'),
            ('sub-769-C-OO-NmM', 'ses-240911111853'),
            ('sub-875-C-OO-NM', 'ses-240911100253'),
        ],
        'APTR': [
            ('sub-AP72-862-901-NmF', 'ses-251210093342'),
            ('sub-AP72-862-901-R1M', 'ses-251210100043'),
            ('sub-Apallial17NmF', 'ses-240816091633'),
            ('sub-Apallial20R1F', 'ses-240909094232'),
        ],
    }
    
    colors = ['#e74c3c', '#3498db', '#27ae60', '#9b59b6']
    
    # First pass: collect all graft data to determine slice alignment
    all_group_data = {}
    graft_slice_counts = []  # For AP and APTR
    
    for group in ['AP', 'APTR']:
        print(f"\n--- Collecting {group} data{shrink_label} ---")
        all_group_data[group] = {}
        
        for (subj, ses), color in zip(all_subjects[group], colors):
            print(f"\nProcessing {subj}...")
            try:
                slice_data = analyze_graft_fod(subj, ses, hemisphere, config, use_shrink=use_shrink)
                all_group_data[group][subj] = {
                    'data': slice_data,
                    'color': color,
                    'session': ses
                }
                graft_slice_counts.append(len(slice_data))
                print(f"  Found {len(slice_data)} slices with data")
            except Exception as e:
                print(f"  ERROR: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    # Determine target number of slices (minimum graft size)
    if not graft_slice_counts:
        print("No graft data found!")
        return
    
    target_n_slices = min(graft_slice_counts)
    print(f"\nTarget slice count (smallest graft): {target_n_slices}")
    
    # Re-index graft data to 0-based (aligned from anterior)
    for group in ['AP', 'APTR']:
        for subj, subj_info in all_group_data[group].items():
            old_data = subj_info['data']
            sorted_slices = sorted(old_data.keys())[:target_n_slices]
            new_data = {i: old_data[s] for i, s in enumerate(sorted_slices)}
            subj_info['data'] = new_data
    
    # Now collect CTRL data using isocortex ROI (same anatomical region as graft)
    print(f"\n--- Collecting CTRL data (isocortex ROI) ---")
    all_group_data['CTRL'] = {}
    ctrl_slice_counts = []
    
    for (subj, ses), color in zip(all_subjects['CTRL'], colors):
        print(f"\nProcessing {subj}...")
        try:
            slice_data = analyze_isocortex_fod(subj, ses, hemisphere, config)
            all_group_data['CTRL'][subj] = {
                'data': slice_data,
                'color': color,
                'session': ses
            }
            ctrl_slice_counts.append(len(slice_data))
            print(f"  Found {len(slice_data)} slices with data")
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Re-index CTRL data to match target slice count (aligned from anterior)
    for subj, subj_info in all_group_data['CTRL'].items():
        old_data = subj_info['data']
        sorted_slices = sorted(old_data.keys())[:target_n_slices]
        new_data = {i: old_data[s] for i, s in enumerate(sorted_slices)}
        subj_info['data'] = new_data
    
    # Find global max WM content for circle normalization
    global_max_wm = 0
    global_max_amplitude = 0  # For radar scaling
    for group_data in all_group_data.values():
        for subj_info in group_data.values():
            for sd in subj_info['data'].values():
                if 'wm_content' in sd:
                    global_max_wm = max(global_max_wm, sd['wm_content'])
                # Track max amplitude for radar scaling
                max_comp = max(sd['ap_sum'], sd['dv_sum'], sd['lr_sum'])
                global_max_amplitude = max(global_max_amplitude, max_comp)
    
    print(f"\nGlobal max WM content: {global_max_wm:.1f}")
    print(f"Global max FOD amplitude component: {global_max_amplitude:.1f}")
    
    # Calculate average values per slice per group
    avg_wm_per_group = {}
    avg_fod_per_group = {}  # Store raw FOD sums for averaging
    
    for group in ['CTRL', 'AP', 'APTR']:
        group_data = all_group_data[group]
        avg_wm = []
        avg_ap, avg_dv, avg_lr = [], [], []
        
        for slice_idx in range(target_n_slices):
            wm_sum, ap_sum, dv_sum, lr_sum = 0, 0, 0, 0
            count = 0
            for subj_info in group_data.values():
                if slice_idx in subj_info['data']:
                    sd = subj_info['data'][slice_idx]
                    if 'wm_content' in sd:
                        wm_sum += sd['wm_content']
                    ap_sum += sd['ap_sum']
                    dv_sum += sd['dv_sum']
                    lr_sum += sd['lr_sum']
                    count += 1
            
            avg_wm.append(wm_sum / count if count > 0 else 0)
            avg_ap.append(ap_sum / count if count > 0 else 0)
            avg_dv.append(dv_sum / count if count > 0 else 0)
            avg_lr.append(lr_sum / count if count > 0 else 0)
        
        avg_wm_per_group[group] = avg_wm
        avg_fod_per_group[group] = {'ap': avg_ap, 'dv': avg_dv, 'lr': avg_lr}
    
    # Create figure with extra column for info/legend
    n_cols = target_n_slices + 3  # +1 for per-subject avg, +1 for group avg, +1 for info panel
    n_rows = 3
    
    fig = plt.figure(figsize=(2.2 * n_cols, 3.5 * n_rows))
    
    # Direction categories for radar
    categories = ['DV', 'AP', 'LR']
    angles = [n / 3.0 * 2 * np.pi for n in range(3)]
    angles += angles[:1]
    
    row_groups = ['CTRL', 'AP', 'APTR']
    group_labels = ['CTRL\n(isocortex)', 'AP\n(graft, no Tx)', 'APTR\n(graft + Tx)']
    group_colors_bg = ['#e8f5e9', '#ffebee', '#e3f2fd']  # Light backgrounds
    
    for row_idx, (group, label, bg_color) in enumerate(zip(row_groups, group_labels, group_colors_bg)):
        group_data = all_group_data[group]
        
        # Plot individual slices
        for col_idx in range(target_n_slices):
            ax_idx = row_idx * n_cols + col_idx + 1
            ax = fig.add_subplot(n_rows, n_cols, ax_idx, projection='polar')
            ax.set_facecolor(bg_color)
            
            # Plot radar for each subject using RAW SUMS (not fractions)
            for subj, subj_info in group_data.items():
                data = subj_info['data']
                color = subj_info['color']
                
                if col_idx in data:
                    sd = data[col_idx]
                    # Use raw amplitude sums, normalized to global max
                    values = [sd['dv_sum'] / global_max_amplitude,
                             sd['ap_sum'] / global_max_amplitude,
                             sd['lr_sum'] / global_max_amplitude]
                    values += values[:1]
                    
                    ax.plot(angles, values, 'o-', linewidth=1.5, color=color,
                           markersize=3, alpha=0.7)
                    ax.fill(angles, values, alpha=0.1, color=color)
            
            # Center circle for average WM content (scaled to global max)
            avg_wm = avg_wm_per_group[group][col_idx]
            # Circle area proportional to WM, so radius ~ sqrt(WM/max)
            circle_radius = 0.25 * np.sqrt(avg_wm / global_max_wm) if global_max_wm > 0 else 0
            
            # Draw transparent circle at center (outline only)
            if circle_radius > 0.005:
                circle_theta = np.linspace(0, 2*np.pi, 50)
                circle_r = np.ones_like(circle_theta) * circle_radius
                ax.plot(circle_theta, circle_r, color='black', 
                       linewidth=1.5, alpha=0.9, zorder=10)
            
            # Styling
            ax.set_theta_zero_location("S")
            ax.set_theta_direction(1)
            ax.set_ylim(0, 1.0)  # Normalized to global max
            # Add radial ticks with units (as fraction of max)
            ax.set_yticks([0.25, 0.5, 0.75])
            if col_idx == 0:  # Only show labels on first column
                ax.set_yticklabels(['25%', '50%', '75%'], fontsize=5)
            else:
                ax.set_yticklabels(['', '', ''], fontsize=5)
            ax.grid(True, alpha=0.3)
            
            if row_idx == 0:
                ax.set_title(f'Slice {col_idx}', fontsize=9, fontweight='bold', pad=8)
            
            if col_idx == 0:
                ax.set_xticks(angles[:-1])
                ax.set_xticklabels(categories, fontsize=7, fontweight='bold')
            else:
                ax.set_xticks([])
        
        # AVERAGE column (second to last)
        ax_avg = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + target_n_slices + 1, projection='polar')
        ax_avg.set_facecolor('#fffde7')  # Light yellow for average
        
        # Plot average radar for each subject (average across all slices)
        for subj, subj_info in group_data.items():
            data = subj_info['data']
            color = subj_info['color']
            
            # Average across all slices for this subject
            subj_dv = np.mean([data[i]['dv_sum'] for i in range(target_n_slices) if i in data])
            subj_ap = np.mean([data[i]['ap_sum'] for i in range(target_n_slices) if i in data])
            subj_lr = np.mean([data[i]['lr_sum'] for i in range(target_n_slices) if i in data])
            
            values = [subj_dv / global_max_amplitude,
                     subj_ap / global_max_amplitude,
                     subj_lr / global_max_amplitude]
            values += values[:1]
            
            ax_avg.plot(angles, values, 'o-', linewidth=1.5, color=color,
                       markersize=3, alpha=0.7)
            ax_avg.fill(angles, values, alpha=0.1, color=color)
        
        # Average WM circle
        avg_wm_all = np.mean(avg_wm_per_group[group])
        circle_radius_avg = 0.25 * np.sqrt(avg_wm_all / global_max_wm) if global_max_wm > 0 else 0
        
        if circle_radius_avg > 0.005:
            circle_theta = np.linspace(0, 2*np.pi, 50)
            circle_r = np.ones_like(circle_theta) * circle_radius_avg
            ax_avg.plot(circle_theta, circle_r, color='black', 
                       linewidth=1.5, alpha=0.9, zorder=10)
        
        # Styling for average
        ax_avg.set_theta_zero_location("S")
        ax_avg.set_theta_direction(1)
        ax_avg.set_ylim(0, 1.0)
        ax_avg.set_yticks([0.25, 0.5, 0.75])
        ax_avg.set_yticklabels(['', '', ''], fontsize=5)
        ax_avg.grid(True, alpha=0.3)
        ax_avg.set_xticks([])
        
        if row_idx == 0:
            ax_avg.set_title('AVG/subj', fontsize=8, fontweight='bold', pad=8, color='darkred')
        
        # GROUP AVERAGE column (average across all mice - single triangle)
        ax_grp = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + target_n_slices + 2, projection='polar')
        ax_grp.set_facecolor('#ffecb3')  # Darker yellow/orange for group average
        
        # Calculate group average (average across all subjects and slices)
        grp_dv = np.mean(avg_fod_per_group[group]['dv'])
        grp_ap = np.mean(avg_fod_per_group[group]['ap'])
        grp_lr = np.mean(avg_fod_per_group[group]['lr'])
        
        grp_values = [grp_dv / global_max_amplitude,
                     grp_ap / global_max_amplitude,
                     grp_lr / global_max_amplitude]
        grp_values += grp_values[:1]
        
        # Single thick line for group average
        ax_grp.plot(angles, grp_values, 'o-', linewidth=3, color='black',
                   markersize=6, alpha=0.9)
        ax_grp.fill(angles, grp_values, alpha=0.3, color=bg_color)
        
        # Group average WM circle
        grp_wm = np.mean(avg_wm_per_group[group])
        circle_radius_grp = 0.25 * np.sqrt(grp_wm / global_max_wm) if global_max_wm > 0 else 0
        
        if circle_radius_grp > 0.005:
            circle_theta = np.linspace(0, 2*np.pi, 50)
            circle_r = np.ones_like(circle_theta) * circle_radius_grp
            ax_grp.plot(circle_theta, circle_r, color='darkred', 
                       linewidth=2, alpha=0.9, zorder=10)
        
        # Styling for group average
        ax_grp.set_theta_zero_location("S")
        ax_grp.set_theta_direction(1)
        ax_grp.set_ylim(0, 1.0)
        ax_grp.set_yticks([0.25, 0.5, 0.75])
        ax_grp.set_yticklabels(['', '', ''], fontsize=5)
        ax_grp.grid(True, alpha=0.3)
        ax_grp.set_xticks([])
        
        if row_idx == 0:
            ax_grp.set_title('GROUP\nAVG', fontsize=8, fontweight='bold', pad=8, color='darkred')
        
        # Info panel (last column)
        ax_info = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + n_cols)
        ax_info.set_facecolor(bg_color)
        ax_info.axis('off')
        
        # Calculate group statistics
        total_wm = np.mean(avg_wm_per_group[group])
        total_ap = np.mean(avg_fod_per_group[group]['ap'])
        total_dv = np.mean(avg_fod_per_group[group]['dv'])
        total_lr = np.mean(avg_fod_per_group[group]['lr'])
        total_fod = total_ap + total_dv + total_lr
        
        info_text = f"{label}\n\n"
        info_text += f"WM (c₀): {total_wm:.0f}\n({100*total_wm/global_max_wm:.0f}% of max)\n\n"
        info_text += f"FOD amp: {total_fod:.0f}\n(a.u., Σ|dir|×amp)\n\n"
        info_text += f"Direction %:\n"
        if total_fod > 0:
            info_text += f"  DV: {100*total_dv/total_fod:.0f}%\n"
            info_text += f"  AP: {100*total_ap/total_fod:.0f}%\n"
            info_text += f"  LR: {100*total_lr/total_fod:.0f}%"
        
        ax_info.text(0.5, 0.5, info_text, transform=ax_info.transAxes,
                    fontsize=9, ha='center', va='center',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add scale reference at bottom
    # Create a small legend area
    ax_legend = fig.add_axes([0.1, 0.02, 0.8, 0.04])
    ax_legend.axis('off')
    
    scale_text = (f"SCALE: Radar = FOD amplitude (max={global_max_amplitude:.0f} a.u.). "
                  f"Circle radius ∝ √(WM/max), max WM c₀={global_max_wm:.0f}. "
                  f"Units: a.u. = Σ(|direction component| × FOD amplitude at direction).")
    ax_legend.text(0.5, 0.5, scale_text, transform=ax_legend.transAxes,
                  fontsize=8, ha='center', va='center',
                  bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.9))
    
    # Color legend for mice
    legend_elements = [plt.Line2D([0], [0], color=c, linewidth=2, marker='o',
                                  label=f'Mouse {i+1}')
                      for i, c in enumerate(colors)]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=8,
              bbox_to_anchor=(0.99, 0.98), title='Subjects')
    
    shrink_title = " (Shrunk ROI -1 voxel)" if use_shrink else ""
    fig.suptitle(f'FOD Orientation & WM Content - {hemisphere.upper()} Hemisphere{shrink_title}\n'
                 f'Radar: FOD amplitude (DV/AP/LR), Circle: WM content (c₀)',
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.06, 1, 0.95])
    
    suffix = '_shrink' if use_shrink else ''
    output_path = os.path.join(output_dir, f'fod_3group_comparison_{hemisphere}{suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_path}")
    plt.close()
    
    return output_path


def analyze_probe_fod(subject, session, config, roi_path):
    """
    Analyze FOD directions within a single-slice probe ROI.
    
    Parameters:
    -----------
    subject : str - subject ID
    session : str - session ID  
    config : dict - configuration
    roi_path : str - path to the probe ROI (already in FOD space)
    
    Returns:
    --------
    dict with {n_voxels, ap_sum, dv_sum, lr_sum, ap_fraction, dv_fraction, lr_fraction, total_amplitude, wm_content}
    """
    deriv = config['paths']['derivatives_flip1']
    base = os.path.join(deriv, subject, session)
    
    # Load FOD
    fod_path = os.path.join(base, 'dwi', 'wm_fod_norm.mif')
    fod_data = load_mif(fod_path)
    print(f"  FOD shape: {fod_data.shape}")
    
    # Load probe ROI
    mask_data = load_mif(roi_path) > 0.5
    n_voxels = mask_data.sum()
    print(f"  Probe ROI voxels: {n_voxels}")
    
    if n_voxels == 0:
        print(f"  WARNING: Probe ROI is empty!")
        return None
    
    # Precompute SH basis for sampling
    n_sample = 100
    sample_dirs = get_sphere_vertices(n_sample)
    sh_basis = compute_sh_basis(sample_dirs, lmax=8)
    
    # Find voxel coordinates
    coords = np.where(mask_data)
    
    # Accumulate weighted sums
    total_ap = 0.0
    total_dv = 0.0
    total_lr = 0.0
    total_wm = 0.0
    n_valid_voxels = 0
    
    for x, y, z in zip(coords[0], coords[1], coords[2]):
        sh_coeffs = fod_data[x, y, z, :]
        
        if np.max(np.abs(sh_coeffs)) < 0.01:
            continue
        
        total_wm += sh_coeffs[0]
        
        amplitudes = sh_basis @ sh_coeffs
        amplitudes = np.maximum(amplitudes, 0)
        
        for d, amp in zip(sample_dirs, amplitudes):
            if amp > 0.01:
                total_ap += np.abs(d[0]) * amp
                total_dv += np.abs(d[1]) * amp
                total_lr += np.abs(d[2]) * amp
        
        n_valid_voxels += 1
    
    if n_valid_voxels == 0:
        print(f"  WARNING: No valid voxels with FOD signal!")
        return None
    
    total_all = total_ap + total_dv + total_lr
    
    result = {
        'n_voxels': n_valid_voxels,
        'ap_sum': total_ap,
        'dv_sum': total_dv,
        'lr_sum': total_lr,
        'ap_fraction': total_ap / total_all if total_all > 0 else 0,
        'dv_fraction': total_dv / total_all if total_all > 0 else 0,
        'lr_fraction': total_lr / total_all if total_all > 0 else 0,
        'total_amplitude': total_all,
        'wm_content': total_wm
    }
    
    print(f"  Results: AP={result['ap_fraction']:.2f}, DV={result['dv_fraction']:.2f}, LR={result['lr_fraction']:.2f}")
    
    return result


def analyze_probe_comparison(config, output_dir):
    """
    Generate single radar plots per mouse using single-slice probe ROIs.
    
    Probe ROIs are manually drawn coronal slices at comparable anatomical locations
    for CTRL, AP, and APTR groups.
    """
    deriv = config['paths']['derivatives_flip1']
    
    # Probe ROI definitions (subject, session, roi_path)
    probe_subjects = {
        'CTRL': {
            'subject': 'sub-129-921-R1F',
            'session': 'ses-260115121231',
            'roi_path': os.path.join(deriv, 'sub-129-921-R1F', 'ses-260115121231', 
                                      'atlas', 'rois', 'fod_space', 
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        },
        'AP': {
            'subject': 'sub-769-C-O5-R1M',
            'session': 'ses-240911115733',
            'roi_path': os.path.join(deriv, 'sub-769-C-O5-R1M', 'ses-240911115733',
                                      'manual_roi', 'fod_space',
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        },
        'APTR': {
            'subject': 'sub-AP72-862-901-NmF',
            'session': 'ses-251210093342',
            'roi_path': os.path.join(deriv, 'sub-AP72-862-901-NmF', 'ses-251210093342',
                                      'manual_roi', 'fod_space',
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        }
    }
    
    # Analyze each probe
    results = {}
    for group, info in probe_subjects.items():
        print(f"\n--- Analyzing {group}: {info['subject']} ---")
        
        if not os.path.exists(info['roi_path']):
            print(f"  ERROR: Probe ROI not found: {info['roi_path']}")
            continue
        
        result = analyze_probe_fod(info['subject'], info['session'], config, info['roi_path'])
        
        if result is not None:
            results[group] = {
                'data': result,
                'subject': info['subject'],
                'session': info['session']
            }
    
    if not results:
        print("ERROR: No probe data collected!")
        return
    
    # Generate comparison figure
    plot_probe_radar(results, config, output_dir)


def plot_probe_radar(results, config, output_dir):
    """
    Generate a figure with one radar plot per group (CTRL, AP, APTR).
    Each radar shows FOD directional content within the single-slice probe ROI.
    Uses RAW FOD amplitude sums (not fractions) normalized to global max, like the other comparison plots.
    Center circle shows WM content (c[0] coefficient sum).
    Includes probe ROI location images below each radar.
    """
    from matplotlib.image import imread
    
    n_groups = len(results)
    
    # Probe ROI location images
    figures_dir = os.path.join(output_dir, 'figures')
    roi_images = {
        'CTRL': os.path.join(figures_dir, 'CTRL_with_1slice_ROI.png'),
        'AP': os.path.join(figures_dir, 'AP_with_1slice_ROI.png'),
        'APTR': os.path.join(figures_dir, 'APTR_with_1slice_ROI.png'),
    }
    
    # Create figure with 2 rows: radar plots on top, ROI images below
    fig = plt.figure(figsize=(5 * (n_groups + 1) + 4, 12))
    
    # Radar categories
    categories = ['Dor-Ven', 'Ant-Post', 'Left-Right']
    angles = [n / 3.0 * 2 * np.pi for n in range(3)]
    angles += angles[:1]  # Close polygon
    
    # Group colors
    group_colors = {
        'CTRL': '#27ae60',   # Green
        'AP': '#e74c3c',     # Red  
        'APTR': '#3498db'    # Blue
    }
    
    group_order = ['CTRL', 'AP', 'APTR']
    
    # Find global max for consistent scaling (using RAW SUMS, not fractions)
    max_amplitude = 0
    max_wm = 0
    for info in results.values():
        data = info['data']
        max_amplitude = max(max_amplitude, data['ap_sum'], data['dv_sum'], data['lr_sum'])
        max_wm = max(max_wm, data['wm_content'])
    
    # Top row: radar plots
    n_cols = n_groups + 2  # 3 groups + summary + explanation
    
    # Plot individual radar per group
    for idx, group in enumerate(group_order):
        if group not in results:
            continue
        
        ax = fig.add_subplot(2, n_cols, idx + 1, projection='polar')
        
        info = results[group]
        data = info['data']
        color = group_colors[group]
        
        # Use RAW amplitude sums normalized to global max (like the other comparison plots)
        values = [data['dv_sum'] / max_amplitude,
                  data['ap_sum'] / max_amplitude,
                  data['lr_sum'] / max_amplitude]
        values += values[:1]  # Close polygon
        
        ax.plot(angles, values, 'o-', linewidth=3, color=color, markersize=10)
        ax.fill(angles, values, alpha=0.3, color=color)
        
        # Draw WM content circle at center
        # Circle radius proportional to sqrt(WM/max) for area scaling
        wm_content = data['wm_content']
        circle_radius = 0.25 * np.sqrt(wm_content / max_wm) if max_wm > 0 else 0
        
        if circle_radius > 0.005:
            circle_theta = np.linspace(0, 2*np.pi, 50)
            circle_r = np.ones_like(circle_theta) * circle_radius
            ax.fill(circle_theta, circle_r, color='gray', alpha=0.4, zorder=5)
            ax.plot(circle_theta, circle_r, color='black', linewidth=1.5, alpha=0.9, zorder=10)
        
        # DV at South (bottom)
        ax.set_theta_zero_location("S")
        ax.set_theta_direction(1)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=10, fontweight='bold')
        ax.set_ylim(0, 1.0)  # Normalized to global max
        ax.set_yticks([0.25, 0.5, 0.75])
        ax.set_yticklabels(['25%', '50%', '75%'], fontsize=7)
        ax.grid(True, alpha=0.3)
        
        # Title with subject info and WM content
        subj_short = info['subject'].replace('sub-', '')[:15]
        wm_pct = 100 * wm_content / max_wm if max_wm > 0 else 0
        ax.set_title(f'{group}\n{subj_short}\n({data["n_voxels"]} vox)', 
                     fontsize=11, fontweight='bold', pad=15)
    
    # Summary plot: all groups overlaid (top row, column 4)
    ax_summary = fig.add_subplot(2, n_cols, n_groups + 1, projection='polar')
    
    for group in group_order:
        if group not in results:
            continue
        
        data = results[group]['data']
        color = group_colors[group]
        
        # Use RAW sums normalized to global max
        values = [data['dv_sum'] / max_amplitude,
                  data['ap_sum'] / max_amplitude,
                  data['lr_sum'] / max_amplitude]
        values += values[:1]
        
        ax_summary.plot(angles, values, 'o-', linewidth=3, color=color, 
                        markersize=10, label=group)
        ax_summary.fill(angles, values, alpha=0.15, color=color)
    
    ax_summary.set_theta_zero_location("S")
    ax_summary.set_theta_direction(1)
    ax_summary.set_xticks(angles[:-1])
    ax_summary.set_xticklabels(categories, fontsize=10, fontweight='bold')
    ax_summary.set_ylim(0, 1.0)
    ax_summary.set_yticks([0.25, 0.5, 0.75])
    ax_summary.set_yticklabels(['25%', '50%', '75%'], fontsize=7)
    ax_summary.grid(True, alpha=0.3)
    ax_summary.set_title('All Groups\nOverlaid', fontsize=11, fontweight='bold', pad=15)
    ax_summary.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0), fontsize=9)
    
    # Add explanation panel (top row, rightmost column)
    ax_text = fig.add_subplot(2, n_cols, n_cols)
    ax_text.axis('off')
    
    # Build explanation text
    explanation = """KEY FINDING:
━━━━━━━━━━━━━━━━━━━━━━━━
APTR shows highest WM 
AND strong DV dominance,
supporting that transplanted 
fibers regenerate in the 
DV direction within graft.

INTERPRETATION:
━━━━━━━━━━━━━━━━━━━━━━━━
• Radar size = FOD amplitude
  (water diffusion anisotropy)
  → larger = more organized fibers

• Center circle = WM content
  (c₀ coefficient from FOD)
  → proxy for myelination
  → larger = more white matter

WM CONTENT (relative):
━━━━━━━━━━━━━━━━━━━━━━━━"""
    
    # Add WM values to explanation
    for group in group_order:
        if group in results:
            wm_pct = 100 * results[group]['data']['wm_content'] / max_wm
            dv_pct = 100 * results[group]['data']['dv_fraction']
            explanation += f"\n  {group:>5}: {wm_pct:5.0f}%  (DV={dv_pct:.0f}%)"
    
    ax_text.text(0.05, 0.95, explanation, transform=ax_text.transAxes,
                 fontsize=9, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', 
                          edgecolor='orange', alpha=0.9))
    
    # Bottom row: Probe ROI location images
    for idx, group in enumerate(group_order):
        if group not in results:
            continue
        
        ax_img = fig.add_subplot(2, n_cols, n_cols + idx + 1)
        
        # Load and display the ROI image
        img_path = roi_images.get(group)
        if img_path and os.path.exists(img_path):
            img = imread(img_path)
            ax_img.imshow(img)
            ax_img.set_title(f'{group} - Probe ROI Location', fontsize=10, fontweight='bold')
        else:
            ax_img.text(0.5, 0.5, f'{group}\nImage not found', 
                       ha='center', va='center', fontsize=12)
            ax_img.set_title(f'{group} - Probe ROI', fontsize=10)
        
        ax_img.axis('off')
    
    # Bottom row: merged explanation for last two columns
    ax_bottom_text = fig.add_subplot(2, n_cols, n_cols + n_groups + 1)
    ax_bottom_text.axis('off')
    
    bottom_explanation = """
ROI IMAGES ABOVE:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Yellow/green overlay shows the 
single coronal slice probe ROI 
used for FOD analysis.

These manually-drawn ROIs sample
comparable anatomical locations
across subjects:
• CTRL: native isocortex
• AP: graft without transplant
• APTR: graft with transplant
"""
    
    ax_bottom_text.text(0.05, 0.95, bottom_explanation, transform=ax_bottom_text.transAxes,
                 fontsize=9, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='lightcyan', 
                          edgecolor='steelblue', alpha=0.9))
    
    # Merge last column of bottom row with explanation
    ax_bottom_last = fig.add_subplot(2, n_cols, 2 * n_cols)
    ax_bottom_last.axis('off')
    
    fig.suptitle('FOD Directionality Analysis - Single-Slice Probe ROIs\n'
                 '(Radar: directional FOD amplitude, normalized to max across groups)', 
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    output_path = os.path.join(output_dir, 'fod_probe_comparison.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {output_path}")
    plt.close()
    
    # Also print summary table
    print("\n" + "=" * 70)
    print("PROBE FOD ANALYSIS SUMMARY")
    print("=" * 70)
    print(f"{'Group':>8} {'Subject':>20} {'Voxels':>8} {'AP%':>8} {'DV%':>8} {'LR%':>8} {'WM':>10}")
    print("-" * 70)
    
    for group in group_order:
        if group not in results:
            continue
        info = results[group]
        data = info['data']
        subj = info['subject'].replace('sub-', '')[:18]
        wm_pct = 100 * data['wm_content'] / max_wm if max_wm > 0 else 0
        print(f"{group:>8} {subj:>20} {data['n_voxels']:>8} "
              f"{data['ap_fraction']*100:>7.1f}% {data['dv_fraction']*100:>7.1f}% {data['lr_fraction']*100:>7.1f}% {wm_pct:>9.0f}%")
    
    print("=" * 70)
    print("WM = White matter content (c₀ coefficient sum), normalized to max")
    print("     → proxy for myelination/axonal density")
    print("Circle size in radar plots is proportional to WM content")
    
    return output_path


def main():
    """Main analysis function."""
    parser = argparse.ArgumentParser(description='Analyze FOD directionality in graft ROI')
    parser.add_argument('--group', type=str, default='APTR', 
                        choices=['APTR', 'AP'],
                        help='Subject group to analyze (default: APTR)')
    parser.add_argument('--polar', action='store_true',
                        help='Also generate polar plots (2D projection)')
    parser.add_argument('--compare', action='store_true',
                        help='Generate 3-group comparison (CTRL vs AP vs APTR)')
    parser.add_argument('--shrink', action='store_true',
                        help='Use shrunk graft ROIs (eroded by 1 voxel) for AP/APTR')
    parser.add_argument('--create-shrink-rois', action='store_true',
                        help='Create shrunk (eroded) graft ROIs for all AP/APTR subjects')
    parser.add_argument('--probe', action='store_true',
                        help='Generate single radar plot per mouse using single-slice probe ROIs')
    args = parser.parse_args()
    
    # Load config (use flip1 pipeline)
    config = load_config(pipeline='flip1')
    
    output_dir = os.path.join(config['paths']['repo_dir'], 'reports')
    os.makedirs(output_dir, exist_ok=True)
    
    # If create-shrink-rois mode, create the eroded ROIs
    if args.create_shrink_rois:
        print("\n" + "=" * 60)
        print("Creating shrunk (eroded) graft ROIs")
        print("=" * 60)
        create_shrunk_graft_rois(config)
        print("\nShrunk ROI creation complete!")
        return
    
    # If compare mode, run 3-group comparison
    if args.compare:
        for hemisphere in ['left', 'right']:
            analyze_all_groups_comparison(hemisphere, config, output_dir, use_shrink=args.shrink)
        shrink_msg = " (with shrunk ROIs)" if args.shrink else ""
        print("\n" + "=" * 60)
        print(f"3-group comparison complete!{shrink_msg}")
        return
    
    # If probe mode, generate single radar plots per mouse using probe ROIs
    if args.probe:
        print("\n" + "=" * 60)
        print("PROBE MODE: Single radar plot per mouse using single-slice probe ROIs")
        print("=" * 60)
        analyze_probe_comparison(config, output_dir)
        return
    
    # Otherwise run single-group analysis
    # Subject definitions
    subjects = {
        'APTR': [
            ('sub-AP72-862-901-NmF', 'ses-251210093342'),
            ('sub-AP72-862-901-R1M', 'ses-251210100043'),
            ('sub-Apallial17NmF', 'ses-240816091633'),
            ('sub-Apallial20R1F', 'ses-240909094232'),
        ],
        'AP': [
            ('sub-286-907-OL-R1F', 'ses-251210115606'),
            ('sub-596-C-OO-NmM', 'ses-250115101115'),
            ('sub-AP-72-885-901-NmM', 'ses-251210105611'),
            ('sub-Apallial19R1M', 'ses-240909091132'),
        ],
    }
    
    # Colors for mice (distinguishable, colorblind-friendly)
    colors = ['#e74c3c', '#3498db', '#27ae60', '#9b59b6']  # Red, Blue, Green, Purple
    
    group_subjects = subjects[args.group]
    
    for hemisphere in ['left', 'right']:
        print(f"\n{'='*60}")
        print(f"Analyzing {args.group} - {hemisphere.upper()} hemisphere")
        print(f"{'='*60}")
        
        all_subject_data = {}
        all_slices = set()
        max_amplitude = 0
        
        # Collect data from all subjects
        for (subj, ses), color in zip(group_subjects, colors):
            print(f"\nProcessing {subj}...")
            try:
                slice_data = analyze_graft_fod(subj, ses, hemisphere, config)
                all_subject_data[subj] = {
                    'data': slice_data,
                    'color': color,
                    'session': ses
                }
                all_slices.update(slice_data.keys())
                
                for sd in slice_data.values():
                    total_amp = sd['total_amplitude']
                    if total_amp > max_amplitude:
                        max_amplitude = total_amp
                
                print(f"  Found {len(slice_data)} slices with data")
            except Exception as e:
                print(f"  ERROR: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        if not all_subject_data:
            print(f"No data found for {hemisphere} hemisphere!")
            continue
        
        # Generate 3D orientation plot
        plot_3d_orientations(all_subject_data, hemisphere, args.group, config, output_dir)
        
        # Generate composite summary
        plot_composite_summary(all_subject_data, hemisphere, args.group, config, output_dir)
        
        # Generate polar plots if requested
        if args.polar:
            generate_polar_plots(all_subject_data, hemisphere, args.group, config, output_dir)
        
        # Save detailed text summary
        save_text_summary(all_subject_data, hemisphere, args.group, output_dir)
    
    print("\n" + "=" * 60)
    print("Analysis complete!")


def save_text_summary(subject_data, hemisphere, group, output_dir):
    """Save detailed text summary with 3D orientation info."""
    summary_path = os.path.join(output_dir, f'fod_details_{group.lower()}_{hemisphere}.txt')
    
    with open(summary_path, 'w') as f:
        f.write(f"FOD 3D Orientation Summary - {group} {hemisphere.upper()} Graft\n")
        f.write("=" * 70 + "\n\n")
        f.write("Axes: X=Anterior-Posterior, Y=Dorsal-Ventral, Z=Left-Right\n")
        f.write("Orientation components are absolute values (antipodal symmetry)\n\n")
        
        for subj, subj_info in subject_data.items():
            f.write(f"\n{subj}:\n")
            f.write("-" * 50 + "\n")
            data = subj_info['data']
            
            f.write(f"{'Slice':>6} {'Voxels':>6} {'TotAmp':>9} {'AP%':>6} {'DV%':>6} {'LR%':>6}\n")
            f.write("-" * 50 + "\n")
            
            for slice_idx in sorted(data.keys()):
                sd = data[slice_idx]
                
                f.write(f"{slice_idx:>6} {sd['n_voxels']:>6} {sd['total_amplitude']:>9.1f} "
                       f"{sd['ap_fraction']*100:>5.1f}% {sd['dv_fraction']*100:>5.1f}% {sd['lr_fraction']*100:>5.1f}%\n")
            
            # Overall summary for this subject
            total_ap = sum(sd['ap_sum'] for sd in data.values())
            total_dv = sum(sd['dv_sum'] for sd in data.values())
            total_lr = sum(sd['lr_sum'] for sd in data.values())
            total_all = total_ap + total_dv + total_lr
            
            if total_all > 0:
                ap_pct = 100 * total_ap / total_all
                dv_pct = 100 * total_dv / total_all
                lr_pct = 100 * total_lr / total_all
            else:
                ap_pct = dv_pct = lr_pct = 0
            
            f.write(f"\nOverall weighted fractions: "
                   f"A-P={ap_pct:.1f}%, "
                   f"D-V={dv_pct:.1f}%, "
                   f"L-R={lr_pct:.1f}%\n")
            
            # Dominant direction
            components = {'A-P': ap_pct, 'D-V': dv_pct, 'L-R': lr_pct}
            dominant = max(components, key=components.get)
            f.write(f"Dominant orientation: {dominant} ({components[dominant]:.1f}%)\n")
    
    print(f"Saved: {summary_path}")


def generate_polar_plots(subject_data, hemisphere, group, config, output_dir):
    """Generate legacy polar plots (2D projection to D-V vs L-R plane).
    
    NOTE: This function requires 'directions' and 'magnitudes' fields which have been
    removed from the new data structure. This is now deprecated in favor of the 
    radar plots in plot_3d_orientations().
    """
    print("  WARNING: Polar plot generation is deprecated (data structure changed)")
    print("  Use the radar plots in the main output instead")
    return


if __name__ == '__main__':
    main()
