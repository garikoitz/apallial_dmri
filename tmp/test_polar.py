#!/usr/bin/env python3
"""Quick test of polar plot generation with real FOD data."""

import numpy as np
import subprocess
import os
import tempfile
import nibabel as nib
import matplotlib.pyplot as plt

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
    """Compute real spherical harmonic basis for given directions."""
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
            
            basis[:, idx] = np.real(Y)  # Ensure real
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

# Test loading one subject
deriv = 'data/derivatives/mrtrix_preproc_flip1'
subj = 'sub-Apallial17NmF'
ses = 'ses-240816091633'
base = os.path.join(deriv, subj, ses)

# Load FOD and mask
fod_data = load_mif(os.path.join(base, 'dwi', 'wm_fod_norm.mif'))
mask_data = load_mif(os.path.join(base, 'manual_roi', 'fod_space', 'roi_Graft_left_flip_fod.mif')) > 0.5

print(f'FOD shape: {fod_data.shape}')
print(f'Mask voxels: {mask_data.sum()}')

# Find A-P slices
mask_slices = np.any(mask_data, axis=(1, 2))
slice_indices = np.where(mask_slices)[0]
print(f'Graft spans slices: {slice_indices.min()} to {slice_indices.max()} ({len(slice_indices)} slices)')

# Compute SH basis
n_sample = 100
sample_dirs = get_sphere_vertices(n_sample)
sh_basis = compute_sh_basis(sample_dirs, lmax=8)
print(f'SH basis shape: {sh_basis.shape}')

# Get a test slice
slice_idx = slice_indices[len(slice_indices)//2]
slice_mask = mask_data[slice_idx, :, :]
y_coords, z_coords = np.where(slice_mask)
print(f'Slice {slice_idx} has {len(y_coords)} voxels')

# Collect all directions and magnitudes for this slice
directions_all = []
magnitudes_all = []

for y, z in zip(y_coords, z_coords):
    sh_coeffs = fod_data[slice_idx, y, z, :]
    
    if np.max(np.abs(sh_coeffs)) < 0.01:
        continue
    
    # Sample FOD
    amplitudes = sh_basis @ sh_coeffs
    amplitudes = np.maximum(amplitudes, 0)
    
    if amplitudes.max() > 0.05:
        # Weight directions by amplitude^2 to emphasize peaks
        weights = amplitudes ** 2
        weighted_dir = np.sum(sample_dirs * weights[:, np.newaxis], axis=0)
        magnitude = np.linalg.norm(weighted_dir)
        
        if magnitude > 0.01:
            direction = weighted_dir / magnitude
            directions_all.append(direction)
            magnitudes_all.append(magnitude)

directions = np.array(directions_all)
magnitudes = np.array(magnitudes_all)
print(f'Found {len(directions)} valid voxels')
print(f'Direction stats: mean norm = {np.linalg.norm(directions, axis=1).mean():.3f}')
print(f'Magnitude range: {magnitudes.min():.4f} to {magnitudes.max():.4f}')

# Create polar histogram
# Project to 2D: Y (axis 1 = D-V) vs Z (axis 2 = L-R)
angles = np.arctan2(directions[:, 1], directions[:, 2])  # Y vs Z

print(f'Angle range: {np.degrees(angles.min()):.1f} to {np.degrees(angles.max()):.1f} degrees')

# Create figure
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Subplot 1: 2D scatter of directions (Y vs Z)
ax1 = axes[0]
sc = ax1.scatter(directions[:, 2], directions[:, 1], c=magnitudes, 
                 cmap='viridis', alpha=0.5, s=10)
ax1.set_xlabel('Z (L-R)')
ax1.set_ylabel('Y (D-V)')
ax1.set_title(f'Slice {slice_idx}: Direction vectors (Y vs Z)')
ax1.set_aspect('equal')
ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax1.axvline(0, color='gray', linestyle='--', alpha=0.5)
plt.colorbar(sc, ax=ax1, label='Magnitude')

# Subplot 2: Polar histogram
ax2 = fig.add_subplot(1, 3, 2, projection='polar')
n_bins = 36
bars, bin_edges = np.histogram(angles, bins=n_bins, weights=magnitudes, range=(-np.pi, np.pi))
if bars.max() > 0:
    bars_norm = bars / bars.max()
else:
    bars_norm = bars

width = (2 * np.pi) / n_bins
ax2.bar(bin_edges[:-1], bars_norm, width=width, color='blue', alpha=0.6, edgecolor='navy')
ax2.set_title(f'Slice {slice_idx}: Polar histogram\n(weighted by magnitude)')
ax2.set_theta_zero_location('N')

# Subplot 3: Average direction arrow
ax3 = fig.add_subplot(1, 3, 3, projection='polar')
mean_dir = np.mean(directions * magnitudes[:, np.newaxis], axis=0)
mean_mag = np.linalg.norm(mean_dir)
mean_dir_norm = mean_dir / (mean_mag + 1e-10)
mean_angle = np.arctan2(mean_dir_norm[1], mean_dir_norm[2])

ax3.arrow(mean_angle, 0, 0, 0.8, head_width=0.1, head_length=0.1, fc='red', ec='darkred', linewidth=2)
ax3.set_title(f'Slice {slice_idx}: Mean direction\n(angle={np.degrees(mean_angle):.1f}°, mag={mean_mag:.3f})')
ax3.set_theta_zero_location('N')
ax3.set_ylim(0, 1)

plt.tight_layout()
plt.savefig('tmp/test_polar.png', dpi=150, bbox_inches='tight')
print('\nSaved: tmp/test_polar.png')
plt.show()
