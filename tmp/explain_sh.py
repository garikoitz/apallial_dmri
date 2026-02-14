#!/usr/bin/env python3
"""Explain what's in the SH file."""

import numpy as np
import subprocess
import tempfile
import nibabel as nib
import os
from scipy.special import sph_harm

# Load one FOD file
deriv = 'data/derivatives/mrtrix_preproc_flip1'
subj = 'sub-Apallial17NmF'
ses = 'ses-240816091633'
fod_path = f'{deriv}/{subj}/{ses}/dwi/wm_fod_norm.mif'

with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
    tmp_path = tmp.name
subprocess.run(['mrconvert', fod_path, tmp_path, '-force', '-quiet'], check=True)
fod = nib.load(tmp_path).get_fdata()
os.remove(tmp_path)

print('FOD SH File Structure:')
print('=' * 50)
print(f'Shape: {fod.shape}')
print(f'  - Dimensions 0,1,2: spatial (30 x 100 x 120 voxels)')
print(f'  - Dimension 3: {fod.shape[3]} SH coefficients (lmax=8)')
print()

# Pick a voxel in the graft
x, y, z = 15, 50, 60
sh = fod[x, y, z, :]

print(f'Example voxel [{x},{y},{z}] SH coefficients:')
print(f'  Total coefficients: {len(sh)}')
print(f'  l=0 (1 coeff):   coef[0] = {sh[0]:.4f}  <- isotropic')
print(f'  l=2 (5 coeffs):  coef[1:6] range: [{sh[1:6].min():.3f}, {sh[1:6].max():.3f}]')
print(f'  l=4 (9 coeffs):  coef[6:15] range: [{sh[6:15].min():.3f}, {sh[6:15].max():.3f}]')
print(f'  l=6 (13 coeffs): coef[15:28] range: [{sh[15:28].min():.3f}, {sh[15:28].max():.3f}]')
print(f'  l=8 (17 coeffs): coef[28:45] range: [{sh[28:45].min():.3f}, {sh[28:45].max():.3f}]')
print()

def sample_fod_direction(sh, direction):
    """Sample FOD at a Cartesian direction [x,y,z]."""
    x, y, z = direction
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / (r + 1e-10))  # polar
    phi = np.arctan2(y, x)  # azimuth
    
    lmax = 8
    idx = 0
    amplitude = 0
    for l in range(0, lmax + 1, 2):
        for m in range(-l, l + 1):
            if m < 0:
                Y = np.sqrt(2) * (-1)**m * np.imag(sph_harm(abs(m), l, phi, theta))
            elif m == 0:
                Y = np.real(sph_harm(0, l, phi, theta))
            else:
                Y = np.sqrt(2) * (-1)**m * np.real(sph_harm(m, l, phi, theta))
            amplitude += sh[idx] * float(np.real(Y))
            idx += 1
    return max(amplitude, 0)

print('What the SH coefficients encode:')
print('-' * 50)
print('The FOD is a function on the sphere: f(direction) -> amplitude')
print('SH coefficients are like Fourier coefficients for spherical data')
print()
print('To get FOD amplitude in any direction, we evaluate:')
print('  amplitude = sum(coef[i] * Y_l^m(direction))')
print()
print('where Y_l^m are spherical harmonic basis functions.')
print()

print('FOD amplitude at different directions for this voxel:')
print('  Direction         Cartesian [A-P, D-V, L-R]    Amplitude')
print('  ' + '-' * 55)

# In voxel space: axis 0 = A-P, axis 1 = D-V, axis 2 = L-R
directions = [
    ('Anterior (+A-P)',  [1, 0, 0]),
    ('Posterior (-A-P)', [-1, 0, 0]),
    ('Dorsal (+D-V)',    [0, 1, 0]),
    ('Ventral (-D-V)',   [0, -1, 0]),
    ('Left (+L-R)',      [0, 0, 1]),
    ('Right (-L-R)',     [0, 0, -1]),
    ('Diagonal A-P/D-V', [0.707, 0.707, 0]),
    ('Diagonal A-P/L-R', [0.707, 0, 0.707]),
]

for name, d in directions:
    amp = sample_fod_direction(sh, d)
    print(f'  {name:20s} {str(d):25s} {amp:.4f}')

print()
print('KEY POINT:')
print('=' * 50)
print('The amplitude IS the FOD magnitude in that direction.')
print('Higher amplitude = more fibers oriented that way.')
print()
print('For our analysis, we should:')
print('1. Sample FOD at many directions')
print('2. For each direction, the AMPLITUDE is the "weight"')
print('3. Compute weighted average: sum(direction * amplitude) / sum(amplitude)')
