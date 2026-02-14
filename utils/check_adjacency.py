#!/usr/bin/env python
"""
Check spatial adjacency of ROI segments.

This script analyzes whether segments appear to overlap visually in MRview
due to being adjacent (touching) vs actually overlapping.
"""
import nibabel as nib
import numpy as np

base = 'data/derivatives/mrtrix_preproc_flip1/sub-129-921-R1F/ses-260115121231/dwi/isocortex_segments_v2'
segments = ['anterior', 'anterior_central', 'central', 'posterior_central', 'posterior']

# Load all left segments and check their spatial extent
print('Spatial extent (axis 0 = A-P) for each segment:')
for seg in segments:
    img = nib.load(f'{base}/isocortex_left_{seg}.nii.gz')
    data = img.get_fdata() > 0
    x_indices = np.where(np.any(data, axis=(1,2)))[0]
    if len(x_indices) > 0:
        print(f'  {seg:20s}: X range [{x_indices.min():3d} - {x_indices.max():3d}] (span: {x_indices.max()-x_indices.min()+1})')

print()
print('Checking if segments are ADJACENT (touching) vs OVERLAPPING:')
rois = {}
for seg in segments:
    img = nib.load(f'{base}/isocortex_left_{seg}.nii.gz')
    rois[seg] = img.get_fdata() > 0

# Simple dilation without scipy
def dilate_3d(arr):
    """Simple 3D dilation by 1 voxel (6-connectivity)."""
    result = arr.copy()
    result[1:, :, :] |= arr[:-1, :, :]
    result[:-1, :, :] |= arr[1:, :, :]
    result[:, 1:, :] |= arr[:, :-1, :]
    result[:, :-1, :] |= arr[:, 1:, :]
    result[:, :, 1:] |= arr[:, :, :-1]
    result[:, :, :-1] |= arr[:, :, 1:]
    return result

for i, seg1 in enumerate(segments[:-1]):
    seg2 = segments[i+1]
    dilated = dilate_3d(rois[seg1])
    adjacent = np.sum(dilated & rois[seg2])
    actual_overlap = np.sum(rois[seg1] & rois[seg2])
    print(f'  {seg1} <-> {seg2}: adjacent={adjacent} voxels, overlap={actual_overlap}')
