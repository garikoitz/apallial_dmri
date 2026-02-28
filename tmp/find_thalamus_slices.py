#!/usr/bin/env python3
"""Find first thalamus slice (axis 0, anterior->posterior) for each subject."""
import nibabel as nib
import numpy as np
import subprocess, tempfile, os

DERIV = 'data/derivatives/mrtrix_preproc_flip1'

subjects = {
    'CTRL': f'{DERIV}/sub-138-1010-L1F/ses-260115101441/atlas/rois/fod_space/roi_thalamus_left.nii.gz',
    'AP': f'{DERIV}/sub-769-C-O5-R1M/ses-240911115733/manual_roi/fod_space/roi_Thalamus_left_flip_fod.mif',
    'APTR': f'{DERIV}/sub-AP72-862-901-NmF/ses-251210093342/manual_roi/fod_space/roi_Thalamus_left_flip_fod.mif',
}

for group, path in subjects.items():
    if path.endswith('.nii.gz') or path.endswith('.nii'):
        img = nib.load(path)
        data = img.get_fdata()
    else:
        tmp_file = tempfile.mktemp(suffix='.nii.gz')
        subprocess.run(['mrconvert', path, tmp_file, '-force'], capture_output=True)
        img = nib.load(tmp_file)
        data = img.get_fdata()
        os.remove(tmp_file)
    
    print(f'{group}: shape={data.shape}')
    nonzero = np.where(data > 0)
    if len(nonzero[0]) > 0:
        slices_axis0 = sorted(set(nonzero[0]))
        print(f'  First thalamus slice (axis 0): {slices_axis0[0]}')
        print(f'  Last thalamus slice (axis 0): {slices_axis0[-1]}')
        print(f'  All thalamus slices: {slices_axis0}')
        print(f'  Total voxels: {len(nonzero[0])}')
    else:
        print('  NO nonzero voxels!')
    print()
