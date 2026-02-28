#!/usr/bin/env python3
"""Find first-thalamus FOD slice (axis 0, anterior→posterior) for ALL subjects."""
import nibabel as nib
import numpy as np
import subprocess, tempfile, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions.config_loader import load_config

cfg = load_config('config1.yaml')
DERIV = cfg['paths']['derivatives_flip1']

groups_cfg = cfg.get('groups', {})

for group in ['CTRL', 'AP', 'APTR']:
    print(f"=== {group} ===")
    for entry in groups_cfg.get(group, []):
        sid = entry['id']
        ses = entry['session']
        if 'test' in sid.lower():
            continue
        subj_dir = os.path.join(DERIV, sid, ses)
        
        # Determine thalamus ROI path
        if group == 'CTRL':
            roi_path = os.path.join(subj_dir, 'atlas', 'rois', 'fod_space', 'roi_thalamus_left.nii.gz')
        else:
            roi_path = os.path.join(subj_dir, 'manual_roi', 'fod_space', 'roi_Thalamus_left_flip_fod.mif')
        
        if not os.path.exists(roi_path):
            print(f"  {sid}: ROI NOT FOUND ({roi_path})")
            continue
        
        # Load data
        if roi_path.endswith('.nii.gz') or roi_path.endswith('.nii'):
            img = nib.load(roi_path)
            data = img.get_fdata()
        else:
            tmp_file = tempfile.mktemp(suffix='.nii.gz')
            subprocess.run(['mrconvert', roi_path, tmp_file, '-force'],
                          capture_output=True)
            img = nib.load(tmp_file)
            data = img.get_fdata()
            os.remove(tmp_file)
        
        nonzero = np.where(data > 0)
        if len(nonzero[0]) > 0:
            slices_ax0 = sorted(set(nonzero[0]))
            first = slices_ax0[0]
            print(f"  {sid} ({ses}): first_thalamus_fod_slice={first}  "
                  f"(range {first}-{slices_ax0[-1]}, shape={data.shape})")
        else:
            print(f"  {sid}: NO nonzero voxels in thalamus ROI!")
