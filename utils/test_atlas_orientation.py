#!/usr/bin/env python3
"""
Quick script to test different atlas orientations without re-running registration.
This creates multiple test versions of the atlas with different axis permutations and flips.
"""

import subprocess
from pathlib import Path
import sys

# Paths
atlas_dir = Path('/Users/glerma/soft/mbrain/local/atlas/allen_ccf')
atlas_downsampled = atlas_dir / 'atlas_template_downsampled.mif'
subject_t2w = Path('/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_t2w_highres.mif')

# Test configurations: (name, axes, flip_axis)
# axes: permutation string like "0,1,2" or "2,1,0"
# flip_axis: which axis to flip after permutation (0, 1, 2, or None)
tests = [
    ("original", "0,1,2", None),
    ("axes_210", "2,1,0", None),
    ("axes_210_flip0", "2,1,0", 0),
    ("axes_210_flip1", "2,1,0", 1),
    ("axes_210_flip2", "2,1,0", 2),
    ("axes_021", "0,2,1", None),
    ("axes_021_flip1", "0,2,1", 1),
    ("axes_102", "1,0,2", None),
    ("axes_102_flip0", "1,0,2", 0),
    ("axes_201", "2,0,1", None),
    ("axes_120", "1,2,0", None),
    ("axes_012_flip0", "0,1,2", 0),
    ("axes_012_flip1", "0,1,2", 1),
    ("axes_012_flip2", "0,1,2", 2),
]

print("Testing different atlas orientations...")
print(f"Source: {atlas_downsampled}")
print(f"Output directory: {atlas_dir}/orientation_tests/")
print()

output_dir = atlas_dir / 'orientation_tests'
output_dir.mkdir(exist_ok=True)

for test_name, axes, flip_axis in tests:
    print(f"Creating {test_name}...")
    
    # Step 1: Apply axis permutation
    temp_file = output_dir / f"{test_name}_temp.mif"
    cmd = [
        'mrconvert',
        str(atlas_downsampled),
        str(temp_file),
        '-axes', axes,
        '-strides', '1,2,3',
        '-force'
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    
    # Step 2: Apply flip if needed
    if flip_axis is not None:
        # Convert to NIfTI for numpy flip
        temp_nii = output_dir / f"{test_name}_temp.nii.gz"
        subprocess.run([
            'mrconvert', str(temp_file), str(temp_nii), '-force'
        ], check=True, capture_output=True)
        
        # Flip with Python
        import nibabel as nib
        import numpy as np
        img = nib.load(str(temp_nii))
        data = np.flip(img.get_fdata(), axis=flip_axis)
        nib.save(nib.Nifti1Image(data, img.affine, img.header), str(temp_nii))
        
        # Convert back to MIF
        output_file = output_dir / f"{test_name}.mif"
        subprocess.run([
            'mrconvert', str(temp_nii), str(output_file),
            '-strides', '1,2,3', '-force'
        ], check=True, capture_output=True)
        
        # Clean up
        temp_nii.unlink()
    else:
        # Just rename
        output_file = output_dir / f"{test_name}.mif"
        temp_file.rename(output_file)
    
    print(f"  ✓ Created {output_file.name}")

print("\n" + "="*80)
print("All test files created!")
print("="*80)
print("\nTo view each test, run:")
print(f"\nmrview {subject_t2w} \\")
print(f"  -overlay.load {output_dir}/TEST_NAME.mif \\")
print(f"  -overlay.opacity 0.5")
print("\nReplace TEST_NAME with one of:")
for test_name, _, _ in tests:
    print(f"  - {test_name}")
print("\nOr view all at once:")
all_overlays = " \\\n  ".join([f"-overlay.load {output_dir}/{name}.mif" for name, _, _ in tests])
print(f"\nmrview {subject_t2w} \\")
print(f"  {all_overlays} \\")
print(f"  -overlay.opacity 0.3")
