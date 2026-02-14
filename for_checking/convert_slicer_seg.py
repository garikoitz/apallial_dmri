#!/usr/bin/env python3
"""
Convert 3D Slicer segmentation NRRD to individual NIfTI ROIs.

This handles .seg.nrrd files from 3D Slicer and extracts each segment
as a separate binary NIfTI mask, properly aligned to the reference space.
"""

import argparse
import numpy as np
import nrrd
import nibabel as nib
from pathlib import Path


def parse_segment_info(header):
    """Extract segment names and label values from NRRD header."""
    segments = {}
    i = 0
    while True:
        name_key = f'Segment{i}_Name'
        label_key = f'Segment{i}_LabelValue'
        
        if name_key not in header:
            break
        
        name = header[name_key]
        label = int(header[label_key])
        segments[label] = name
        i += 1
    
    return segments


def nrrd_to_nifti_affine(header, adjust_origin=True):
    """Convert NRRD space directions and origin to NIfTI affine matrix.
    
    Args:
        header: NRRD header dictionary
        adjust_origin: If True, shift origin by half voxel to account for 
                      1-based (Slicer) vs 0-based (Python) indexing
    """
    # Get space directions (3x3 matrix stored as list of 3 vectors)
    space_dirs = np.array(header['space directions'])
    
    # Get origin
    origin = np.array(header['space origin'])
    
    # Build 4x4 affine matrix BEFORE any adjustments
    affine = np.eye(4)
    affine[:3, :3] = space_dirs.T  # Transpose for column-major
    affine[:3, 3] = origin
    
    # Handle LPS to RAS conversion
    space = header.get('space', 'left-posterior-superior')
    if 'left-posterior-superior' in space.lower() or 'lps' in space.lower():
        # LPS to RAS: flip x and y
        affine[0, :] *= -1
        affine[1, :] *= -1
    
    # Adjust origin for 1-based vs 0-based indexing AFTER coordinate flip
    # Slicer uses 1-based, NIfTI uses 0-based
    # Shift origin by one full voxel along each axis
    if adjust_origin:
        # Get voxel sizes from the affine matrix (after any flips)
        voxel_sizes = np.sqrt(np.sum(affine[:3, :3]**2, axis=0))
        # Shift by one voxel in the positive direction
        affine[:3, 3] = affine[:3, 3] + affine[:3, :3] @ np.array([1, 1, 1])
    
    return affine


def convert_slicer_seg_to_nifti(input_nrrd, output_dir=None, reference_image=None):
    """
    Convert 3D Slicer segmentation to individual NIfTI ROIs.
    
    Args:
        input_nrrd: Path to .seg.nrrd file
        output_dir: Output directory (default: same as input)
        reference_image: Optional reference NIfTI to match geometry
    
    Returns:
        List of created NIfTI files
    """
    input_path = Path(input_nrrd)
    
    if output_dir is None:
        output_dir = input_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading: {input_path}")
    
    # Read NRRD
    data, header = nrrd.read(str(input_path))
    
    print(f"  Shape: {data.shape}")
    print(f"  Data type: {data.dtype}")
    print(f"  Unique labels: {np.unique(data)}")
    
    # Parse segment information
    segments = parse_segment_info(header)
    print(f"\n  Segments found:")
    for label, name in segments.items():
        count = np.sum(data == label)
        print(f"    Label {label}: '{name}' ({count} voxels)")
    
    # Get affine transformation
    affine = nrrd_to_nifti_affine(header)
    
    # If reference image provided, use its affine
    if reference_image:
        ref_path = Path(reference_image)
        if ref_path.exists():
            ref_img = nib.load(str(ref_path))
            print(f"\n  Using reference geometry from: {ref_path.name}")
            # Note: We keep the data in its native space but could resample if needed
    
    # Create combined labelmap NIfTI
    combined_path = output_dir / f"{input_path.stem.replace('.seg', '')}_labels.nii.gz"
    combined_img = nib.Nifti1Image(data.astype(np.uint8), affine)
    nib.save(combined_img, str(combined_path))
    print(f"\n  Saved combined labelmap: {combined_path.name}")
    
    # Create individual binary ROIs
    created_files = [combined_path]
    
    for label, name in segments.items():
        if label == 0:
            continue  # Skip background
        
        # Create binary mask
        binary_mask = (data == label).astype(np.uint8)
        
        # Clean up name for filename
        safe_name = name.replace(' ', '_').replace('/', '-')
        roi_path = output_dir / f"roi_{safe_name}.nii.gz"
        
        roi_img = nib.Nifti1Image(binary_mask, affine)
        nib.save(roi_img, str(roi_path))
        
        print(f"  Saved ROI: {roi_path.name}")
        created_files.append(roi_path)
    
    return created_files


def main():
    parser = argparse.ArgumentParser(
        description="Convert 3D Slicer segmentation NRRD to NIfTI ROIs"
    )
    parser.add_argument("input", help="Input .seg.nrrd file")
    parser.add_argument("-o", "--output-dir", help="Output directory")
    parser.add_argument("-r", "--reference", help="Reference NIfTI for geometry")
    parser.add_argument("--visualize", action="store_true", 
                       help="Open mrview after conversion")
    
    args = parser.parse_args()
    
    files = convert_slicer_seg_to_nifti(
        args.input, 
        args.output_dir, 
        args.reference
    )
    
    if args.visualize and files:
        import subprocess
        # Find a suitable background image
        input_dir = Path(args.input).parent.parent
        fa_candidates = list(input_dir.glob("dwi/fa.mif")) + list(input_dir.glob("dwi/*fa*.nii.gz"))
        
        if fa_candidates:
            bg_image = fa_candidates[0]
        else:
            bg_image = files[0]  # Use the labelmap itself
        
        # Build mrview command
        cmd = ["mrview", str(bg_image)]
        
        # Add ROIs as overlays
        for f in files[1:]:  # Skip combined labelmap for overlay
            cmd.extend(["-overlay.load", str(f), "-overlay.opacity", "0.5"])
        
        print(f"\nLaunching mrview...")
        subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
