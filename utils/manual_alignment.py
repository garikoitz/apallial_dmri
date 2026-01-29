#!/usr/bin/env python3
"""
Manual alignment helper for atlas-to-subject registration.

This script helps create an initial transform by:
1. Downsampling atlas to subject resolution
2. Allowing manual adjustment of translation/rotation
3. Generating initial transform for ANTs
"""

import argparse
import subprocess
import numpy as np
from pathlib import Path


def run_command(cmd, description):
    """Execute a shell command."""
    print(f"\n{description}")
    print(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        return False
    return True


def downsample_atlas(atlas_path, output_path, target_voxel_size=(0.1, 0.1, 0.3)):
    """Downsample atlas to match subject resolution approximately."""
    cmd = [
        'mrgrid', str(atlas_path),
        'regrid', str(output_path),
        '-voxel', ','.join(map(str, target_voxel_size)),
        '-force'
    ]
    return run_command(cmd, "Downsampling atlas to subject resolution...")


def apply_rotation(input_path, output_path, axis='z', angle_degrees=90):
    """
    Apply rotation around specified axis using axis permutation and flipping.
    
    Axes:
    - x: left-right (sagittal axis) - use to convert axial↔coronal
    - y: anterior-posterior (coronal axis) - use to convert axial↔sagittal  
    - z: superior-inferior (axial axis) - use to convert coronal↔sagittal
    
    For axial→coronal (subject sees dorsal-ventral instead of anterior-posterior):
    Use X-axis rotation
    """
    angle = angle_degrees % 360
    axis = axis.lower()
    
    if angle == 0:
        # No rotation, just copy
        cmd = ['mrconvert', str(input_path), str(output_path), '-force']
    elif axis == 'x':
        # Rotation around X-axis (left-right): changes Y-Z relationship
        if angle == 90 or angle == -270:
            # Axial→Coronal: swap Y,Z
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '0,2,1',  # Keep X, swap Y and Z
                '-strides', '1,2,-3',  # Flip Z
                '-force'
            ]
        elif angle == 270 or angle == -90:
            # Coronal→Axial: swap Y,Z opposite direction
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '0,2,1',  # Keep X, swap Y and Z
                '-strides', '1,-2,3',  # Flip Y
                '-force'
            ]
        elif angle == 180 or angle == -180:
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-strides', '1,-2,-3',  # Flip Y and Z
                '-force'
            ]
        else:
            print(f"Warning: Rotation angle {angle_degrees}° not supported.")
            return False
    elif axis == 'y':
        # Rotation around Y-axis (anterior-posterior): changes X-Z relationship
        if angle == 90 or angle == -270:
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '2,1,0',  # Swap X and Z, keep Y
                '-force'
            ]
        elif angle == 270 or angle == -90:
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '2,1,0',  # Swap X and Z, keep Y  
                '-force'
            ]
        elif angle == 180 or angle == -180:
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-strides', '-1,2,-3',  # Flip X and Z
                '-force'
            ]
        else:
            print(f"Warning: Rotation angle {angle_degrees}° not supported.")
            return False
    elif axis == 'z':
        # Rotation around Z-axis (superior-inferior): changes X-Y relationship
        if angle == 270 or angle == -90:
            # 90° clockwise: swap X,Y and flip X
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '1,0,2',  # Swap X and Y
                '-strides', '-1,2,3',  # Flip first axis
                '-force'
            ]
        elif angle == 180 or angle == -180:
            # 180°: flip both X and Y
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-strides', '-1,-2,3',  # Flip X and Y
                '-force'
            ]
        elif angle == 90 or angle == -270:
            # 90° counter-clockwise: swap X,Y and flip Y
            cmd = [
                'mrconvert', str(input_path), str(output_path),
                '-axes', '1,0,2',  # Swap X and Y
                '-strides', '1,-2,3',  # Flip second axis
                '-force'
            ]
        else:
            print(f"Warning: Rotation angle {angle_degrees}° not supported.")
            return False
    else:
        print(f"Warning: Unknown axis {axis}. Use 'x', 'y', or 'z'.")
        return False
    
    return run_command(cmd, f"Applying {angle_degrees}° rotation around {axis.upper()}-axis...")


def create_initial_alignment(atlas_path, subject_path, output_dir, 
                            rotation_axis='z', rotation_angle=0, flip_axis=None):
    """
    Create initial alignment of atlas to subject.
    
    Returns path to downsampled and aligned atlas.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Downsample atlas
    downsampled = output_dir / 'atlas_downsampled.mif'
    print("\n" + "=" * 80)
    print("Step 1: Downsampling atlas to subject resolution")
    print("=" * 80)
    if not downsample_atlas(atlas_path, downsampled):
        return None
    
    # Step 2: Apply rotation if needed
    rotated = output_dir / 'atlas_rotated.mif'
    print("\n" + "=" * 80)
    print("Step 2: Applying rotation")
    print("=" * 80)
    if not apply_rotation(downsampled, rotated, rotation_axis, rotation_angle):
        return None
    
    # Step 3: Apply flip if needed
    aligned = output_dir / 'atlas_manually_aligned.mif'
    if flip_axis:
        print("\n" + "=" * 80)
        print(f"Step 3: Flipping along {flip_axis}-axis")
        print("=" * 80)
        
        # Determine which stride to flip based on axis
        if flip_axis == 'x':
            strides = '-1,2,3'
        elif flip_axis == 'y':
            strides = '1,-2,3'
        elif flip_axis == 'z':
            strides = '1,2,-3'
        else:
            print(f"Warning: Unknown flip axis {flip_axis}")
            return None
        
        cmd = ['mrconvert', str(rotated), str(aligned), '-strides', strides, '-force']
        if not run_command(cmd, f"Flipping along {flip_axis}-axis..."):
            return None
    else:
        # No flip needed, just rename
        cmd = ['mrconvert', str(rotated), str(aligned), '-force']
        run_command(cmd, "Finalizing alignment...")
    
    print("\n" + "=" * 80)
    print("Manual alignment complete!")
    print("=" * 80)
    print(f"\nDownsampled atlas: {downsampled}")
    if rotation_angle != 0:
        print(f"Rotated atlas: {rotated}")
    print(f"Aligned atlas: {aligned}")
    print(f"\nView in mrview:")
    print(f"mrview {subject_path} -overlay.load {aligned} -overlay.opacity 0.5")
    
    return aligned


def main():
    parser = argparse.ArgumentParser(
        description='Manual alignment helper for atlas registration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Just downsample
  %(prog)s --downsample-only
  
  # Convert axial to coronal (rotate around X-axis)
  %(prog)s --axis x --rotate 90
  
  # Rotate in-plane 90 degrees clockwise
  %(prog)s --axis z --rotate -90
  
  # View result and iterate
  mrview <subject> -overlay.load <atlas_manually_aligned.mif> -overlay.opacity 0.5
        """
    )
    
    parser.add_argument('--atlas', type=str,
                       default='/Users/glerma/soft/mbrain/local/atlas/allen_ccf/average_template_25.nii.gz',
                       help='Path to atlas template')
    parser.add_argument('--subject', type=str,
                       default='/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_t2w_highres.mif',
                       help='Path to subject T2w image')
    parser.add_argument('--output-dir', type=str,
                       default='/Users/glerma/soft/mbrain/local/atlas/allen_ccf',
                       help='Output directory for aligned files')
    
    parser.add_argument('--axis', type=str, default='z', choices=['x', 'y', 'z'],
                       help='Rotation axis: x (L-R, for axial↔coronal), y (A-P), z (S-I, for in-plane)')
    parser.add_argument('--rotate', type=int, default=0,
                       help='Rotation angle in degrees (90, 180, 270 or -90). Negative = opposite direction')
    parser.add_argument('--flip', type=str, default=None, choices=['x', 'y', 'z'],
                       help='Flip along specified axis after rotation (x=L-R, y=A-P, z=S-I/D-V)')
    
    parser.add_argument('--downsample-only', action='store_true',
                       help='Only downsample, no rotation')
    
    args = parser.parse_args()
    
    # Run alignment
    result = create_initial_alignment(
        args.atlas, args.subject, args.output_dir,
        rotation_axis=args.axis if not args.downsample_only else 'z',
        rotation_angle=args.rotate if not args.downsample_only else 0,
        flip_axis=args.flip if not args.downsample_only else None
    )
    
    if result:
        print("\n✓ Success! Now view the alignment and adjust parameters as needed.")
    else:
        print("\n✗ Failed to create initial alignment")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
