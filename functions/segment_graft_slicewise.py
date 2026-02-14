#!/usr/bin/env python3
"""
Segment graft ROI into sub-regions for tractography analysis.

Supports two segmentation modes:
- 6 segments (3x2): anterior/central/posterior x dorsal/ventral
- 9 segments (3x3): anterior/central/posterior x dorsal/medial/ventral

Input: Graft ROI from Slicer in T2w_highres space (manual_roi/flipped/)
Output: Segmented ROIs in both T2w space and resampled to FOD space

Usage: 
    # 6-segment mode (3x2 grid) - default
    python functions/segment_graft_slicewise.py --subject sub-596-C-OO-NmM --session ses-240911103349

    # 9-segment mode (3x3 grid)
    python functions/segment_graft_slicewise.py --subject sub-596-C-OO-NmM --session ses-240911103349 --segments 9

    # Left hemisphere only, 9 segments
    python functions/segment_graft_slicewise.py --subject sub-596-C-OO-NmM --session ses-240911103349 --segments 9 --hemisphere left
"""

import nibabel as nib
import numpy as np
import os
import subprocess
import argparse

# Set pipeline to flip1
os.environ['MBRAIN_PIPELINE'] = 'flip1'

from config_loader import load_config

# Load config
config = load_config()

# Parse arguments
parser = argparse.ArgumentParser(description="Segment graft ROI into sub-regions (6 or 9 segments)")
parser.add_argument("--subject", required=True, help="Subject ID (e.g., sub-596-C-OO-NmM)")
parser.add_argument("--session", required=True, help="Session ID (e.g., ses-240911103349)")
parser.add_argument("--hemisphere", choices=["left", "right", "both"], default="both",
                    help="Which hemisphere to process (default: both)")
parser.add_argument("--segments", type=int, choices=[6, 9], default=6,
                    help="Number of segments: 6 (3x2 grid) or 9 (3x3 grid). Default: 6")
args = parser.parse_args()

subject = args.subject
session = args.session
num_segments = args.segments

# Use flip1 derivatives folder from config
derivatives_dir = config['paths']['derivatives_dir']
base_dir = os.path.join(derivatives_dir, subject, session)
dwi_dir = os.path.join(base_dir, "dwi")

# Output directories - T2w space and FOD space
output_dir_t2w = os.path.join(dwi_dir, "graft_segments_t2w")
output_dir_fod = os.path.join(dwi_dir, "graft_segments_fod")
os.makedirs(output_dir_t2w, exist_ok=True)
os.makedirs(output_dir_fod, exist_ok=True)

# FOD template for resampling (flip1 pipeline uses wm_fod_norm.mif)
fod_template = os.path.join(dwi_dir, "wm_fod_norm.mif")
if not os.path.exists(fod_template):
    print(f"ERROR: FOD template not found: {fod_template}")
    print("Cannot resample to FOD space. Exiting.")
    exit(1)

# Determine which hemispheres to process
hemispheres = ["left", "right"] if args.hemisphere == "both" else [args.hemisphere]

# Segment names based on mode
AP_SECTIONS = ["anterior", "central", "posterior"]
if num_segments == 6:
    DV_SECTIONS = ["dorsal", "ventral"]
    prefix = "seg6"
    print(f"\n=== 6-SEGMENT MODE (3x2 grid: anterior/central/posterior x dorsal/ventral) ===")
else:
    DV_SECTIONS = ["dorsal", "medial", "ventral"]
    prefix = "seg9"
    print(f"\n=== 9-SEGMENT MODE (3x3 grid: anterior/central/posterior x dorsal/medial/ventral) ===")

print(f"Subject: {subject}")
print(f"Session: {session}")
print(f"Output prefix: {prefix}_")


def resample_to_fod(input_path, output_path, template_path):
    """Resample a NIfTI file to FOD space using nearest-neighbor interpolation."""
    cmd = f"mrgrid '{input_path}' regrid '{output_path}' -template '{template_path}' -interp nearest -force -quiet"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.returncode == 0


for hemi in hemispheres:
    print(f"\n{'='*60}")
    print(f"Processing {hemi} hemisphere")
    print(f"{'='*60}")
    
    # Input ROI (flip1: manual_roi/roi_Graft_{hemi}_flip.mif)
    roi_path = os.path.join(base_dir, "manual_roi", f"roi_Graft_{hemi}_flip.mif")
    
    if not os.path.exists(roi_path):
        print(f"  WARNING: ROI not found: {roi_path}")
        continue
    
    # Convert to nifti temporarily
    tmp_nii = os.path.join(output_dir_t2w, f"tmp_graft_{hemi}.nii.gz")
    subprocess.run(["mrconvert", roi_path, tmp_nii, "-force", "-quiet"], check=True)
    
    # Load the ROI
    img = nib.load(tmp_nii)
    data = img.get_fdata()
    affine = img.affine
    
    # Find voxels in the ROI
    coords = np.where(data > 0)
    if len(coords[0]) == 0:
        print(f"  WARNING: No voxels found in {hemi} graft ROI")
        os.remove(tmp_nii)
        continue
    
    x_coords = coords[0]  # Anterior-Posterior (axis 0)
    y_coords = coords[1]  # Dorsal-Ventral (axis 1)
    z_coords = coords[2]  # Medial-Lateral (axis 2)
    
    print(f"  Total voxels: {len(x_coords)}")
    print(f"  X range (ant-post): {x_coords.min()} - {x_coords.max()}")
    print(f"  Y range (dor-vent): {y_coords.min()} - {y_coords.max()}")
    print(f"  Z range (med-lat): {z_coords.min()} - {z_coords.max()}")
    
    # Divide X (anterior-posterior) into thirds
    x_min, x_max = x_coords.min(), x_coords.max()
    x_range = x_max - x_min + 1
    x_third = x_range / 3
    
    x_ant_max = x_min + x_third      # boundary between posterior and central
    x_cent_max = x_min + 2 * x_third  # boundary between central and anterior
    
    print(f"\n  X-axis (Anterior-Posterior) divisions:")
    print(f"    Posterior: X < {x_ant_max:.1f}  (LOW X = back of brain)")
    print(f"    Central:   {x_ant_max:.1f} <= X < {x_cent_max:.1f}")
    print(f"    Anterior:  X >= {x_cent_max:.1f}  (HIGH X = nose)")
    
    # Initialize output arrays for all segments
    segments = {}
    for ap in AP_SECTIONS:
        for dv in DV_SECTIONS:
            segments[f"{ap}_{dv}"] = np.zeros_like(data)
    
    # Process slice by slice in Z (medial-lateral direction)
    unique_z = np.unique(z_coords)
    print(f"\n  Processing {len(unique_z)} Z slices...")
    
    for z in unique_z:
        # Get voxels in this Z slice
        z_mask = (z_coords == z)
        slice_x = x_coords[z_mask]
        slice_y = y_coords[z_mask]
        
        if len(slice_y) == 0:
            continue
        
        # Divide Y (dorsal-ventral) for THIS SLICE
        y_min_slice = slice_y.min()
        y_max_slice = slice_y.max()
        y_range_slice = y_max_slice - y_min_slice + 1
        
        if num_segments == 6:
            # 2-way split: dorsal/ventral at midpoint
            y_mid = (y_min_slice + y_max_slice) / 2
        else:
            # 3-way split: dorsal/medial/ventral in thirds
            y_third = y_range_slice / 3
            y_dorsal_max = y_min_slice + y_third
            y_medial_max = y_min_slice + 2 * y_third
        
        # Assign each voxel to appropriate segment
        for i, (x, y) in enumerate(zip(slice_x, slice_y)):
            # Determine anterior/central/posterior
            # FLIP1 AXIS 0 (X): HIGH index = ANTERIOR (nose), LOW index = POSTERIOR (back)
            # Verified: Graft X=24.1 > Thalamus X=12.4 (Graft is anterior to Thalamus)
            if x < x_ant_max:
                ap_section = "posterior"
            elif x < x_cent_max:
                ap_section = "central"
            else:
                ap_section = "anterior"
            
            # Determine dorsal/ventral (or dorsal/medial/ventral)
            # Note: In flip1 orientation, LOW Y = dorsal, HIGH Y = ventral
            if num_segments == 6:
                if y < y_mid:
                    dv_section = "dorsal"
                else:
                    dv_section = "ventral"
            else:
                if y < y_dorsal_max:
                    dv_section = "dorsal"
                elif y < y_medial_max:
                    dv_section = "medial"
                else:
                    dv_section = "ventral"
            
            segment_name = f"{ap_section}_{dv_section}"
            segments[segment_name][int(x), int(y), int(z)] = 1
    
    # Save each segment in T2w space and resample to FOD space
    print(f"\n  Segment voxel counts for {hemi}:")
    total_voxels = 0
    for ap in AP_SECTIONS:
        for dv in DV_SECTIONS:
            seg_name = f"{ap}_{dv}"
            seg_data = segments[seg_name]
            count = int(np.sum(seg_data))
            total_voxels += count
            
            if count > 0:
                # Output filename with prefix - T2w in nii.gz, FOD in mif
                out_filename_t2w = f"{prefix}_graft_{hemi}_{seg_name}.nii.gz"
                out_filename_fod = f"{prefix}_graft_{hemi}_{seg_name}_fod.mif"
                
                # Save in T2w space
                out_path_t2w = os.path.join(output_dir_t2w, out_filename_t2w)
                out_img = nib.Nifti1Image(seg_data.astype(np.float32), affine)
                nib.save(out_img, out_path_t2w)
                
                # Resample to FOD space as .mif
                out_path_fod = os.path.join(output_dir_fod, out_filename_fod)
                if resample_to_fod(out_path_t2w, out_path_fod, fod_template):
                    print(f"    {out_filename_fod}: {count} voxels ✓ (T2w + FOD)")
                else:
                    print(f"    {out_filename_fod}: {count} voxels ✓ (T2w only, FOD resample failed)")
            else:
                print(f"    {prefix}_graft_{hemi}_{seg_name}: 0 voxels ⚠ EMPTY")
    
    print(f"  Total voxels saved: {total_voxels}")
    
    # Clean up temp file
    os.remove(tmp_nii)

print(f"\n{'='*60}")
print(f"DONE!")
print(f"{'='*60}")
print(f"\nOutputs saved to:")
print(f"  T2w space: {output_dir_t2w}/{prefix}_graft_*.nii.gz")
print(f"  FOD space: {output_dir_fod}/{prefix}_graft_*.nii.gz")
print(f"\n=== VISUALIZATION ===")
print(f"mrview {dwi_dir}/t2w_highres_flip.mif \\")
print(f"  -overlay.load {output_dir_t2w}/{prefix}_graft_left_anterior_dorsal.nii.gz -overlay.opacity 0.6 -overlay.colour 1,0,0 &")
print(f"\n=== TRACTOGRAPHY (use FOD space ROIs) ===")
print(f"# After running this, update tractography_segments.py to use the new prefix")
