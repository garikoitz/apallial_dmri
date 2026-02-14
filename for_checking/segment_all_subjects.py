#!/usr/bin/env python3
"""
Segment graft (AP/APTR) or isocortex (CTRL) into 6 anatomical regions:
- anterior/central/posterior (X axis divided into thirds)
- dorsal/ventral (Y axis, per-slice midpoint)

For AP/APTR: uses manual_roi/flipped/roi_Graft_{hemi}_flip.mif
For CTRL: uses atlas/rois/roi_isocortex_{hemi}.nii.gz
"""

import os
import subprocess
import numpy as np
import nibabel as nib

# Base path
base_path = "/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc"

# Subject definitions by group
subjects = {
    "APTR": [
        ("sub-AP72-862-901-NmF", "ses-251210093342"),
        ("sub-AP72-862-901-R1M", "ses-251210100043"),
        ("sub-Apallial17NmF", "ses-240816091633"),
        ("sub-Apallial20R1F", "ses-240909094232"),
    ],
    "AP": [
        ("sub-596-C-OO-NmM", "ses-240911103349"),
        ("sub-769-C-O5-R1M", "ses-240911115733"),
        ("sub-769-C-OO-NmM", "ses-240911111853"),
        ("sub-875-C-OO-NM", "ses-240911100253"),
    ],
    "CTRL": [
        ("sub-129-921-R1F", "ses-260115121231"),
        ("sub-129-921-R1L1F", "ses-260116085433"),
        ("sub-129-1014-L1M", "ses-260115113814"),
        ("sub-138-1010-L1F", "ses-260115101441"),
    ],
}


def segment_roi(data, affine, output_dir, roi_name, hemi):
    """Segment an ROI into 6 parts: ant/cent/post x dorsal/ventral"""
    
    # Find voxels in the ROI
    coords = np.where(data > 0)
    
    if len(coords[0]) == 0:
        print(f"    WARNING: No voxels found in {roi_name} {hemi}!")
        return {}
    
    x_coords = coords[0]  # Anterior-Posterior
    y_coords = coords[1]  # Dorsal-Ventral
    z_coords = coords[2]  # Medial-Lateral (left-right slices)
    
    print(f"    Total voxels: {len(x_coords)}")
    print(f"    X range (ant-post): {x_coords.min()} - {x_coords.max()}")
    print(f"    Y range (dor-vent): {y_coords.min()} - {y_coords.max()}")
    print(f"    Z range (med-lat): {z_coords.min()} - {z_coords.max()}")
    
    # Divide X (anterior-posterior) into thirds
    x_min, x_max = x_coords.min(), x_coords.max()
    x_range = x_max - x_min + 1
    x_third = x_range / 3
    
    x_post_max = x_min + x_third  # Low X = posterior
    x_cent_max = x_min + 2 * x_third
    
    print(f"    Posterior: X < {x_post_max:.1f}")
    print(f"    Central: {x_post_max:.1f} <= X < {x_cent_max:.1f}")
    print(f"    Anterior: X >= {x_cent_max:.1f}")
    
    # Initialize output arrays
    segments = {
        "anterior_dorsal": np.zeros_like(data),
        "anterior_ventral": np.zeros_like(data),
        "central_dorsal": np.zeros_like(data),
        "central_ventral": np.zeros_like(data),
        "posterior_dorsal": np.zeros_like(data),
        "posterior_ventral": np.zeros_like(data),
    }
    
    # Process slice by slice in Z (medial-lateral direction)
    unique_z = np.unique(z_coords)
    
    for z in unique_z:
        # Get voxels in this Z slice
        z_mask = (z_coords == z)
        slice_x = x_coords[z_mask]
        slice_y = y_coords[z_mask]
        
        if len(slice_y) == 0:
            continue
        
        # Find dorsal/ventral midpoint for this slice
        y_mid = (slice_y.min() + slice_y.max()) / 2
        
        # Assign each voxel to appropriate segment
        for i, (x, y) in enumerate(zip(slice_x, slice_y)):
            # Determine anterior/central/posterior
            if x < x_post_max:
                ap_section = "posterior"  # Low X = posterior
            elif x < x_cent_max:
                ap_section = "central"
            else:
                ap_section = "anterior"  # High X = anterior
            
            # Determine dorsal/ventral (Y axis: higher = more ventral)
            if y >= y_mid:
                dv_section = "ventral"  # High Y = ventral
            else:
                dv_section = "dorsal"  # Low Y = dorsal
            
            segment_name = f"{ap_section}_{dv_section}"
            segments[segment_name][int(x), int(y), int(z)] = 1
    
    # Save each segment
    results = {}
    print(f"    Segment voxel counts:")
    for seg_name, seg_data in segments.items():
        count = int(np.sum(seg_data))
        print(f"      {seg_name}: {count} voxels")
        
        if count > 0:
            out_path = os.path.join(output_dir, f"{roi_name}_{hemi}_{seg_name}.nii.gz")
            out_img = nib.Nifti1Image(seg_data.astype(np.float32), affine)
            nib.save(out_img, out_path)
            results[seg_name] = out_path
    
    return results


def process_graft_subject(sub, ses, group):
    """Process AP/APTR subject with graft ROI from manual_roi/flipped/"""
    print(f"\n{'='*60}")
    print(f"Processing {group}: {sub}/{ses} (GRAFT)")
    print(f"{'='*60}")
    
    session_dir = os.path.join(base_path, sub, ses)
    dwi_dir = os.path.join(session_dir, "dwi")
    
    # Output directory
    output_dir = os.path.join(dwi_dir, "graft_segments_v2")
    os.makedirs(output_dir, exist_ok=True)
    
    for hemi in ["left", "right"]:
        print(f"\n  === Processing {hemi} hemisphere ===")
        
        # Input ROI (from flipped folder - correct orientation)
        roi_path = os.path.join(session_dir, "manual_roi", "flipped", f"roi_Graft_{hemi}_flip.mif")
        
        if not os.path.exists(roi_path):
            print(f"    WARNING: ROI not found: {roi_path}")
            continue
        
        # Convert to nifti temporarily
        tmp_nii = os.path.join(output_dir, f"tmp_graft_{hemi}.nii.gz")
        subprocess.run(["mrconvert", roi_path, tmp_nii, "-force"], check=True, 
                      capture_output=True)
        
        # Load the ROI
        img = nib.load(tmp_nii)
        data = img.get_fdata()
        affine = img.affine
        
        # Segment it
        segment_roi(data, affine, output_dir, "graft", hemi)
        
        # Clean up temp file
        os.remove(tmp_nii)
    
    print(f"\n  Segments saved to: {output_dir}")


def process_isocortex_subject(sub, ses):
    """Process CTRL subject with isocortex ROI from atlas/rois/"""
    print(f"\n{'='*60}")
    print(f"Processing CTRL: {sub}/{ses} (ISOCORTEX)")
    print(f"{'='*60}")
    
    session_dir = os.path.join(base_path, sub, ses)
    dwi_dir = os.path.join(session_dir, "dwi")
    
    # Output directory
    output_dir = os.path.join(dwi_dir, "isocortex_segments_v2")
    os.makedirs(output_dir, exist_ok=True)
    
    for hemi in ["left", "right"]:
        print(f"\n  === Processing {hemi} hemisphere ===")
        
        # Input ROI (from atlas - already in correct orientation)
        roi_path = os.path.join(dwi_dir, "atlas", "rois", f"roi_isocortex_{hemi}.nii.gz")
        
        if not os.path.exists(roi_path):
            print(f"    WARNING: ROI not found: {roi_path}")
            continue
        
        # Load the ROI
        img = nib.load(roi_path)
        data = img.get_fdata()
        affine = img.affine
        
        # Segment it
        segment_roi(data, affine, output_dir, "isocortex", hemi)
    
    print(f"\n  Segments saved to: {output_dir}")


def main():
    print("="*60)
    print("BATCH SEGMENTATION: ALL SUBJECTS")
    print("="*60)
    
    # Process APTR subjects (graft)
    print("\n\n" + "#"*60)
    print("# APTR GROUP (Graft)")
    print("#"*60)
    for sub, ses in subjects["APTR"]:
        process_graft_subject(sub, ses, "APTR")
    
    # Process AP subjects (graft)
    print("\n\n" + "#"*60)
    print("# AP GROUP (Graft)")
    print("#"*60)
    for sub, ses in subjects["AP"]:
        process_graft_subject(sub, ses, "AP")
    
    # Process CTRL subjects (isocortex)
    print("\n\n" + "#"*60)
    print("# CTRL GROUP (Isocortex)")
    print("#"*60)
    for sub, ses in subjects["CTRL"]:
        process_isocortex_subject(sub, ses)
    
    print("\n\n" + "="*60)
    print("SEGMENTATION COMPLETE!")
    print("="*60)


if __name__ == "__main__":
    main()
