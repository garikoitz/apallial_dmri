#!/usr/bin/env python3
"""
Create spherical ROIs from voxel coordinates and run tractography.
"""
import subprocess
import os
import yaml
import numpy as np
import nibabel as nib

def create_sphere_roi(template_path, center_voxel, radius_mm, output_path):
    """Create a spherical ROI centered at given voxel coordinates."""
    # Load template to get dimensions and affine
    img = nib.load(template_path)
    data = img.get_fdata()
    affine = img.affine
    header = img.header
    
    # Get voxel dimensions in mm
    voxel_dims = np.abs(np.diag(affine)[:3])
    print(f"  Voxel dimensions: {voxel_dims}")
    
    # Create empty mask
    mask = np.zeros(data.shape[:3], dtype=np.uint8)
    
    # Convert radius to voxels for each dimension
    radius_voxels = radius_mm / voxel_dims
    print(f"  Radius in voxels: {radius_voxels}")
    
    # Create sphere
    cx, cy, cz = center_voxel
    
    # Create coordinate grids
    x = np.arange(mask.shape[0])
    y = np.arange(mask.shape[1])
    z = np.arange(mask.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Calculate normalized distance from center (ellipsoid accounting for voxel size)
    dist = np.sqrt(
        ((X - cx) * voxel_dims[0])**2 + 
        ((Y - cy) * voxel_dims[1])**2 + 
        ((Z - cz) * voxel_dims[2])**2
    )
    
    # Set voxels within radius to 1
    mask[dist <= radius_mm] = 1
    
    print(f"  Sphere center: [{cx}, {cy}, {cz}], radius: {radius_mm}mm")
    print(f"  Non-zero voxels: {np.sum(mask)}")
    
    # Save
    out_img = nib.Nifti1Image(mask, affine, header)
    nib.save(out_img, output_path)
    print(f"  Saved: {output_path}")
    
    return output_path

def run_tractography(dwi_dir, seed_roi, target_roi, output_name, nstreamlines=5000):
    """Run tckgen tractography between two ROIs."""
    fod_file = os.path.join(dwi_dir, "fod_wm.mif")
    output_tck = os.path.join(dwi_dir, f"tracks_{output_name}.tck")
    
    # Check FOD exists
    if not os.path.exists(fod_file):
        print(f"  ERROR: FOD not found: {fod_file}")
        return None
    
    cmd = [
        "tckgen", fod_file, output_tck,
        "-algorithm", "iFOD2",
        "-seed_image", seed_roi,
        "-include", target_roi,
        "-select", str(nstreamlines),
        "-maxlength", "25",
        "-minlength", "1",
        "-cutoff", "0.05",
        "-force"
    ]
    
    print(f"  Running: tckgen -> {output_name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        return None
    
    # Check output
    if os.path.exists(output_tck):
        # Get track count
        info_result = subprocess.run(["tckinfo", output_tck], capture_output=True, text=True)
        for line in info_result.stdout.split('\n'):
            if 'count' in line.lower():
                print(f"  {line.strip()}")
        return output_tck
    return None

def main():
    # Load config
    with open("/Users/glerma/soft/mbrain/config.yaml", 'r') as f:
        config = yaml.safe_load(f)
    
    spherical_rois = config.get('spherical_rois', {})
    radius_mm = spherical_rois.get('radius_mm', 0.4)
    
    base_dir = "/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc"
    
    # Test subjects (one per group)
    test_subjects = [
        # (group, subject, session, roi_name, coords, seed_type)
        ("CTRL", "sub-129-921-R1F", "ses-260115121231", "sphere_frontal_left", [34, 108, 95], "thalamus_left"),
        ("AP", "sub-596-C-OO-NmM", "ses-240911103349", "sphere_frontal_right", [39, 110, 102], "graft"),
        ("APTR", "sub-AP72-862-901-NmF", "ses-251210093342", "sphere_frontal_right", [42, 108, 79], "graft"),
    ]
    
    generated_tracks = []
    
    for group, subject, session, roi_name, coords, seed_type in test_subjects:
        print(f"\n{'='*60}")
        print(f"Processing: {subject} ({group})")
        print(f"{'='*60}")
        
        dwi_dir = os.path.join(base_dir, subject, session, "dwi")
        template_file = os.path.join(dwi_dir, "t2w_highres_flip.mif")
        
        # Convert to nifti for nibabel
        template_nii = os.path.join(dwi_dir, "t2w_highres_flip_tmp.nii.gz")
        if not os.path.exists(template_nii):
            subprocess.run(["mrconvert", template_file, template_nii, "-force"], 
                          capture_output=True)
        
        # Create spherical ROIs directory
        sphere_roi_dir = os.path.join(dwi_dir, "spherical_rois")
        os.makedirs(sphere_roi_dir, exist_ok=True)
        
        # Create sphere ROI
        sphere_path = os.path.join(sphere_roi_dir, f"{roi_name}.nii.gz")
        print(f"\nCreating sphere ROI: {roi_name}")
        create_sphere_roi(template_nii, coords, radius_mm, sphere_path)
        
        # Get seed ROI based on group
        if group == "CTRL":
            # Use atlas thalamus
            seed_roi = os.path.join(dwi_dir, "atlas/rois/roi_thalamus_left.nii.gz")
        else:
            # Use manual graft ROI
            manual_roi_dir = os.path.join(base_dir, subject, session, "manual_roi")
            seed_roi = os.path.join(manual_roi_dir, "flipped/roi_Graft_left_flip.mif")
        
        if not os.path.exists(seed_roi):
            print(f"  WARNING: Seed ROI not found: {seed_roi}")
            continue
        
        # Run tractography: seed -> frontal sphere
        track_name = f"Frontal_{roi_name.replace('sphere_frontal_', '')}"
        tck_file = run_tractography(dwi_dir, seed_roi, sphere_path, track_name)
        
        if tck_file:
            generated_tracks.append((subject, session, tck_file))
            print(f"  SUCCESS: Generated {tck_file}")
        
        # Clean up temp file
        if os.path.exists(template_nii):
            os.remove(template_nii)
    
    print(f"\n{'='*60}")
    print("SUMMARY: Generated tracks for MATLAB cleaning")
    print(f"{'='*60}")
    for subj, sess, tck in generated_tracks:
        print(f"  {subj}: {os.path.basename(tck)}")
    
    return generated_tracks

if __name__ == "__main__":
    main()
