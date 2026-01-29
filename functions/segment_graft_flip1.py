#!/usr/bin/env python3
"""
Segment graft ROI into sub-regions for flip1 pipeline.
Creates both 6-segment (3x2) and 9-segment (3x3) versions.

Usage:
    python functions/segment_graft_flip1.py --subject sub-Apallial17NmF --session ses-240816091633
"""

import nibabel as nib
import numpy as np
import subprocess
import os
import argparse

# FLIP1 AXIS ORIENTATION:
# AXIS 0: ANTERIOR-POSTERIOR (low=anterior, high=posterior)
# AXIS 1: DORSAL-VENTRAL (low=dorsal, high=ventral)
# AXIS 2: LEFT-RIGHT (medial-lateral)

def segment_graft_6(data, affine, hemi):
    """Create 6 segments (3x2): anterior/central/posterior x dorsal/ventral"""
    nonzero = np.where(data > 0)
    if len(nonzero[0]) == 0:
        return {}
    
    ap_min, ap_max = nonzero[0].min(), nonzero[0].max()
    dv_min, dv_max = nonzero[1].min(), nonzero[1].max()
    
    ap_range = ap_max - ap_min + 1
    ap_third = ap_range // 3
    ap_boundaries = [ap_min, ap_min + ap_third, ap_min + 2*ap_third, ap_max + 1]
    
    dv_range = dv_max - dv_min + 1
    dv_half = dv_range // 2
    dv_boundaries = [dv_min, dv_min + dv_half, dv_max + 1]
    
    segments = {}
    ap_names = ["anterior", "central", "posterior"]
    dv_names = ["dorsal", "ventral"]
    
    for i, ap_name in enumerate(ap_names):
        for j, dv_name in enumerate(dv_names):
            seg_name = f"graft_{hemi}_{ap_name}_{dv_name}"
            mask = np.zeros_like(data)
            
            ap_start, ap_end = ap_boundaries[i], ap_boundaries[i+1]
            dv_start, dv_end = dv_boundaries[j], dv_boundaries[j+1]
            
            for x in range(ap_start, ap_end):
                for y in range(dv_start, dv_end):
                    mask[x, y, :] = data[x, y, :]
            
            voxel_count = np.sum(mask > 0)
            segments[seg_name] = (mask, voxel_count)
    
    return segments


def segment_graft_9(data, affine, hemi):
    """Create 9 segments (3x3): anterior/central/posterior x dorsal/medial/ventral"""
    nonzero = np.where(data > 0)
    if len(nonzero[0]) == 0:
        return {}
    
    ap_min, ap_max = nonzero[0].min(), nonzero[0].max()
    dv_min, dv_max = nonzero[1].min(), nonzero[1].max()
    
    ap_range = ap_max - ap_min + 1
    ap_third = ap_range // 3
    ap_boundaries = [ap_min, ap_min + ap_third, ap_min + 2*ap_third, ap_max + 1]
    
    dv_range = dv_max - dv_min + 1
    dv_third = dv_range // 3
    dv_boundaries = [dv_min, dv_min + dv_third, dv_min + 2*dv_third, dv_max + 1]
    
    segments = {}
    ap_names = ["anterior", "central", "posterior"]
    dv_names = ["dorsal", "medial", "ventral"]
    
    for i, ap_name in enumerate(ap_names):
        for j, dv_name in enumerate(dv_names):
            seg_name = f"graft_{hemi}_{ap_name}_{dv_name}"
            mask = np.zeros_like(data)
            
            ap_start, ap_end = ap_boundaries[i], ap_boundaries[i+1]
            dv_start, dv_end = dv_boundaries[j], dv_boundaries[j+1]
            
            for x in range(ap_start, ap_end):
                for y in range(dv_start, dv_end):
                    mask[x, y, :] = data[x, y, :]
            
            voxel_count = np.sum(mask > 0)
            segments[seg_name] = (mask, voxel_count)
    
    return segments


def main():
    parser = argparse.ArgumentParser(description="Segment graft ROI for flip1 pipeline")
    parser.add_argument("--subject", required=True, help="Subject ID")
    parser.add_argument("--session", required=True, help="Session ID")
    args = parser.parse_args()
    
    subject = args.subject
    session = args.session
    
    # Paths for flip1 pipeline
    base_dir = f"/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1/{subject}/{session}"
    manual_roi_dir = os.path.join(base_dir, "manual_roi")
    output_dir = os.path.join(base_dir, "dwi", "graft_segments_t2w")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Segmenting Graft ROIs for {subject}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}")
    
    # Process both hemispheres
    for hemi in ["left", "right"]:
        print(f"\n--- Processing {hemi} hemisphere ---")
        
        roi_mif = os.path.join(manual_roi_dir, f"roi_Graft_{hemi}_flip.mif")
        if not os.path.exists(roi_mif):
            print(f"  WARNING: ROI not found: {roi_mif}")
            continue
        
        roi_nii = os.path.join(output_dir, f"temp_graft_{hemi}.nii.gz")
        subprocess.run(f"mrconvert '{roi_mif}' '{roi_nii}' -force -quiet", shell=True)
        
        img = nib.load(roi_nii)
        data = img.get_fdata()
        affine = img.affine
        
        total_voxels = int(np.sum(data > 0))
        print(f"  Total voxels in Graft_{hemi}: {total_voxels}")
        
        # 6-segment version
        print(f"\n  6 segments (3x2 grid):")
        segments_6 = segment_graft_6(data, affine, hemi)
        for seg_name, (mask, count) in segments_6.items():
            out_path = os.path.join(output_dir, f"seg6_{seg_name}.nii.gz")
            nib.save(nib.Nifti1Image(mask.astype(np.float32), affine), out_path)
            print(f"    {seg_name}: {int(count)} voxels")
        
        # 9-segment version
        print(f"\n  9 segments (3x3 grid):")
        segments_9 = segment_graft_9(data, affine, hemi)
        for seg_name, (mask, count) in segments_9.items():
            out_path = os.path.join(output_dir, f"seg9_{seg_name}.nii.gz")
            nib.save(nib.Nifti1Image(mask.astype(np.float32), affine), out_path)
            print(f"    {seg_name}: {int(count)} voxels")
        
        os.remove(roi_nii)
    
    print(f"\n{'='*60}")
    print(f"Done! Files saved to: {output_dir}")
    print(f"{'='*60}")
    
    # List files
    print("\nCreated files:")
    for f in sorted(os.listdir(output_dir)):
        print(f"  {f}")


if __name__ == "__main__":
    main()
