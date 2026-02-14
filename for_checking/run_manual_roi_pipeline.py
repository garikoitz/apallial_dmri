#!/usr/bin/env python3
"""
Complete pipeline for manual ROI processing and tractography.

Steps:
1. Convert 3D Slicer .seg.nrrd files to NIfTI ROIs
2. Flip ROIs to match T2w_highres_flip orientation
3. Transform ROIs to DWI space using existing coregistration
4. Generate tractography from graft to target regions

Usage:
    python run_manual_roi_pipeline.py                    # Process all subjects with NRRD files
    python run_manual_roi_pipeline.py --subject sub-XXX  # Process specific subject
    python run_manual_roi_pipeline.py --step convert     # Run only conversion step
"""

import subprocess
import argparse
from pathlib import Path
import numpy as np
import sys


def find_subjects_with_nrrd(deriv_dir):
    """Find all subjects with .seg.nrrd files in manual_roi folder."""
    subjects = []
    for nrrd_file in deriv_dir.glob("*/*/manual_roi/*.seg.nrrd"):
        subject = nrrd_file.parent.parent.parent.name
        session = nrrd_file.parent.parent.name
        subjects.append({
            'subject': subject,
            'session': session,
            'nrrd_file': nrrd_file,
            'manual_roi_dir': nrrd_file.parent,
            'dwi_dir': nrrd_file.parent.parent / 'dwi'
        })
    return subjects


def convert_nrrd_to_nifti(nrrd_file, output_dir):
    """Convert 3D Slicer segmentation NRRD to individual NIfTI ROIs."""
    print(f"\n  Converting: {nrrd_file.name}")
    
    # Import here to avoid loading if not needed
    import nrrd
    import nibabel as nib
    
    # Read NRRD
    data, header = nrrd.read(str(nrrd_file))
    
    print(f"    Shape: {data.shape}")
    print(f"    Unique labels: {np.unique(data)}")
    
    # Parse segment information
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
    
    print(f"    Segments found:")
    for label, name in segments.items():
        count = np.sum(data == label)
        print(f"      Label {label}: '{name}' ({count} voxels)")
    
    # Build affine matrix from NRRD header
    space_dirs = np.array(header['space directions'])
    origin = np.array(header['space origin'])
    
    affine = np.eye(4)
    affine[:3, :3] = space_dirs.T
    affine[:3, 3] = origin
    
    # LPS to RAS conversion
    space = header.get('space', 'left-posterior-superior')
    if 'left-posterior-superior' in space.lower() or 'lps' in space.lower():
        affine[0, :] *= -1
        affine[1, :] *= -1
    
    # Adjust origin for indexing
    affine[:3, 3] = affine[:3, 3] + affine[:3, :3] @ np.array([1, 1, 1])
    
    # Create individual binary ROIs
    created_files = []
    
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
        
        print(f"    ✓ Saved: {roi_path.name}")
        created_files.append(roi_path)
    
    return created_files


def flip_roi(roi_file, output_file):
    """Flip ROI to match T2w_highres_flip orientation."""
    # Step 1: Axis permutation (2,1,0 for 3D)
    temp_file = output_file.parent / f"{output_file.stem}_perm.mif"
    
    cmd_perm = [
        "mrconvert", str(roi_file), str(temp_file),
        "-axes", "2,1,0", "-strides", "1,2,3", "-force"
    ]
    subprocess.run(cmd_perm, check=True, capture_output=True)
    
    # Step 2: Axis flip
    cmd_flip = [
        "mrtransform", str(temp_file), str(output_file),
        "-flip", "1", "-force"
    ]
    subprocess.run(cmd_flip, check=True, capture_output=True)
    
    # Clean up temp file
    temp_file.unlink()
    
    return output_file


def transform_roi_to_dwi(roi_flipped, transform_mat, roi_dwi_output, dwi_reference):
    """Transform ROI from T2w highres space to DWI space."""
    cmd = [
        "mrtransform",
        str(roi_flipped),
        str(roi_dwi_output),
        "-linear", str(transform_mat),
        "-template", str(dwi_reference),
        "-interp", "nearest",
        "-force"
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return roi_dwi_output


def run_tractography(fod_file, seed_roi, include_roi, output_tck, select=5000, maxlength=15):
    """Run tractography with iFOD2 algorithm."""
    cmd = [
        "tckgen",
        str(fod_file),
        str(output_tck),
        "-algorithm", "iFOD2",
        "-seed_image", str(seed_roi),
        "-include", str(include_roi),
        "-select", str(select),
        "-maxlength", str(maxlength),
        "-force"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"    ⚠ tckgen warning/error: {result.stderr[:200] if result.stderr else 'unknown'}")
    
    # Check track count
    info_result = subprocess.run(
        ["tckinfo", str(output_tck)],
        capture_output=True, text=True
    )
    
    track_count = 0
    for line in info_result.stdout.split('\n'):
        if 'count:' in line.lower():
            track_count = int(line.split(':')[1].strip())
            break
    
    return track_count


def run_sift(tracks_input, fod_file, tracks_output, term_number=50):
    """Apply SIFT filtering."""
    cmd = [
        "tcksift",
        str(tracks_input),
        str(fod_file),
        str(tracks_output),
        "-term_number", str(term_number),
        "-force"
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return tracks_output


def ensure_transform_exists(dwi_dir):
    """Ensure anat2dwi transform exists (invert dwi2anat if needed)."""
    anat2dwi = dwi_dir / "anat2dwi_mrtrix.txt"
    dwi2anat = dwi_dir / "dwi2anat_mrtrix.txt"
    
    if anat2dwi.exists():
        return anat2dwi
    
    if dwi2anat.exists():
        print("  Inverting dwi2anat transform...")
        subprocess.run([
            "transformcalc",
            str(dwi2anat),
            "invert",
            str(anat2dwi)
        ], check=True)
        return anat2dwi
    
    return None


def process_subject(subject_info, steps=['all'], dry_run=False):
    """Process a single subject through the pipeline."""
    subject = subject_info['subject']
    session = subject_info['session']
    nrrd_file = subject_info['nrrd_file']
    manual_roi_dir = subject_info['manual_roi_dir']
    dwi_dir = subject_info['dwi_dir']
    
    print(f"\n{'='*80}")
    print(f"Processing: {subject} / {session}")
    print(f"{'='*80}")
    print(f"  NRRD file: {nrrd_file.name}")
    
    if dry_run:
        print("  [DRY RUN - no changes will be made]")
        return True
    
    # Check prerequisites
    fod_file = dwi_dir / "fod_wm.mif"
    dwi_reference = dwi_dir / "dwi_flip_gradcheck_denoise_bias_anat.mif"
    
    if not fod_file.exists():
        print(f"  ❌ FOD file not found: {fod_file}")
        print("     → Need to run preprocessing first")
        return False
    
    if not dwi_reference.exists():
        # Try alternative
        dwi_reference = dwi_dir / "fa.mif"
        if not dwi_reference.exists():
            print(f"  ❌ DWI reference not found")
            return False
    
    # Create output directories
    flipped_dir = manual_roi_dir / "flipped"
    roi_dwi_dir = manual_roi_dir / "dwi"
    flipped_dir.mkdir(exist_ok=True)
    roi_dwi_dir.mkdir(exist_ok=True)
    
    run_all = 'all' in steps
    
    # ===== STEP 1: Convert NRRD to NIfTI =====
    if run_all or 'convert' in steps:
        print(f"\n--- Step 1: Convert NRRD to NIfTI ---")
        try:
            roi_files = convert_nrrd_to_nifti(nrrd_file, manual_roi_dir)
            print(f"  ✓ Created {len(roi_files)} ROI files")
        except Exception as e:
            print(f"  ❌ Conversion failed: {e}")
            return False
    else:
        # Find existing ROI files
        roi_files = list(manual_roi_dir.glob("roi_*.nii.gz"))
        print(f"\n--- Step 1: Skipped (found {len(roi_files)} existing ROI files) ---")
    
    if not roi_files:
        roi_files = list(manual_roi_dir.glob("roi_*.nii.gz"))
    
    if not roi_files:
        print("  ❌ No ROI files found")
        return False
    
    # ===== STEP 2: Flip ROIs =====
    if run_all or 'flip' in steps:
        print(f"\n--- Step 2: Flip ROIs to match T2w orientation ---")
        flipped_files = []
        for roi_file in roi_files:
            roi_name = roi_file.stem.replace('.nii', '')
            output_file = flipped_dir / f"{roi_name}_flip.mif"
            try:
                flip_roi(roi_file, output_file)
                flipped_files.append(output_file)
                print(f"    ✓ {output_file.name}")
            except Exception as e:
                print(f"    ❌ Failed to flip {roi_file.name}: {e}")
        print(f"  ✓ Flipped {len(flipped_files)} ROIs")
    else:
        flipped_files = list(flipped_dir.glob("roi_*_flip.mif"))
        print(f"\n--- Step 2: Skipped (found {len(flipped_files)} flipped ROIs) ---")
    
    # ===== STEP 3: Transform to DWI space =====
    if run_all or 'transform' in steps:
        print(f"\n--- Step 3: Transform ROIs to DWI space ---")
        
        transform_mat = ensure_transform_exists(dwi_dir)
        if not transform_mat:
            print("  ❌ No transform matrix found. Run coregistration first.")
            return False
        
        print(f"  Using transform: {transform_mat.name}")
        
        flipped_files = list(flipped_dir.glob("roi_*_flip.mif"))
        dwi_files = []
        for flipped_roi in flipped_files:
            roi_name = flipped_roi.stem.replace('_flip', '')
            output_file = roi_dwi_dir / f"{roi_name}_dwi.mif"
            try:
                transform_roi_to_dwi(flipped_roi, transform_mat, output_file, dwi_reference)
                dwi_files.append(output_file)
                print(f"    ✓ {output_file.name}")
            except Exception as e:
                print(f"    ❌ Failed to transform {flipped_roi.name}: {e}")
        print(f"  ✓ Transformed {len(dwi_files)} ROIs to DWI space")
    else:
        dwi_files = list(roi_dwi_dir.glob("roi_*_dwi.mif"))
        print(f"\n--- Step 3: Skipped (found {len(dwi_files)} DWI-space ROIs) ---")
    
    # ===== STEP 4: Tractography =====
    if run_all or 'tractography' in steps:
        print(f"\n--- Step 4: Generate tractography ---")
        
        # Find available ROIs in DWI space
        dwi_rois = {f.stem.replace('roi_', '').replace('_dwi', ''): f 
                    for f in roi_dwi_dir.glob("roi_*_dwi.mif")}
        
        print(f"  Available ROIs: {list(dwi_rois.keys())}")
        
        # Define tract configurations based on available ROIs
        # Look for graft ROIs and their targets
        tract_configs = []
        
        for hemi in ['left', 'right']:
            graft_key = f"Graft_{hemi}"
            if graft_key not in dwi_rois:
                continue
            
            # Define targets for graft tractography
            targets = [
                (f"Striatum_{hemi}", 15),
                (f"Thalamus_{hemi}", 15),
                (f"Palecortex_{hemi}", 15),
                (f"cerebellum_{hemi}", 30),
            ]
            
            for target, maxlen in targets:
                if target in dwi_rois:
                    tract_configs.append((hemi, graft_key, target, maxlen))
        
        if not tract_configs:
            print("  ⚠ No valid graft→target pairs found")
            print("    Looking for alternative tractography options...")
            
            # Try striatum-cortex or other pairs
            for hemi in ['left', 'right']:
                striatum = f"Striatum_{hemi}"
                cortex = f"Palecortex_{hemi}"
                if striatum in dwi_rois and cortex in dwi_rois:
                    tract_configs.append((hemi, striatum, cortex, 15))
        
        if not tract_configs:
            print("  ❌ No tractography configurations possible with available ROIs")
            return True  # Not a failure, just no tractography possible
        
        print(f"\n  Running tractography for {len(tract_configs)} tract(s):")
        
        for hemi, seed_name, target_name, maxlen in tract_configs:
            seed_roi = dwi_rois[seed_name]
            target_roi = dwi_rois[target_name]
            
            # Clean names for output file
            seed_short = seed_name.replace(f"_{hemi}", "")
            target_short = target_name.replace(f"_{hemi}", "")
            
            tracks_raw = dwi_dir / f"tracks_{seed_short}_{target_short}_{hemi}.tck"
            tracks_sift = dwi_dir / f"tracks_{seed_short}_{target_short}_{hemi}_sift50.tck"
            
            print(f"\n    {seed_name} → {target_name}:")
            
            try:
                count = run_tractography(fod_file, seed_roi, target_roi, tracks_raw,
                                        select=5000, maxlength=maxlen)
                print(f"      Raw tracks: {count}")
                
                if count > 0:
                    run_sift(tracks_raw, fod_file, tracks_sift, term_number=50)
                    print(f"      SIFT50: ✓")
                else:
                    print(f"      SIFT: skipped (no tracks)")
            except Exception as e:
                print(f"      ❌ Failed: {e}")
        
        print(f"\n  ✓ Tractography complete")
    else:
        print(f"\n--- Step 4: Skipped ---")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Complete pipeline for manual ROI processing and tractography"
    )
    parser.add_argument("--subject", "-s", help="Process specific subject (e.g., sub-XXX)")
    parser.add_argument("--session", help="Process specific session")
    parser.add_argument("--step", action="append", choices=['convert', 'flip', 'transform', 'tractography', 'all'],
                       help="Run specific step(s). Default: all")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done without executing")
    parser.add_argument("--list", action="store_true", help="List subjects with NRRD files")
    
    args = parser.parse_args()
    
    base_dir = Path("/Users/glerma/soft/mbrain")
    deriv_dir = base_dir / "data" / "derivatives" / "mrtrix_preproc"
    
    # Find all subjects with NRRD files
    subjects = find_subjects_with_nrrd(deriv_dir)
    
    if args.list:
        print(f"\nSubjects with .seg.nrrd files:\n")
        for s in subjects:
            print(f"  {s['subject']}/{s['session']}")
            print(f"    NRRD: {s['nrrd_file'].name}")
        print(f"\nTotal: {len(subjects)} subjects")
        return
    
    # Filter by subject/session if specified
    if args.subject:
        subjects = [s for s in subjects if s['subject'] == args.subject]
        if args.session:
            subjects = [s for s in subjects if s['session'] == args.session]
    
    if not subjects:
        print("No subjects found with .seg.nrrd files")
        return
    
    print(f"\n{'#'*80}")
    print(f"# Manual ROI Pipeline")
    print(f"# Subjects to process: {len(subjects)}")
    print(f"{'#'*80}")
    
    steps = args.step if args.step else ['all']
    
    # Process each subject
    results = []
    for subj_info in subjects:
        success = process_subject(subj_info, steps=steps, dry_run=args.dry_run)
        results.append((subj_info['subject'], subj_info['session'], success))
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    for subject, session, success in results:
        status = "✓" if success else "❌"
        print(f"  {status} {subject}/{session}")
    
    successful = sum(1 for _, _, s in results if s)
    print(f"\nCompleted: {successful}/{len(results)} subjects")


if __name__ == "__main__":
    main()
