#!/usr/bin/env python3
"""
Main analysis script for mbrain - Mouse Brain fMRI and DTI Analysis

This script coordinates the BIDS conversion and subsequent MRtrix3 analysis.

Usage:
    python main.py --bids              # Only run BIDS conversion
    python main.py --preprocess        # Only run preprocessing
    python main.py --all               # Run both (default)
    python main.py --bids --dry-run    # Preview BIDS conversion
    python main.py --preprocess --steps denoise unring  # Run specific preprocessing steps
"""

import argparse
import os
from pathlib import Path
from utils import FlywheelToBIDSConverter, DWIPreprocessor, Config


def run_bids_conversion(dry_run=False):
    """Run BIDS conversion step."""
    print("\n" + "=" * 80)
    print("BIDS Conversion")
    print("=" * 80)
    print(f"Source: {Config.FLYWHEEL_DIR}")
    print(f"Target: {Config.BIDS_DIR}")
    print(f"Dry run: {dry_run}")
    print()
    
    # Initialize converter
    converter = FlywheelToBIDSConverter()
    
    # Run conversion
    stats = converter.convert(dry_run=dry_run)
    
    # Print summary
    print("\nConversion Summary:")
    print(f"  Subjects processed: {stats['subjects']}")
    print(f"  Sessions processed: {stats['sessions']}")
    print(f"  Acquisitions processed: {stats['acquisitions']}")
    print(f"  Files copied: {stats['files_copied']}")
    
    if stats['errors']:
        print(f"\n  Errors encountered: {len(stats['errors'])}")
        for error in stats['errors']:
            print(f"    - {error}")
    
    print(f"\nBIDS dataset: {Config.BIDS_DIR}")
    
    return stats


def run_preprocessing(steps=None, 
                      first_subject_only=False, 
                      dry_run=False):
    """Run DWI preprocessing step."""
    print("\n" + "=" * 80)
    print("DWI Preprocessing with MRtrix3")
    print("=" * 80)
    
    if dry_run:
        print("\n🔍 DRY RUN MODE - Commands will be displayed but not executed")
    if first_subject_only:
        print("\n⚠️  Running FIRST SUBJECT ONLY (testing mode)")

    # Orientation correction step: convert and flip all highres T2w and DWI files before anything else
    import subprocess
    from pathlib import Path
    print("\n=== Orientation Correction: mrconvert + mrtransform for all highres T2w and DWI files ===")
    bids_dir = Path(Config.BIDS_DIR)
    derivatives_dir = Path(Config.DERIVATIVES_DIR) / "mrtrix_preproc"

    # Process all highres T2w files
    anat_files = list(bids_dir.glob("sub-*/ses-*/anat/*_acq-highres_T2w.nii.gz"))
    if first_subject_only:
        anat_files = anat_files[:1]

    for anat_file in anat_files:
        parts = anat_file.parts
        sub = [p for p in parts if p.startswith("sub-")][0]
        ses = [p for p in parts if p.startswith("ses-")][0]
        out_dir = derivatives_dir / sub / ses / "anat"
        out_dir.mkdir(parents=True, exist_ok=True)
        mif1 = out_dir / "highres_T2w_perm_2_1_0.mif"
        base = anat_file.name.replace(".nii.gz", "").replace("_acq-highres_T2w", "")
        mif2 = out_dir / f"{base}_highres_T2w_flip.mif"
        # Step 1: mrconvert with axes permutation
        cmd1 = [
            "mrconvert", str(anat_file), str(mif1),
            "-axes", "2,1,0", "-strides", "1,2,3", "-force"
        ]
        print(" ".join(cmd1))
        subprocess.run(cmd1, check=True)
        # Step 2: mrtransform with axis flip
        cmd2 = [
            "mrtransform", str(mif1), str(mif2),
            "-flip", "1", "-force"
        ]
        print(" ".join(cmd2))
        subprocess.run(cmd2, check=True)
        if os.path.exists(mif1):
            os.remove(mif1)
        print(f"Converted and reoriented: {anat_file} -> {mif2}")

    # Process all DWI files (with bvec/bval)
    dwi_files = list(bids_dir.glob("sub-*/ses-*/dwi/*_dwi.nii.gz"))
    if first_subject_only:
        dwi_files = dwi_files[:1]

    for dwi_file in dwi_files:
        parts = dwi_file.parts
        sub = [p for p in parts if p.startswith("sub-")][0]
        ses = [p for p in parts if p.startswith("ses-")][0]
        out_dir = derivatives_dir / sub / ses / "dwi"
        out_dir.mkdir(parents=True, exist_ok=True)
        # Find corresponding bvec/bval
        base = dwi_file.name.replace(".nii.gz", "")
        bvec = dwi_file.parent / (base + ".bvec")
        bval = dwi_file.parent / (base + ".bval")
        mif1 = out_dir / "dwi_perm_2_1_0.mif"
        mif2 = out_dir /  (base + "_dwi_flip.mif")
        # Step 1: mrconvert with axes permutation and fslgrad (use -axes 2,1,0,3 for 4D DWI)
        cmd1 = [
            "mrconvert", str(dwi_file), str(mif1),
            "-fslgrad", str(bvec), str(bval),
            "-axes", "2,1,0,3", "-strides", "1,2,3,4", "-force"
        ]
        print(" ".join(cmd1))
        subprocess.run(cmd1, check=True)
        # Step 2: mrtransform with axis flip (add -reorient_fod no for multi-volume DWI)
        cmd2 = [
            "mrtransform", str(mif1), str(mif2),
            "-flip", "1", "-reorient_fod", "no", "-force"
        ]
        print(" ".join(cmd2))
        subprocess.run(cmd2, check=True)
        if os.path.exists(mif1):
            os.remove(mif1)
        print(f"Converted and reoriented: {dwi_file} -> {mif2}")

    if steps:
        print(f"\nRunning steps: {', '.join(steps)}")
    else:
        print("\nPreprocessing steps:")
        print("  0. Convert NIfTI to MIF with gradient information")
        print("  1. Gradient orientation check (dwigradcheck)")
        print("  2. Denoising (Marchenko-Pastur PCA)")
        print("  3. Gibbs ringing artifact removal")
        print("  4. Brain mask estimation")
        print("  5. Bias field correction (ANTs)")
        print("  6. Two-stage coregistration:")
        print("       • Stage 1: DWI → lowres T2w (matching resolution)")
        print("       • Stage 2: lowres → highres T2w (for ROI drawing)")
        print("       → Creates transform chain for ROI mapping: highres → lowres → DWI")
        print("  7. Atlas-based coregistration (Allen Mouse Brain CCF):")
        print("       • Step 7.1: N4 bias field correction on T2w")
        print("       • Step 7.2: Brain extraction using Allen template prior")
        print("       • Step 7.3: Registers Allen CCF template to subject space")
        print("       • Extracts and transforms cortex, striatum, thalamus ROIs")
        print("       → ROIs available in native DWI space for tractography")
    print()

    # Initialize preprocessor
    preprocessor = DWIPreprocessor()

    # Ensure gradcheck and all subsequent steps use dwi_flip.mif as the input
    # (Assumes DWIPreprocessor uses the correct file; if not, pass the path or update logic accordingly)
    # If you need to explicitly pass the file, you can add an argument or set an attribute here, e.g.:
    # preprocessor.input_dwi = str(mif2)  # where mif2 is dwi_flip.mif

    # Run preprocessing
    preproc_stats = preprocessor.preprocess_all_subjects(steps=steps,
                                derivatives_dir = derivatives_dir,
                                first_subject_only=first_subject_only, 
                                dry_run=dry_run, 
                                input_dwi_name="_dwi_flip.mif")

    # Print summary
    print("\nPreprocessing Summary:")
    print(f"  Subjects processed: {preproc_stats['subjects_processed']}")
    print(f"  Steps completed: {len(preproc_stats['steps_completed'])}")

    if preproc_stats['errors']:
        print(f"\n  Errors encountered: {len(preproc_stats['errors'])}")
        for error in preproc_stats['errors']:
            print(f"    - {error}")

    print(f"\nPreprocessed data: {Config.DERIVATIVES_DIR}/mrtrix_preproc")

    return preproc_stats


def run_tracking(steps=None, first_subject_only=False, dry_run=False):
    """Run tractography analysis."""
    print("\n" + "=" * 80)
    print("Tractography Analysis with MRtrix3")
    print("=" * 80)
    
    if dry_run:
        print("\n🔍 DRY RUN MODE - Commands will be displayed but not executed")
    if first_subject_only:
        print("\n⚠️  Running FIRST SUBJECT ONLY (testing mode)")
    
    if steps:
        print(f"\nRunning steps: {', '.join(steps)}")
    else:
        print("\nTracking steps:")
        print("  8. Response function estimation (dwi2response)")
        print("  9. Fiber orientation distribution (dwi2fod)")
        print("  10. Tractography (tckgen)")
    print()
    
    # Initialize preprocessor (handles both preprocessing and tracking)
    preprocessor = DWIPreprocessor()
    
    # Run tracking
    tracking_stats = preprocessor.preprocess_all_subjects(steps=steps, first_subject_only=first_subject_only, dry_run=dry_run)
    
    # Print summary
    print("\nTracking Summary:")
    print(f"  Subjects processed: {tracking_stats['subjects_processed']}")
    print(f"  Steps completed: {len(tracking_stats['steps_completed'])}")
    
    if tracking_stats['errors']:
        print(f"\n  Errors encountered: {len(tracking_stats['errors'])}")
        for error in tracking_stats['errors']:
            print(f"    - {error}")
    
    print(f"\nTracking results: {Config.DERIVATIVES_DIR}/mrtrix_preproc")
    
    return tracking_stats


def main():
    """Main analysis pipeline."""
    parser = argparse.ArgumentParser(
        description='mbrain - Mouse Brain fMRI and DTI Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --bids                          # Only convert to BIDS
  %(prog)s --preprocess                    # Only run preprocessing
  %(prog)s --tracking                      # Only run tractography
  %(prog)s --all                           # Run all steps (BIDS + preprocessing + tracking)
  %(prog)s --bids --dry-run                # Preview BIDS conversion
  %(prog)s --preprocess --steps denoise unring  # Run specific preprocessing steps
  %(prog)s --tracking --steps dwi2fod tckgen    # Run specific tracking steps
        """
    )
    
    # Main operation modes
    parser.add_argument('--bids', action='store_true',
                        help='Run BIDS conversion')
    parser.add_argument('--preprocess', action='store_true',
                        help='Run DWI preprocessing')
    parser.add_argument('--tracking', action='store_true',
                        help='Run tractography analysis')
    parser.add_argument('--all', action='store_true',
                        help='Run all steps (BIDS conversion + preprocessing + tracking)')
    
    # Options
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview operations without executing (show commands for preprocessing, skip file copy for BIDS)')
    
    # Preprocessing options
    parser.add_argument('--steps', nargs='+', 
                        choices=['gradcheck', 'denoise', 'unring', 'mask', 'bias_correct', 'coreg', 'atlas_coreg', 'dwi2response', 'dwi2fod', 'tckgen'],
                        help='Specific preprocessing or tracking steps to run. '
                             'Preprocessing: gradcheck, denoise, unring, mask, bias_correct, coreg, atlas_coreg. '
                             'Tracking: dwi2response, dwi2fod, tckgen')
    parser.add_argument('--first-subject', action='store_true',
                        help='Process only the first subject (for testing)')
    
    args = parser.parse_args()
    
    # Default to --all if no mode specified
    if not (args.bids or args.preprocess or args.tracking or args.all):
        args.all = True
    
    # Validate steps based on mode
    preprocess_steps = ['gradcheck', 'denoise', 'unring', 'mask', 'bias_correct', 'coreg', 'atlas_coreg']
    tracking_steps = ['dwi2response', 'dwi2fod', 'tckgen']
    
    if args.steps:
        if args.preprocess and not args.tracking:
            # Check that all steps are preprocessing steps
            invalid = [s for s in args.steps if s not in preprocess_steps]
            if invalid:
                parser.error(f"Invalid preprocessing steps: {invalid}. Use --tracking for: {tracking_steps}")
        elif args.tracking and not args.preprocess:
            # Check that all steps are tracking steps
            invalid = [s for s in args.steps if s not in tracking_steps]
            if invalid:
                parser.error(f"Invalid tracking steps: {invalid}. Use --preprocess for: {preprocess_steps}")
    
    print("=" * 80)
    print("mbrain - Mouse Brain Analysis Pipeline")
    print("=" * 80)
    
    # Run requested steps
    if args.all or args.bids:
        run_bids_conversion(dry_run=args.dry_run)
    
    if args.all or args.preprocess:
        # Filter steps for preprocessing
        preprocess_steps_to_run = None
        if args.steps and not args.all:
            preprocess_steps_to_run = [s for s in args.steps if s in preprocess_steps]
            if not preprocess_steps_to_run:
                preprocess_steps_to_run = None  # Run all preprocessing steps
        run_preprocessing(steps=preprocess_steps_to_run, first_subject_only=args.first_subject, dry_run=args.dry_run)
    
    if args.all or args.tracking:
        # Filter steps for tracking
        tracking_steps_to_run = None
        if args.steps and not args.all:
            tracking_steps_to_run = [s for s in args.steps if s in tracking_steps]
            if not tracking_steps_to_run:
                tracking_steps_to_run = None  # Run all tracking steps
        run_tracking(steps=tracking_steps_to_run, first_subject_only=args.first_subject, dry_run=args.dry_run)
    
    # Print next steps
    print("\n" + "=" * 80)
    print("Pipeline Complete!")
    print("=" * 80)
    print("\nNext steps:")
    print("  1. Validate preprocessed data with mrview")
    print("  2. Validate tractography results with mrview")
    print("  3. Apply SIFT filtering (tcksift/tcksift2)")
    print("  4. Generate tract density images (tckmap)")
    print("  5. Extract tract metrics and build connectomes")


if __name__ == "__main__":
    main()

'''
dwidenoise -extent 5,5,5 -noise noise.mif -estimator Exp2 \
/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed.nii.gz \
/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised.nii.gz \
-force

dwi2mask /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised.nii.gz \
/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised_mask.nii.gz \
-force





cmd = [ "tckgen -quiet ", "-algorithm", ts.algorithm, ...
                        "-select ", ts.select, ...
                        "-seed_image", roi1, "-seed_unidirectional", ...
                        "-include", roi2, ...
                        roi3, roi4, roiexc1, roiexc2, ...
                        "-angle", ts.tckangle, "-cutoff", ts.cutoff, ...
                        "-minlength", ts.tckminlen, "-maxlength", ts.tckmaxlen, ...
                        "-stop", ...
                        input_file, ...
                        ts.fpath];
    
cmd_str = ['tckgen -force ' ...
                      '-algo ' algo  ' ' optionalStr ' ' ...
                      '-backtrack -crop_at_gmwmi -info ' ...
                      '-seed_gmwmi ' files.gmwmi ' ' ...
                      '-act ' files.tt5 ' ' ...
                      '-angle ' num2str(ET_angleValues(na)) ' ' ...
                      '-select ' num2str(ET_numberFibers) ' ' ...
                      '-maxlength ' num2str(ET_maxlength(na)) ' ' ...
                      input_file ' ' ...
                      fgFileNameWithDir{na}];

cmd_str = ['tckgen -force ' ...
                    '-algo ' algo  ' ' ...
                    optionalStr ...
                    '-seed_image ' roi ' ' ...
                    '-mask ' mask ' ' ...
                    '-angle ' num2str(ET_angleValues(na)) ' ' ...
                    '-select ' num2str(ET_numberFibers) ' ' ...
                    '-maxlength ' num2str(ET_maxlength(na)) ' ' ...
                    input_file ' ' ...
                    fgFileNameWithDir_presift{na}];
cmd_sift = ['tcksift -force ' ...
                    '-term_number ' num2str(sift_nFibers) ' ' ...
                    fgFileNameWithDir_presift{na} ' ' ...
                    input_file ' ' ...
                    fgFileNameWithDir{na}];    


tckgen \
  '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/' \
'ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_fod_wm.mif' \
  '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/' \
'ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_tracks.tck' \
  -algorithm iFOD2 \
  -seed_image \
    '/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/' \
    'ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_mask.mif' \
  -select 10000000 \
  -step 0.02 \
  -angle 60 \
  -minlength 10 \
  -maxlength 50 \
  -cutoff 0.05 \
  -force


  antsRegistration --verbose 1 --dimensionality 3 --float 0 \
                   --collapse-output-transforms 1 --output \
                    [ /Users/glerma/soft/mbx

antsRegistrationSyN.sh \
  -d 3 \
  -t b \ glerma@MBPM2:~/soft/mbrain:$ 
  -f \RegistrationSyN.sh \
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/allen_template_\
  -m \ \
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/sub-596-C-OO-Nm\
  -o \sers/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/allen_template_\
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/sub-596-C-OO-Nm_
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/atlas/sub-596-C-OO-Nm\

mrdegibbs \
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised.nii.gz \
    /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised_unringed.nii.gz \
    -force

dwi2mask \
/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised_unringed.nii.gz \
/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-596-C-OO-NmM/ses-240911103349/dwi/sub-596-C-OO-NmM_ses-240911103349_dwi_preprocessed_denoised_unringed.nii.gz
                '-force'

                



ANAT
        cmd1 = [
            "mrconvert", str(anat_file), str(mif1),
            "-axes", "2,1,0", "-strides", "1,2,3", "-force"
        ]
        print(" ".join(cmd1))
        subprocess.run(cmd1, check=True)
        # Step 2: mrtransform with axis flip
        cmd2 = [
            "mrtransform", str(mif1), str(mif2),
            "-flip", "1", "-force"



DWI
                
# 0. Convert original DWI to .mif (if starting from NIfTI)
mrconvert {input_dwi_nii} {basename}_raw.mif -fslgrad {bvecs} {bvals} -force

# 0b. (Optional) Flip axes if needed (example: flip x axis)
mrconvert {basename}_raw.mif {basename}_dwi_flip.mif -coord 0 -1:0 -force

# 1. Gradient orientation check
dwigradcheck {basename}_dwi_flip.mif -export_grad_mrtrix raw.b -force
mrconvert {basename}_dwi_flip.mif -grad raw.b {basename}_gradchecked.mif -force

# 2. Denoising
dwidenoise {basename}_gradchecked.mif {basename}_denoised.mif -noise {basename}_noise.mif -force

# 3. Remove Gibbs ringing artifacts
mrdegibbs {basename}_denoised.mif {basename}_unringed.mif -force

# 4. Brain mask estimation
dwi2mask {basename}_unringed.mif {basename}_mask.mif -force

# 5. Bias field correction
dwibiascorrect ants {basename}_unringed.mif {basename}_preprocessed.mif -mask {basename}_mask.mif -force

# 6. Convert to NIfTI
mrconvert {basename}_preprocessed.mif {basename}_preprocessed.nii.gz -force


ANAT
        cmd1 = [
            "mrconvert", str(anat_file), str(mif1),
            "-axes", "2,1,0", "-strides", "1,2,3", "-force"
        ]
        print(" ".join(cmd1))
        subprocess.run(cmd1, check=True)
        # Step 2: mrtransform with axis flip
        cmd2 = [
            "mrtransform", str(mif1), str(mif2),
            "-flip", "1", "-force"



DWI
                
# 0. Convert original DWI to .mif (if starting from NIfTI)
mrconvert {input_dwi_nii} {basename}_raw.mif -fslgrad {bvecs} {bvals} -force

# 0b. (Optional) Flip axes if needed (example: flip x axis)
mrconvert {basename}_raw.mif {basename}_dwi_flip.mif -coord 0 -1:0 -force

# 1. Gradient orientation check
dwigradcheck {basename}_dwi_flip.mif -export_grad_mrtrix raw.b -force
mrconvert {basename}_dwi_flip.mif -grad raw.b {basename}_gradchecked.mif -force

# 2. Denoising
dwidenoise {basename}_gradchecked.mif {basename}_denoised.mif -noise {basename}_noise.mif -force

# 3. Remove Gibbs ringing artifacts
mrdegibbs {basename}_denoised.mif {basename}_unringed.mif -force

# 4. Brain mask estimation
dwi2mask {basename}_unringed.mif {basename}_mask.mif -force

# 5. Bias field correction
dwibiascorrect ants {basename}_unringed.mif {basename}_preprocessed.mif -mask {basename}_mask.mif -force

# 6. Convert to NIfTI
mrconvert {basename}_preprocessed.mif {basename}_preprocessed.nii.gz -force

GROUP AP: [sub-596-C-OO-NmM, sub-769-C-O5-R1M, sub-769-C-OO-NmM, sub-875-C-OO-NM]

GROUP APTR: [sub-AP72-862-901-NmF, sub-AP72-862-901-R1M, sub-Apallial17NmF, sub-Apallial20R1F]

GROUP CTRL: [sub-129-921-R1F,sub-129-921-R1L1F,sub-129-1014-L1M,sub-138-1010-L1F]


(1) Isocortex
(2) Piriform cortex + Entorhinal area + Piriform amygdalar area + Postpiriform transition area
(3) Striatum dorsal region + Striatum ventral region
(4) Thalamus
(5) Cerebellum

cd /Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/sub-129-921-R1F/ses-260115121231/dwi && \
mrview t2w_highres_flip.mif \
-overlay.load atlas/rois/roi_bla_left.nii.gz -overlay.opacity 0.5 \
-overlay.load atlas/rois/roi_bla_right.nii.gz -overlay.opacity 0.5 \
-overlay.load atlas/rois/roi_aca_left.nii.gz -overlay.opacity 0.5 \
-overlay.load atlas/rois/roi_aca_right.nii.gz -overlay.opacity 0.5 \
-tractography.load tracks_BLA_ACA_left_top100.tck \
-tractography.load tracks_BLA_ACA_right_top100.tck &

'''


