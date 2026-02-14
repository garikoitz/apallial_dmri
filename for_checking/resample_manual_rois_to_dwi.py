#!/usr/bin/env python3
"""
Resample manual ROIs from T2w highres space to DWI space.

Uses nearest-neighbor interpolation to preserve discrete label values.
No interpolation - labels remain exactly 0, 1, 2, etc.
"""

import subprocess
import argparse
from pathlib import Path


def resample_roi_to_dwi(roi_path, dwi_template, output_path):
    """
    Resample ROI to DWI space using nearest-neighbor interpolation.
    
    Args:
        roi_path: Path to ROI NIfTI file
        dwi_template: Path to DWI file to use as template geometry
        output_path: Path for output resampled ROI
    """
    cmd = [
        "mrtransform",
        str(roi_path),
        str(output_path),
        "-template", str(dwi_template),
        "-interp", "nearest",  # Preserve discrete label values
        "-force"
    ]
    
    print(f"Resampling: {roi_path.name}")
    print(f"  Template: {dwi_template.name}")
    print(f"  Output: {output_path.name}")
    
    subprocess.run(cmd, check=True)
    
    # Verify label values are preserved
    result = subprocess.run(
        ["mrstats", str(output_path), "-output", "min", "-output", "max"],
        capture_output=True, text=True
    )
    if result.returncode == 0:
        values = result.stdout.strip().split()
        print(f"  ✓ Label range: {values[0]} to {values[1]}")


def main():
    parser = argparse.ArgumentParser(
        description="Resample manual ROIs to DWI space with nearest-neighbor interpolation"
    )
    parser.add_argument("--subject", "-s", default="sub-Apallial17NmF",
                       help="Subject ID")
    parser.add_argument("--session", default="ses-240816091633",
                       help="Session ID")
    parser.add_argument("--visualize", action="store_true",
                       help="Open mrview after resampling")
    
    args = parser.parse_args()
    
    # Paths
    base_dir = Path(f"/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc/{args.subject}/{args.session}")
    manual_roi_dir = base_dir / "manual_roi"
    dwi_dir = base_dir / "dwi"
    
    # Input ROIs
    roi_files = list(manual_roi_dir.glob("roi_*.nii.gz"))
    labelmap = manual_roi_dir / "Segmentation_6_labels.nii.gz"
    
    # DWI template (use FA as template - same geometry as FOD)
    dwi_template = dwi_dir / "fa.mif"
    
    if not dwi_template.exists():
        print(f"❌ DWI template not found: {dwi_template}")
        return
    
    print("="*60)
    print("Resampling Manual ROIs to DWI Space")
    print("="*60)
    print(f"Subject: {args.subject}")
    print(f"Session: {args.session}")
    print(f"Template: {dwi_template.name}")
    print()
    
    # Resample individual ROIs
    resampled_files = []
    for roi_file in roi_files:
        output_file = dwi_dir / f"{roi_file.stem}_dwi.nii.gz"
        resample_roi_to_dwi(roi_file, dwi_template, output_file)
        resampled_files.append(output_file)
        print()
    
    # Resample combined labelmap
    if labelmap.exists():
        output_labelmap = dwi_dir / "manual_rois_dwi_labels.nii.gz"
        resample_roi_to_dwi(labelmap, dwi_template, output_labelmap)
        resampled_files.append(output_labelmap)
        print()
    
    print("="*60)
    print("✓ Resampling Complete")
    print("="*60)
    print("\nResampled ROIs (in DWI space):")
    for f in resampled_files:
        print(f"  {f.name}")
    
    # Visualize if requested
    if args.visualize:
        print("\nLaunching mrview...")
        cmd = [
            "mrview",
            str(dwi_template),
            "-mode", "2"
        ]
        
        # Add individual ROIs as overlays
        for roi_file in resampled_files:
            if "labels" not in roi_file.name:  # Skip combined labelmap
                cmd.extend([
                    "-overlay.load", str(roi_file),
                    "-overlay.opacity", "0.5"
                ])
        
        subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
