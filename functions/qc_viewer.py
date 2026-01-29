#!/usr/bin/env python3
"""
QC Viewer - View key outputs for all processed subjects.

Usage:
    python qc_viewer.py                    # View all subjects one by one
    python qc_viewer.py --subject sub-596  # View specific subject
    python qc_viewer.py --list             # List status of all subjects
    python qc_viewer.py --all-at-once      # Open all subjects simultaneously
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Dict
from config_loader import load_config


class QCViewer:
    """QC viewer for processed subjects."""
    
    def __init__(self):
        config = load_config()
        self.base_dir = Path(config['paths']['repo_dir'])
        self.deriv_dir = Path(config['paths']['derivatives_dir'])
    
    def get_subjects(self) -> List[str]:
        """Get list of all subjects in derivatives."""
        subjects = []
        for item in self.deriv_dir.iterdir():
            if item.is_dir() and item.name.startswith("sub-"):
                subjects.append(item.name)
        return sorted(subjects)
    
    def get_session(self, subject: str) -> Optional[str]:
        """Get first session for a subject."""
        subj_dir = self.deriv_dir / subject
        for item in subj_dir.iterdir():
            if item.is_dir() and item.name.startswith("ses-"):
                return item.name
        return None
    
    def get_subject_files(self, subject: str) -> Dict[str, Path]:
        """Get paths to key files for a subject."""
        session = self.get_session(subject)
        if not session:
            return {}
        
        out_dir = self.deriv_dir / subject / session / "dwi"
        
        files = {
            "t2w": out_dir / "t2w_highres_flip.nii.gz",
            "isocortex": out_dir / "atlas" / "rois" / "roi_isocortex.nii.gz",
            "thalamus": out_dir / "atlas" / "rois" / "roi_thalamus.nii.gz",
            "fa": out_dir / "fa.mif",
            "mask": out_dir / "mask.mif",
            "dwi": out_dir / "dwi_flip_gradcheck_denoise_bias_anat.mif",
        }
        
        # Find tracks file (try different naming conventions)
        track_names = [
            "tracks_thalamus_cortex_siftratio0.1.tck",
            "tracks_thalamus_cortex_sift50.tck",
            "tracks_thalamus_cortex_siftratio1.tck",
            "tracks_thalamus_cortex.tck",
        ]
        
        for track_name in track_names:
            track_path = out_dir / track_name
            if track_path.exists():
                files["tracks"] = track_path
                break
        
        return files
    
    def check_subject_status(self, subject: str) -> Dict[str, bool]:
        """Check which key files exist for a subject."""
        files = self.get_subject_files(subject)
        
        status = {}
        for key, path in files.items():
            status[key] = path.exists() if path else False
        
        return status
    
    def get_fa_stats(self, subject: str) -> Optional[Dict[str, float]]:
        """Get FA statistics for a subject."""
        files = self.get_subject_files(subject)
        
        fa_path = files.get("fa")
        mask_path = files.get("mask")
        
        if not fa_path or not fa_path.exists():
            return None
        if not mask_path or not mask_path.exists():
            return None
        
        try:
            result = subprocess.run(
                f"mrstats {fa_path} -mask {mask_path} -output mean,std,min,max",
                shell=True,
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                values = result.stdout.strip().split()
                if len(values) >= 4:
                    return {
                        "mean": float(values[0]),
                        "std": float(values[1]),
                        "min": float(values[2]),
                        "max": float(values[3]),
                    }
        except:
            pass
        
        return None
    
    def list_all_subjects(self):
        """Print status of all subjects."""
        subjects = self.get_subjects()
        
        print("=" * 60)
        print("QC STATUS - ALL SUBJECTS")
        print("=" * 60)
        
        for subject in subjects:
            session = self.get_session(subject)
            status = self.check_subject_status(subject)
            
            print(f"\n--- {subject} / {session or 'NO SESSION'} ---")
            
            if not session:
                print("  ✗ No session found")
                continue
            
            # Check key files
            t2w_ok = "✓" if status.get("t2w") else "✗"
            dwi_ok = "✓" if status.get("dwi") else "✗"
            fa_ok = "✓" if status.get("fa") else "✗"
            rois_ok = "✓" if (status.get("isocortex") and status.get("thalamus")) else "✗"
            tracks_ok = "✓" if status.get("tracks") else "✗"
            
            print(f"  T2w processed:  {t2w_ok}")
            print(f"  DWI processed:  {dwi_ok}")
            print(f"  FA computed:    {fa_ok}")
            print(f"  ROIs warped:    {rois_ok}")
            print(f"  Tractography:   {tracks_ok}")
            
            # FA stats
            fa_stats = self.get_fa_stats(subject)
            if fa_stats:
                print(f"  Mean FA:        {fa_stats['mean']:.4f} ± {fa_stats['std']:.4f}")
        
        print("\n" + "=" * 60)
    
    def view_subject(self, subject: str, wait: bool = True) -> bool:
        """Open mrview for a subject."""
        files = self.get_subject_files(subject)
        status = self.check_subject_status(subject)
        
        # Check required files
        required = ["t2w", "isocortex", "thalamus", "tracks"]
        missing = [f for f in required if not status.get(f)]
        
        if missing:
            print(f"  ✗ Missing: {', '.join(missing)}")
            return False
        
        print(f"  ✓ All files present")
        print(f"  Opening mrview...")
        
        # Build mrview command
        cmd = [
            "mrview", str(files["t2w"]),
            "-overlay.load", str(files["isocortex"]),
            "-overlay.load", str(files["thalamus"]),
            "-tractography.load", str(files["tracks"]),
        ]
        
        if wait:
            # Run and wait for it to close
            subprocess.run(cmd)
        else:
            # Run in background
            subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        return True
    
    def view_all_subjects(self, wait_between: bool = True):
        """View all subjects one by one."""
        subjects = self.get_subjects()
        
        print("=" * 50)
        print("QC VIEWER - Thalamus-Cortex Tractography")
        print("=" * 50)
        
        for i, subject in enumerate(subjects):
            session = self.get_session(subject)
            
            print(f"\n--- [{i+1}/{len(subjects)}] {subject} / {session} ---")
            
            if not session:
                print("  ✗ No session found")
                continue
            
            success = self.view_subject(subject, wait=wait_between)
            
            if wait_between and success and i < len(subjects) - 1:
                try:
                    input("\nPress ENTER for next subject (Ctrl+C to quit)...")
                except KeyboardInterrupt:
                    print("\n\nStopped by user.")
                    break
        
        print("\nDone!")
    
    def view_all_at_once(self):
        """Open all subjects simultaneously."""
        subjects = self.get_subjects()
        
        print("Opening all subjects...")
        
        opened = 0
        for subject in subjects:
            session = self.get_session(subject)
            if not session:
                continue
            
            print(f"  {subject}...", end=" ")
            if self.view_subject(subject, wait=False):
                opened += 1
                print("opened")
            else:
                print("skipped")
        
        print(f"\nOpened {opened} subjects.")


def main():
    parser = argparse.ArgumentParser(
        description="QC Viewer for processed subjects",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python qc_viewer.py --list              # List status of all subjects
  python qc_viewer.py                     # View all subjects one by one
  python qc_viewer.py --subject sub-596   # View specific subject
  python qc_viewer.py --all-at-once       # Open all subjects simultaneously
        """
    )
    
    parser.add_argument(
        "--subject", "-s",
        help="View specific subject (partial match supported)"
    )
    parser.add_argument(
        "--list", "-l",
        action="store_true",
        help="List status of all subjects (no viewer)"
    )
    parser.add_argument(
        "--all-at-once", "-a",
        action="store_true",
        help="Open all subjects simultaneously"
    )
    parser.add_argument(
        "--base-dir", "-b",
        default="/Users/glerma/soft/mbrain",
        help="Base directory for the project"
    )
    
    args = parser.parse_args()
    
    viewer = QCViewer(base_dir=args.base_dir)
    
    if args.list:
        viewer.list_all_subjects()
    elif args.subject:
        # Find matching subject
        subjects = viewer.get_subjects()
        matches = [s for s in subjects if args.subject in s]
        
        if not matches:
            print(f"No subjects matching '{args.subject}'")
            sys.exit(1)
        elif len(matches) > 1:
            print(f"Multiple matches: {matches}")
            print("Please be more specific.")
            sys.exit(1)
        else:
            subject = matches[0]
            session = viewer.get_session(subject)
            print(f"--- {subject} / {session} ---")
            viewer.view_subject(subject, wait=True)
    elif args.all_at_once:
        viewer.view_all_at_once()
    else:
        viewer.view_all_subjects()


if __name__ == "__main__":
    main()
