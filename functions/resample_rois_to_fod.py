#!/usr/bin/env python3
"""
Resample ROIs and segments from T2w space to FOD space for tractography (flip1 pipeline).

For AP/APTR subjects (manual ROIs):
    - Resamples manual ROIs from manual_roi/ to manual_roi/fod_space/
    - Resamples graft segments from dwi/graft_segments_t2w/ to dwi/graft_segments_fod/

For CTRL subjects (atlas ROIs):
    - Resamples atlas ROIs to fod_space
    - Resamples isocortex segments to fod_space

Usage:
    python resample_rois_to_fod.py --subject sub-Apallial17NmF
    python resample_rois_to_fod.py --group APTR
    python resample_rois_to_fod.py --group CTRL --force
"""

import argparse
import subprocess
import os
from pathlib import Path
from typing import Optional, Tuple, List

# Force flip1 config
os.environ['MBRAIN_PIPELINE'] = 'flip1'
from config_loader import load_config


def get_subjects(cfg, group=None, subject=None):
    """Get list of (subject_id, session, group) tuples based on filters."""
    groups = cfg.get('groups', {})
    subjects = []
    
    if subject:
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session'], g))
                    return subjects
        raise ValueError(f"Subject {subject} not found in any group")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found. Available: {list(groups.keys())}")
        for s in groups[group]:
            subjects.append((s['id'], s['session'], group))
    
    return subjects


def run_cmd(cmd, description=None):
    """Run shell command with optional description."""
    if description:
        print(f"    {description}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[:500]}")
        return False
    return True


def resample_roi(input_file, output_file, template_file):
    """Resample a single ROI to match template space using nearest neighbor."""
    cmd = f"mrgrid '{input_file}' regrid '{output_file}' -template '{template_file}' -interp nearest -force -quiet"
    return run_cmd(cmd)


def get_roi_voxel_count(roi_file):
    """Get count of non-zero voxels in ROI."""
    result = subprocess.run(
        f"mrstats '{roi_file}' -mask '{roi_file}' -ignorezero -output count",
        shell=True, capture_output=True, text=True
    )
    if result.returncode == 0:
        try:
            return int(float(result.stdout.strip()))
        except:
            return 0
    return 0


def process_manual_rois(subject, session, deriv, template, force):
    """Resample manual ROIs for AP/APTR subjects."""
    
    roi_in_dir = deriv / subject / session / "manual_roi"
    roi_out_dir = roi_in_dir / "fod_space"
    
    roi_files = list(roi_in_dir.glob("roi_*_flip.mif"))
    if not roi_files:
        print(f"    WARNING: No manual ROIs found in {roi_in_dir}")
        return False
    
    roi_out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n  Manual ROIs ({len(roi_files)} files):")
    success_count = 0
    for roi_file in sorted(roi_files):
        out_name = roi_file.stem + "_fod.mif"
        out_file = roi_out_dir / out_name
        
        if out_file.exists() and not force:
            count = get_roi_voxel_count(out_file)
            print(f"    + {roi_file.name} -> {out_name} (exists, {count} voxels)")
            success_count += 1
        elif resample_roi(roi_file, out_file, template):
            count = get_roi_voxel_count(out_file)
            print(f"    + {roi_file.name} -> {out_name} ({count} voxels)")
            success_count += 1
        else:
            print(f"    X {roi_file.name} FAILED")
    
    return success_count == len(roi_files)


def process_graft_segments(subject, session, deriv, template, force):
    """Resample graft segments for AP/APTR subjects."""
    
    seg_in_dir = deriv / subject / session / "dwi" / "graft_segments_t2w"
    seg_out_dir = deriv / subject / session / "dwi" / "graft_segments_fod"
    
    if not seg_in_dir.exists():
        print(f"    WARNING: No graft segments found in {seg_in_dir}")
        return True
    
    seg_files = list(seg_in_dir.glob("seg*.nii.gz"))
    if not seg_files:
        print(f"    WARNING: No segment files found in {seg_in_dir}")
        return True
    
    seg_out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n  Graft Segments ({len(seg_files)} files):")
    success_count = 0
    for seg_file in sorted(seg_files):
        out_name = seg_file.stem.replace('.nii', '') + "_fod.mif"
        out_file = seg_out_dir / out_name
        
        if out_file.exists() and not force:
            print(f"    + {seg_file.name} (exists)")
            success_count += 1
        elif resample_roi(seg_file, out_file, template):
            print(f"    + {seg_file.name}")
            success_count += 1
        else:
            print(f"    X {seg_file.name} FAILED")
    
    return success_count == len(seg_files)


def process_atlas_rois(subject, session, deriv, template, force):
    """Resample atlas ROIs for CTRL subjects."""
    
    atlas_dir = deriv / subject / session / "dwi" / "atlas" / "rois"
    out_dir = atlas_dir / "fod_space"
    
    if not atlas_dir.exists():
        print(f"    WARNING: No atlas ROIs found in {atlas_dir}")
        return False
    
    roi_files = [f for f in atlas_dir.glob("roi_*.nii.gz") if "fod_space" not in str(f)]
    if not roi_files:
        print(f"    WARNING: No ROI files found in {atlas_dir}")
        return False
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n  Atlas ROIs ({len(roi_files)} files):")
    success_count = 0
    for roi_file in sorted(roi_files):
        out_file = out_dir / roi_file.name
        
        if out_file.exists() and not force:
            print(f"    + {roi_file.name} (exists)")
            success_count += 1
        elif resample_roi(roi_file, out_file, template):
            print(f"    + {roi_file.name}")
            success_count += 1
        else:
            print(f"    X {roi_file.name} FAILED")
    
    return success_count > 0


def process_isocortex_segments(subject, session, deriv, template, force):
    """Resample isocortex segments for CTRL subjects."""
    
    seg_dir = deriv / subject / session / "dwi" / "isocortex_segments_v2"
    out_dir = seg_dir / "fod_space"
    
    if not seg_dir.exists():
        print(f"    WARNING: No isocortex segments found in {seg_dir}")
        return True
    
    seg_files = [f for f in seg_dir.glob("isocortex_*.nii.gz") if "fod_space" not in str(f)]
    if not seg_files:
        print(f"    WARNING: No segment files found in {seg_dir}")
        return True
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n  Isocortex Segments ({len(seg_files)} files):")
    success_count = 0
    for seg_file in sorted(seg_files):
        out_file = out_dir / seg_file.name
        
        if out_file.exists() and not force:
            print(f"    + {seg_file.name} (exists)")
            success_count += 1
        elif resample_roi(seg_file, out_file, template):
            print(f"    + {seg_file.name}")
            success_count += 1
        else:
            print(f"    X {seg_file.name} FAILED")
    
    return success_count == len(seg_files)


def process_subject(subject, session, group, cfg, force=False):
    """Resample ROIs to FOD space for a single subject."""
    paths = cfg.get('paths', {})
    deriv = Path(paths['derivatives_dir'])
    dwi_dir = deriv / subject / session / "dwi"
    
    print(f"\n{'='*60}")
    print(f"Processing: {subject}/{session} (Group: {group})")
    print(f"{'='*60}")
    
    template = dwi_dir / "wm_fod_norm.mif"
    if not template.exists():
        print(f"  ERROR: FOD template not found: {template}")
        print(f"  Run dwi_fod.py first to generate FOD")
        return False
    
    success = True
    
    if group in ['AP', 'APTR']:
        if not process_manual_rois(subject, session, deriv, template, force):
            success = False
        if not process_graft_segments(subject, session, deriv, template, force):
            success = False
    else:
        if not process_atlas_rois(subject, session, deriv, template, force):
            success = False
        if not process_isocortex_segments(subject, session, deriv, template, force):
            success = False
    
    if success:
        print(f"\n  Done!")
    else:
        print(f"\n  Some resampling failed")
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description='Resample ROIs and segments to FOD space (flip1 pipeline)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python resample_rois_to_fod.py --subject sub-Apallial17NmF
    python resample_rois_to_fod.py --group APTR
    python resample_rois_to_fod.py --group CTRL --force

Output directories:
    AP/APTR: manual_roi/fod_space/, dwi/graft_segments_fod/
    CTRL:    dwi/atlas/rois/fod_space/, dwi/isocortex_segments_v2/fod_space/
        """
    )
    parser.add_argument('--group', type=str,
                        help='Process all subjects in group (CTRL, AP, APTR)')
    parser.add_argument('--subject', type=str,
                        help='Process specific subject')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration even if output exists')
    
    args = parser.parse_args()
    
    cfg = load_config()
    
    print(f"Pipeline: flip1")
    print(f"Derivatives: {cfg['paths']['derivatives_dir']}")
    
    if not args.group and not args.subject:
        parser.error("Must specify --group or --subject")
    
    subjects = get_subjects(cfg, group=args.group, subject=args.subject)
    
    print(f"Processing {len(subjects)} subjects...")
    
    results = []
    for subject_id, session, group in subjects:
        success = process_subject(subject_id, session, group, cfg, force=args.force)
        results.append((subject_id, success))
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for subj, success in results:
        status = "+" if success else "X"
        print(f"  {status} {subj}")


if __name__ == "__main__":
    main()
