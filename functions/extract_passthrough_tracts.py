#!/usr/bin/env python3
"""
Extract passthrough streamlines from whole-brain tractograms.

This script extracts streamlines that pass through BOTH a seed segment ROI 
AND a target ROI from whole-brain tractograms. Uses tckedit with -include flags.

For CTRL subjects:
    - Seed: 5 isocortex segments (anterior, anterior_central, central, posterior_central, posterior)
    - Target: thalamus, striatum, cerebellum, olfactory_cortex

For AP/APTR subjects:
    - Seed: 3 dorsal graft segments (anterior_dorsal, central_dorsal, posterior_dorsal)
    - Modes: 6seg and 9seg
    - Target: Thalamus, Striatum, Palecortex, cerebellum

Output: dwi/wholebrain/passthrough_tracts/
        dwi/wholebrain/passthrough_tracts_shrink/ (when --shrink is used)

Usage:
    python functions/extract_passthrough_tracts.py                    # All subjects, left hemisphere
    python functions/extract_passthrough_tracts.py --hemi right       # All subjects, right hemisphere
    python functions/extract_passthrough_tracts.py --group CTRL       # CTRL only
    python functions/extract_passthrough_tracts.py --subject sub-129-921-R1F
    python functions/extract_passthrough_tracts.py --shrink           # Use shrunk segments (AP/APTR only)
    python functions/extract_passthrough_tracts.py --create-shrink-segments  # Create shrunk segments
    python functions/extract_passthrough_tracts.py --dry-run          # Preview
"""

import argparse
import subprocess
import os
import sys
import tempfile
from pathlib import Path
from datetime import datetime

# Set pipeline before importing config_loader
os.environ['MBRAIN_PIPELINE'] = 'flip1'

from config_loader import load_config, get_paths, get_subject_groups


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session, group) tuples based on filters."""
    groups = get_subject_groups()
    subjects = []
    
    if subject:
        for g, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session'], g))
                    return subjects
        raise ValueError(f"Subject {subject} not found")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found. Available: {list(groups.keys())}")
        for s in groups[group]:
            subjects.append((s['id'], s['session'], group))
    else:
        for g, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session'], g))
    
    return subjects


def create_shrunk_segments(deriv_dir):
    """
    Create shrunk (eroded by 1 voxel) versions of ALL graft segments for AP/APTR.
    Saves to dwi/graft_segments_fod/shrink/ folder.
    
    Uses scipy.ndimage.binary_erosion with a 3D ball structuring element.
    """
    import nibabel as nib
    import numpy as np
    from scipy import ndimage
    
    print("\n" + "=" * 70)
    print("CREATING SHRUNK GRAFT SEGMENTS")
    print("=" * 70)
    
    groups = get_subject_groups()
    
    # 3D ball structuring element (radius 1 = 6-connectivity erosion)
    struct = ndimage.generate_binary_structure(3, 1)
    
    # All 9 segments
    all_segments = [
        "anterior_dorsal", "anterior_medial", "anterior_ventral",
        "central_dorsal", "central_medial", "central_ventral",
        "posterior_dorsal", "posterior_medial", "posterior_ventral"
    ]
    
    for group in ['AP', 'APTR']:
        print(f"\n--- Creating shrunk segments for {group} ---")
        
        for subj_info in groups[group]:
            subj = subj_info['id']
            ses = subj_info['session']
            
            seg_dir = deriv_dir / subj / ses / "dwi" / "graft_segments_fod"
            shrink_dir = seg_dir / "shrink"
            shrink_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"\n  {subj}:")
            
            for hemi in ['left', 'right']:
                for segment in all_segments:
                    seg_name = f"seg9_graft_{hemi}_{segment}_fod.mif"
                    seg_path = seg_dir / seg_name
                    
                    if not seg_path.exists():
                        continue
                    
                    # Load segment via mrconvert
                    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
                        tmp_nii = tmp.name
                    
                    subprocess.run(['mrconvert', str(seg_path), tmp_nii, '-force', '-quiet'], 
                                   check=True, capture_output=True)
                    
                    img = nib.load(tmp_nii)
                    data = img.get_fdata() > 0.5
                    original_voxels = data.sum()
                    
                    if original_voxels == 0:
                        os.remove(tmp_nii)
                        continue
                    
                    # Erode by 1 voxel
                    eroded = ndimage.binary_erosion(data, structure=struct)
                    eroded_voxels = eroded.sum()
                    
                    # Save shrunk segment
                    out_name = f"seg9_graft_{hemi}_{segment}_fod_shrink.mif"
                    out_path = shrink_dir / out_name
                    
                    # Save as nifti then convert to mif
                    eroded_img = nib.Nifti1Image(eroded.astype(np.float32), img.affine, img.header)
                    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp2:
                        tmp_out = tmp2.name
                    nib.save(eroded_img, tmp_out)
                    
                    subprocess.run(['mrconvert', tmp_out, str(out_path), '-force', '-quiet'],
                                   check=True, capture_output=True)
                    
                    # Cleanup
                    os.remove(tmp_nii)
                    os.remove(tmp_out)
                    
                    pct = 100 * eroded_voxels / original_voxels if original_voxels > 0 else 0
                    print(f"    {hemi} {segment}: {int(original_voxels)} -> {int(eroded_voxels)} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("SHRUNK SEGMENTS CREATED SUCCESSFULLY")
    print("=" * 70)


def run_tckedit(wholebrain_tck, output_tck, seed_roi, target_roi, dry_run=False, extra_includes=None):
    """
    Run tckedit to extract passthrough streamlines.
    
    Args:
        wholebrain_tck: Path to whole-brain tractogram
        output_tck: Path for output passthrough tract
        seed_roi: Path to seed ROI (segment)
        target_roi: Path to target ROI
        dry_run: If True, print command without executing
        extra_includes: Optional list of additional ROI paths to add as -include
                        (e.g., thalamus as waypoint for cerebellum tracts)
        
    Returns:
        int: Number of streamlines in output, or -1 on error
    """
    cmd = [
        "tckedit",
        str(wholebrain_tck),
        str(output_tck),
        "-include", str(seed_roi),
    ]
    # Add extra includes (e.g., thalamus waypoint for cerebellum tracts)
    if extra_includes:
        for inc in extra_includes:
            cmd.extend(["-include", str(inc)])
    cmd.extend(["-include", str(target_roi), "-force"])
    
    if dry_run:
        print(f"        [DRY RUN] {' '.join(cmd)}")
        return 0
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode != 0:
            print(f"        ERROR: tckedit failed: {result.stderr[:200]}")
            return -1
        
        # Get streamline count
        info_result = subprocess.run(
            ["tckinfo", str(output_tck), "-count"],
            capture_output=True, text=True
        )
        for line in info_result.stdout.split('\n'):
            if 'count:' in line and 'total_count' not in line:
                count = int(line.split()[-1])
                return count
        return 0
        
    except subprocess.TimeoutExpired:
        print(f"        ERROR: tckedit timed out (>600s)")
        return -1
    except Exception as e:
        print(f"        ERROR: {e}")
        return -1


def process_ctrl_subject(subject_id, session, deriv_dir, hemi="left", dry_run=False):
    """
    Process CTRL subject: extract isocortex -> target passthrough tracts.
    
    CTRL uses 5 A-P isocortex segments and Allen atlas ROIs for targets.
    """
    base_dir = deriv_dir / subject_id / session
    dwi_dir = base_dir / "dwi"
    wholebrain_tck = dwi_dir / "wholebrain" / "wholebrain_1000k.tck"
    out_dir = dwi_dir / "wholebrain" / "passthrough_tracts"
    seg_dir = dwi_dir / "isocortex_segments_v2"
    roi_dir = base_dir / "atlas" / "rois" / "fod_space"
    
    if not wholebrain_tck.exists():
        print(f"    ERROR: Wholebrain tractogram not found: {wholebrain_tck}")
        return {}
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 5 isocortex segments (A-P only, no D-V split)
    segments = ["anterior", "anterior_central", "central", "posterior_central", "posterior"]
    
    # 4 targets (using Allen atlas naming)
    targets = ["thalamus", "striatum", "cerebellum", "olfactory_cortex"]
    
    results = {}
    total = len(segments) * len(targets)
    current = 0
    
    for segment in segments:
        seed_roi = seg_dir / f"isocortex_{hemi}_{segment}.nii.gz"
        
        if not seed_roi.exists():
            print(f"    WARNING: Seed ROI not found: {seed_roi.name}")
            continue
        
        for target in targets:
            current += 1
            target_roi = roi_dir / f"roi_{target}_{hemi}.nii.gz"
            output_tck = out_dir / f"isocortex_{hemi}_{segment}_to_{target}_{hemi}.tck"
            
            key = f"isocortex_{hemi}_{segment}_to_{target}_{hemi}"
            
            if not target_roi.exists():
                print(f"    [{current}/{total}] {key}: TARGET NOT FOUND")
                results[key] = -1
                continue
            
            # For cerebellum tracts, add thalamus as waypoint include
            extra_includes = None
            if target == "cerebellum":
                thalamus_roi = roi_dir / f"roi_thalamus_{hemi}.nii.gz"
                if thalamus_roi.exists():
                    extra_includes = [thalamus_roi]
                    print(f"    [{current}/{total}] {key} (+thalamus include)...", end=" ", flush=True)
                else:
                    print(f"    [{current}/{total}] {key}: WARNING thalamus ROI not found, extracting without")
            
            if output_tck.exists() and not dry_run:
                # Get existing count
                info_result = subprocess.run(
                    ["tckinfo", str(output_tck), "-count"],
                    capture_output=True, text=True
                )
                for line in info_result.stdout.split('\n'):
                    if 'count:' in line and 'total_count' not in line:
                        count = int(line.split()[-1])
                        print(f"    [{current}/{total}] {key}: SKIP (exists, {count} streamlines)")
                        results[key] = count
                        break
                continue
            
            if extra_includes is None:
                print(f"    [{current}/{total}] {key}...", end=" ", flush=True)
            count = run_tckedit(wholebrain_tck, output_tck, seed_roi, target_roi, dry_run, extra_includes=extra_includes)
            
            if count >= 0:
                print(f"{count} streamlines")
                results[key] = count
            else:
                results[key] = -1
    
    return results


def process_apaptr_subject(subject_id, session, group, deriv_dir, hemi="left", dry_run=False, use_shrink=False):
    """
    Process AP or APTR subject: extract graft dorsal -> target passthrough tracts.
    
    AP/APTR uses 3 dorsal graft segments (from seg9) and manual ROIs.
    
    Args:
        use_shrink: If True, use shrunk segment ROIs from graft_segments_fod/shrink/
                    and save outputs to passthrough_tracts_shrink/ folder
    """
    base_dir = deriv_dir / subject_id / session
    dwi_dir = base_dir / "dwi"
    wholebrain_tck = dwi_dir / "wholebrain" / "wholebrain_1000k.tck"
    
    # Use different output folder for shrink mode
    out_folder = "passthrough_tracts_shrink" if use_shrink else "passthrough_tracts"
    out_dir = dwi_dir / "wholebrain" / out_folder
    
    # Use shrunk segments from shrink subfolder if use_shrink is True
    if use_shrink:
        seg_dir = dwi_dir / "graft_segments_fod" / "shrink"
        if not seg_dir.exists():
            print(f"    ERROR: Shrunk segments folder not found: {seg_dir}")
            print(f"    Run with --create-shrink-segments first!")
            return {}
    else:
        seg_dir = dwi_dir / "graft_segments_fod"
    
    roi_dir = base_dir / "manual_roi" / "fod_space"
    
    if not wholebrain_tck.exists():
        print(f"    ERROR: Wholebrain tractogram not found: {wholebrain_tck}")
        return {}
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 3 dorsal segments only (A-P positions × dorsal)
    dorsal_segments = ["anterior_dorsal", "central_dorsal", "posterior_dorsal"]
    
    # 4 targets (using manual ROI naming - note capitalization)
    targets = ["Thalamus", "Striatum", "Palecortex", "cerebellum"]
    
    # Use seg9 only (all subjects have seg9, not all have seg6)
    seg_modes = [9]
    
    results = {}
    total = len(seg_modes) * len(dorsal_segments) * len(targets)
    current = 0
    
    # Segment file suffix differs for shrink mode
    seg_suffix = "_shrink.mif" if use_shrink else ".mif"
    
    for seg_mode in seg_modes:
        for segment in dorsal_segments:
            seg_name = f"seg{seg_mode}_graft_{hemi}_{segment}_fod{seg_suffix}"
            seed_roi = seg_dir / seg_name
            
            if not seed_roi.exists():
                print(f"    WARNING: Seed ROI not found: {seed_roi.name}")
                continue
            
            for target in targets:
                current += 1
                target_roi = roi_dir / f"roi_{target}_{hemi}_flip_fod.mif"
                output_tck = out_dir / f"seg{seg_mode}_graft_{hemi}_{segment}_to_{target}_{hemi}.tck"
                
                key = f"seg{seg_mode}_graft_{hemi}_{segment}_to_{target}_{hemi}"
                
                if not target_roi.exists():
                    print(f"    [{current}/{total}] {key}: TARGET NOT FOUND")
                    results[key] = -1
                    continue
                
                # For cerebellum tracts, add thalamus as waypoint include
                extra_includes = None
                if target == "cerebellum":
                    thalamus_roi = roi_dir / f"roi_Thalamus_{hemi}_flip_fod.mif"
                    if thalamus_roi.exists():
                        extra_includes = [thalamus_roi]
                    else:
                        print(f"    [{current}/{total}] {key}: WARNING thalamus ROI not found for cerebellum waypoint")
                
                if output_tck.exists() and not dry_run:
                    # Get existing count
                    info_result = subprocess.run(
                        ["tckinfo", str(output_tck), "-count"],
                        capture_output=True, text=True
                    )
                    for line in info_result.stdout.split('\n'):
                        if 'count:' in line and 'total_count' not in line:
                            count = int(line.split()[-1])
                            print(f"    [{current}/{total}] {key}: SKIP (exists, {count} streamlines)")
                            results[key] = count
                            break
                    continue
                
                inc_label = " (+thalamus)" if extra_includes else ""
                print(f"    [{current}/{total}] {key}{inc_label}...", end=" ", flush=True)
                count = run_tckedit(wholebrain_tck, output_tck, seed_roi, target_roi, dry_run, extra_includes=extra_includes)
                
                if count >= 0:
                    print(f"{count} streamlines")
                    results[key] = count
                else:
                    results[key] = -1
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract passthrough streamlines from whole-brain tractograms',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--group', type=str, choices=['CTRL', 'AP', 'APTR'],
                        help='Process specific group only')
    parser.add_argument('--subject', type=str, help='Process specific subject only')
    parser.add_argument('--hemi', type=str, choices=['left', 'right'], default='left',
                        help='Hemisphere to process (default: left)')
    parser.add_argument('--dry-run', action='store_true', help='Preview commands without executing')
    parser.add_argument('--shrink', action='store_true',
                        help='Use shrunk graft segments (eroded by 1 voxel) - AP/APTR only')
    parser.add_argument('--create-shrink-segments', action='store_true',
                        help='Create shrunk versions of all graft segments for AP/APTR')
    
    args = parser.parse_args()
    
    # Load configuration
    cfg = load_config()
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_flip1'])
    
    # Create shrunk segments mode
    if args.create_shrink_segments:
        create_shrunk_segments(deriv_dir)
        return
    
    # Shrink mode only applies to AP/APTR
    if args.shrink and args.group == 'CTRL':
        print("ERROR: --shrink option only applies to AP/APTR groups (CTRL has no graft)")
        sys.exit(1)
    
    shrink_label = " (SHRUNK SEGMENTS)" if args.shrink else ""
    print("=" * 70)
    print(f"PASSTHROUGH TRACT EXTRACTION{shrink_label}")
    print(f"Hemisphere: {args.hemi.upper()}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    # Get subjects to process
    subjects = get_subjects(group=args.group, subject=args.subject)
    
    # If shrink mode with no group specified, filter to AP/APTR only
    if args.shrink and not args.group:
        subjects = [(s, ses, g) for (s, ses, g) in subjects if g in ['AP', 'APTR']]
    
    print(f"\nProcessing {len(subjects)} subjects ({args.hemi} hemisphere)...")
    
    all_results = {}
    
    for subject_id, session, group in subjects:
        print(f"\n{'='*70}")
        print(f"Subject: {subject_id} ({group})")
        print(f"Session: {session}")
        print("=" * 70)
        
        if group == "CTRL":
            results = process_ctrl_subject(subject_id, session, deriv_dir, args.hemi, args.dry_run)
        else:  # AP or APTR
            results = process_apaptr_subject(subject_id, session, group, deriv_dir, args.hemi, args.dry_run, args.shrink)
        
        all_results[subject_id] = results
        
        # Summary for this subject
        if results:
            success = sum(1 for v in results.values() if v >= 0)
            total_streamlines = sum(v for v in results.values() if v > 0)
            print(f"\n    Summary: {success}/{len(results)} tracts extracted, {total_streamlines} total streamlines")
    
    # Final summary
    print("\n" + "=" * 70)
    print("EXTRACTION COMPLETE")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    
    for subj, results in all_results.items():
        if results:
            success = sum(1 for v in results.values() if v >= 0)
            total = sum(v for v in results.values() if v > 0)
            print(f"  {subj}: {success} tracts, {total} streamlines")


if __name__ == "__main__":
    main()
