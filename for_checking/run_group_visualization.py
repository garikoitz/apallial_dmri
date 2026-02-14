#!/usr/bin/env python3
"""
Complete group visualization pipeline.

Generates, processes, transforms, and visualizes all tracts for groups:
- CTRL: Isocortex->Striatum, Isocortex->Thalamus (from atlas ROIs)
- AP/APTR: Graft->Striatum, Graft->Thalamus (from manual ROIs)

Creates 3 mrview windows (one per group) with:
- 16 tracts per window (4 subjects × 2 tract types × 2 hemispheres)
- T2w_highres as background
- ROIs overlaid on template

Usage:
    python run_group_visualization.py                    # Full pipeline
    python run_group_visualization.py --skip-generate    # Skip tract generation
    python run_group_visualization.py --skip-sift        # Skip SIFT2
    python run_group_visualization.py --skip-transform   # Skip transform to template
    python run_group_visualization.py --visualize-only   # Only open mrview
"""

import argparse
import subprocess
import os
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups, get_tract_definitions
)


def get_group_tracts(group_name):
    """Get tract definitions applicable to a group."""
    tract_defs = get_tract_definitions()
    applicable = []
    for t in tract_defs:
        groups = t.get('groups', [])
        if not groups or group_name in groups:
            applicable.append(t)
    return applicable


def run_command(cmd, description, dry_run=False):
    """Run a command with logging."""
    print(f"  {description}...")
    if dry_run:
        print(f"    [DRY RUN] {' '.join(cmd[:6])}...")
        return True
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    ERROR: {result.stderr[:300]}")
        return False
    return True


def generate_tracts_for_group(group_name, dry_run=False):
    """Generate all tracts for a group."""
    print(f"\n{'='*60}")
    print(f"Generating tracts for group: {group_name}")
    print(f"{'='*60}")
    
    groups = get_subject_groups()
    tracts = get_group_tracts(group_name)
    paths = get_paths()
    cfg = load_config()
    params = cfg.get('tractography', {})
    
    for subject in groups[group_name]:
        subject_id = subject['id']
        session = subject['session']
        
        deriv_dir = Path(paths['derivatives_dir'])
        dwi_dir = deriv_dir / subject_id / session / "dwi"
        manual_roi_dir = deriv_dir / subject_id / session / "manual_roi"
        
        fod_file = dwi_dir / "fod_wm.mif"
        mask_file = dwi_dir / "mask.mif"
        
        print(f"\n  Subject: {subject_id}")
        
        for tract in tracts:
            roi_source = tract.get('roi_source', 'atlas')
            
            for hemi in ['left', 'right']:
                output_file = dwi_dir / f"tracks_{tract['name']}_{hemi}.tck"
                
                # Get ROI paths
                if roi_source == 'manual':
                    seed_roi = manual_roi_dir / f"roi_{tract['seed']}_{hemi}.nii.gz"
                    include_roi = manual_roi_dir / f"roi_{tract['target']}_{hemi}.nii.gz"
                else:
                    roi_dir = dwi_dir / "atlas" / "rois" / "fod_space"
                    if not roi_dir.exists():
                        roi_dir = dwi_dir / "atlas" / "rois"
                    seed_roi = roi_dir / f"roi_{tract['seed'].lower()}_{hemi}.nii.gz"
                    if not seed_roi.exists():
                        seed_roi = roi_dir / f"roi_{tract['seed'].lower()}_{hemi}_masked.nii.gz"
                    include_roi = roi_dir / f"roi_{tract['target'].lower()}_{hemi}.nii.gz"
                    if not include_roi.exists():
                        include_roi = roi_dir / f"roi_{tract['target'].lower()}_{hemi}_masked.nii.gz"
                
                # Skip if already exists
                if output_file.exists() and not dry_run:
                    info_result = subprocess.run(
                        ["tckinfo", "-count", str(output_file)],
                        capture_output=True, text=True
                    )
                    count = int(info_result.stdout.strip()) if info_result.returncode == 0 else 0
                    if count >= 1000:
                        print(f"    ✓ {tract['name']}_{hemi} exists ({count} tracks)")
                        continue
                
                # Check inputs exist
                if not seed_roi.exists():
                    print(f"    ✗ SKIP {tract['name']}_{hemi}: seed ROI not found: {seed_roi.name}")
                    continue
                if not include_roi.exists():
                    print(f"    ✗ SKIP {tract['name']}_{hemi}: include ROI not found: {include_roi.name}")
                    continue
                
                cmd = [
                    "tckgen",
                    str(fod_file),
                    str(output_file),
                    "-algorithm", params.get('algorithm', 'iFOD2'),
                    "-select", str(params.get('select', 5000)),
                    "-step", str(params.get('step_size', 0.02)),
                    "-angle", str(params.get('angle', 60)),
                    "-minlength", str(params.get('min_length', 1.5)),
                    "-maxlength", str(params.get('max_length', 15)),
                    "-cutoff", str(params.get('cutoff', 0.05)),
                    "-seed_image", str(seed_roi),
                    "-include", str(include_roi),
                ]
                
                if mask_file.exists():
                    cmd.extend(["-mask", str(mask_file)])
                cmd.append("-force")
                
                run_command(cmd, f"tckgen {tract['name']}_{hemi}", dry_run)


def run_sift2_for_group(group_name, dry_run=False):
    """Run SIFT2 on all tracts for a group."""
    print(f"\n{'='*60}")
    print(f"Running SIFT2 for group: {group_name}")
    print(f"{'='*60}")
    
    groups = get_subject_groups()
    tracts = get_group_tracts(group_name)
    paths = get_paths()
    
    for subject in groups[group_name]:
        subject_id = subject['id']
        session = subject['session']
        
        deriv_dir = Path(paths['derivatives_dir'])
        dwi_dir = deriv_dir / subject_id / session / "dwi"
        
        fod_file = dwi_dir / "fod_wm.mif"
        
        print(f"\n  Subject: {subject_id}")
        
        for tract in tracts:
            for hemi in ['left', 'right']:
                input_file = dwi_dir / f"tracks_{tract['name']}_{hemi}.tck"
                weights_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_sift2_weights.csv"
                
                if not input_file.exists():
                    continue
                
                if weights_file.exists():
                    print(f"    ✓ {tract['name']}_{hemi} SIFT2 weights exist")
                    continue
                
                cmd = [
                    "tcksift2",
                    str(input_file),
                    str(fod_file),
                    str(weights_file),
                    "-force"
                ]
                
                run_command(cmd, f"tcksift2 {tract['name']}_{hemi}", dry_run)


def clean_tracts_with_matlab(group_name, dry_run=False):
    """Clean tracts using MATLAB AFQ_removeFiberOutliers."""
    print(f"\n{'='*60}")
    print(f"MATLAB Cleaning tracts for group: {group_name}")
    print(f"{'='*60}")
    
    cfg = load_config()
    matlab_cfg = cfg.get('matlab_cleaning', {})
    paths = get_paths()
    groups = get_subject_groups()
    tracts = get_group_tracts(group_name)
    
    matlab_bin = matlab_cfg.get('matlab_bin', '/Applications/MATLAB_R2025b.app/bin/matlab')
    afq_path = matlab_cfg.get('afq_path', '/Users/glerma/soft/gitlab/rtp2-pipeline')
    max_dist = matlab_cfg.get('max_dist', 3)
    max_len = matlab_cfg.get('max_len', 3)
    num_nodes = matlab_cfg.get('num_nodes', 100)
    core_method = matlab_cfg.get('core_method', 'mean')
    verbose = matlab_cfg.get('verbose', 0)
    max_iter = matlab_cfg.get('max_iter', 10)
    show_viz = matlab_cfg.get('show_viz', 0)
    
    for subject in groups[group_name]:
        subject_id = subject['id']
        session = subject['session']
        
        deriv_dir = Path(paths['derivatives_dir'])
        dwi_dir = deriv_dir / subject_id / session / "dwi"
        
        print(f"\n  Subject: {subject_id}")
        
        for tract in tracts:
            for hemi in ['left', 'right']:
                input_file = dwi_dir / f"tracks_{tract['name']}_{hemi}.tck"
                cleaned_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned.tck"
                sf_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_SF.tck"
                
                if not input_file.exists():
                    print(f"    ✗ SKIP {tract['name']}_{hemi}: input not found")
                    continue
                
                if cleaned_file.exists():
                    print(f"    ✓ {tract['name']}_{hemi}_cleaned exists")
                    continue
                
                if dry_run:
                    print(f"    [DRY RUN] MATLAB clean {tract['name']}_{hemi}")
                    continue
                
                # Build MATLAB command
                matlab_cmd = f"""
addpath(genpath('{afq_path}'));
dwidir = '{dwi_dir}';
tckfilename = '{input_file.name}';
fg = fgRead(fullfile(dwidir, tckfilename));
[fg_clean, ~, SF] = AFQ_removeFiberOutliers(fg, {max_dist}, {max_len}, {num_nodes}, '{core_method}', {verbose}, {max_iter}, {show_viz});
fg_SF = dtiNewFiberGroup([fg.name '_SF'], [0 0 155], [], [], SF.fibers);
AFQ_fgWrite(fg_clean, fullfile(dwidir, '{cleaned_file.name}'), 'tck');
AFQ_fgWrite(fg_SF, fullfile(dwidir, '{sf_file.name}'), 'tck');
fprintf('Cleaned: %s -> %d fibers (removed %d)\\n', tckfilename, length(fg_clean.fibers), length(fg.fibers) - length(fg_clean.fibers));
exit;
"""
                
                cmd = [
                    matlab_bin,
                    "-nodisplay", "-nosplash", "-r",
                    matlab_cmd.strip().replace('\n', ' ')
                ]
                
                print(f"    Cleaning {tract['name']}_{hemi}...")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                
                if result.returncode != 0 or not cleaned_file.exists():
                    print(f"      ERROR: MATLAB cleaning failed")
                    print(f"      {result.stderr[:200] if result.stderr else 'No error output'}")
                else:
                    # Get track count
                    info_result = subprocess.run(
                        ["tckinfo", "-count", str(cleaned_file)],
                        capture_output=True, text=True
                    )
                    count = info_result.stdout.strip() if info_result.returncode == 0 else "?"
                    print(f"      ✓ Cleaned -> {count} fibers")


def extract_top_tracts_for_group(group_name, top_n=200, dry_run=False):
    """Extract top N tracts by SIFT2 weight for a group (from cleaned tracts)."""
    print(f"\n{'='*60}")
    print(f"Extracting top {top_n} tracts for group: {group_name}")
    print(f"{'='*60}")
    
    import numpy as np
    
    groups = get_subject_groups()
    tracts = get_group_tracts(group_name)
    paths = get_paths()
    cfg = load_config()
    fod_file_name = "fod_wm.mif"
    
    for subject in groups[group_name]:
        subject_id = subject['id']
        session = subject['session']
        
        deriv_dir = Path(paths['derivatives_dir'])
        dwi_dir = deriv_dir / subject_id / session / "dwi"
        fod_file = dwi_dir / fod_file_name
        
        print(f"\n  Subject: {subject_id}")
        
        for tract in tracts:
            for hemi in ['left', 'right']:
                # Prefer cleaned tracts
                cleaned_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned.tck"
                raw_file = dwi_dir / f"tracks_{tract['name']}_{hemi}.tck"
                
                # Use cleaned if available, otherwise fall back to raw
                if cleaned_file.exists():
                    input_file = cleaned_file
                    output_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned_top{top_n}.tck"
                    weights_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned_sift2_weights.csv"
                    base_name = f"{tract['name']}_{hemi}_cleaned"
                else:
                    input_file = raw_file
                    output_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_top{top_n}.tck"
                    weights_file = dwi_dir / f"tracks_{tract['name']}_{hemi}_sift2_weights.csv"
                    base_name = f"{tract['name']}_{hemi}"
                    print(f"    ⚠️  Using uncleaned tracts for {base_name}")
                
                if not input_file.exists():
                    continue
                
                # Run SIFT2 on cleaned tract if weights don't exist
                if not weights_file.exists():
                    if not fod_file.exists():
                        print(f"    ✗ SKIP {base_name}: FOD not found")
                        continue
                    print(f"    Running SIFT2 on {base_name}...")
                    cmd = ["tcksift2", str(input_file), str(fod_file), str(weights_file), "-force"]
                    if not dry_run:
                        subprocess.run(cmd, capture_output=True, check=True)
                
                if output_file.exists():
                    print(f"    ✓ {base_name}_top{top_n} exists")
                    continue
                
                if dry_run:
                    print(f"    [DRY RUN] Extract top {top_n} for {base_name}")
                    continue
                
                # Load weights and get top N indices
                try:
                    weights = np.loadtxt(weights_file)
                    if len(weights) < top_n:
                        top_n_actual = len(weights)
                    else:
                        top_n_actual = top_n
                    top_indices = np.argsort(weights)[-top_n_actual:]
                    
                    # Save indices to temp file
                    indices_file = dwi_dir / f"temp_indices_{tract['name']}_{hemi}.txt"
                    np.savetxt(str(indices_file), top_indices, fmt='%d')
                    
                    # Extract using tckedit
                    cmd = [
                        "tckedit",
                        str(input_file),
                        str(output_file),
                        "-include_indices", str(indices_file),
                        "-force"
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    # Clean up
                    indices_file.unlink()
                    
                    if result.returncode == 0:
                        print(f"    ✓ Extracted {top_n_actual} tracts for {tract['name']}_{hemi}")
                    else:
                        print(f"    ✗ ERROR: {result.stderr[:200]}")
                except Exception as e:
                    print(f"    ✗ ERROR extracting {tract['name']}_{hemi}: {e}")


def transform_tracts_to_template(group_name, top_n=200, dry_run=False):
    """Transform all top-N tracts to group template."""
    print(f"\n{'='*60}")
    print(f"Transforming tracts to {group_name} template")
    print(f"{'='*60}")
    
    groups = get_subject_groups()
    tracts = get_group_tracts(group_name)
    paths = get_paths()
    
    deriv_dir = Path(paths['derivatives_dir'])
    template_dir = deriv_dir / "group_templates" / group_name
    output_dir = template_dir / "tracts_warped"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    template_file = template_dir / f"{group_name}_template.nii.gz"
    
    for subject in groups[group_name]:
        subject_id = subject['id']
        session = subject['session']
        
        dwi_dir = deriv_dir / subject_id / session / "dwi"
        transforms_dir = template_dir / "transforms" / subject_id
        
        # Check for warp files
        warp_file = transforms_dir / "subject_to_template_Warp.nii.gz"
        if not warp_file.exists():
            # Try alternate naming
            warp_file = transforms_dir / "final_1Warp.nii.gz"
        
        if not warp_file.exists():
            print(f"\n  SKIP {subject_id}: No warp file found")
            continue
        
        print(f"\n  Subject: {subject_id}")
        
        # Create deformation field if needed
        deform_file = output_dir / f"{subject_id}_deformation.mif"
        
        if not deform_file.exists() and not dry_run:
            # Convert ANTs warp to MRtrix deformation
            try:
                # ANTs warp is 4D displacement field
                disp_file = str(deform_file).replace('.mif', '_disp.mif')
                
                cmd = ["mrconvert", str(warp_file), "-axes", "0,1,2,4", disp_file, "-force"]
                subprocess.run(cmd, check=True, capture_output=True)
                
                cmd = ["warpconvert", disp_file, "displacement2deformation", str(deform_file), "-force"]
                subprocess.run(cmd, check=True, capture_output=True)
                
                # Clean up intermediate
                Path(disp_file).unlink()
                print(f"    ✓ Created deformation field")
            except Exception as e:
                print(f"    ✗ Error creating deformation field: {e}")
                continue
        
        # Transform each tract
        for tract in tracts:
            for hemi in ['left', 'right']:
                # Priority order: cleaned_top > raw_top > cleaned > raw
                input_file = None
                for candidate in [
                    dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned_top{top_n}.tck",
                    dwi_dir / f"tracks_{tract['name']}_{hemi}_top{top_n}.tck",
                    dwi_dir / f"tracks_{tract['name']}_{hemi}_cleaned.tck",
                    dwi_dir / f"tracks_{tract['name']}_{hemi}.tck",
                ]:
                    if candidate.exists():
                        input_file = candidate
                        break
                
                if input_file is None:
                    continue
                
                output_file = output_dir / f"{subject_id}_tracks_{tract['name']}_{hemi}.tck"
                
                if output_file.exists():
                    print(f"    ✓ {tract['name']}_{hemi} already warped")
                    continue
                
                if dry_run:
                    print(f"    [DRY RUN] Transform {tract['name']}_{hemi}")
                    continue
                
                cmd = [
                    "tcktransform",
                    str(input_file),
                    str(deform_file),
                    str(output_file),
                    "-force"
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode == 0:
                    print(f"    ✓ Transformed {tract['name']}_{hemi}")
                else:
                    print(f"    ✗ Error: {result.stderr[:150]}")


def create_template_rois(group_name, dry_run=False):
    """Create average/combined ROIs on the group template."""
    print(f"\n{'='*60}")
    print(f"Creating template ROIs for {group_name}")
    print(f"{'='*60}")
    
    # For now, just print what would be done
    # This would warp individual ROIs to template and average them
    print(f"  [SKIP] Template ROI creation not yet implemented")
    print(f"  ROIs can be viewed on individual subject space")


def open_group_visualization(group_name, dry_run=False):
    """Open mrview with all tracts for a group."""
    print(f"\n{'='*60}")
    print(f"Opening visualization for {group_name}")
    print(f"{'='*60}")
    
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    template_dir = deriv_dir / "group_templates" / group_name
    tracts_dir = template_dir / "tracts_warped"
    
    template_file = template_dir / f"{group_name}_template.nii.gz"
    
    if not template_file.exists():
        print(f"  ERROR: Template not found: {template_file}")
        return
    
    # Find all tract files
    tract_files = sorted(tracts_dir.glob("*.tck"))
    
    if not tract_files:
        print(f"  WARNING: No warped tracts found in {tracts_dir}")
        return
    
    print(f"  Found {len(tract_files)} tract files")
    
    # Build mrview command
    cmd = ["mrview", str(template_file)]
    
    # Add tract overlays
    for tck in tract_files:
        cmd.extend(["-tractography.load", str(tck)])
    
    if dry_run:
        print(f"  [DRY RUN] Would open mrview with {len(tract_files)} tracts")
        print(f"  Command: {' '.join(cmd[:8])}...")
    else:
        print(f"  Launching mrview...")
        subprocess.Popen(cmd)


def main():
    parser = argparse.ArgumentParser(
        description='Complete group tract visualization pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group only')
    parser.add_argument('--skip-generate', action='store_true', help='Skip tract generation')
    parser.add_argument('--skip-sift', action='store_true', help='Skip SIFT2')
    parser.add_argument('--skip-clean', action='store_true', help='Skip MATLAB cleaning')
    parser.add_argument('--skip-extract', action='store_true', help='Skip top-N extraction')
    parser.add_argument('--skip-transform', action='store_true', help='Skip transformation')
    parser.add_argument('--visualize-only', action='store_true', help='Only open mrview')
    parser.add_argument('--top-n', type=int, default=200, help='Number of top tracts to extract')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    cfg = load_config()
    matlab_cfg = cfg.get('matlab_cleaning', {})
    cleaning_required = matlab_cfg.get('required', True)
    
    groups_to_process = [args.group] if args.group else ['CTRL', 'AP', 'APTR']
    
    print("="*70)
    print("GROUP TRACT VISUALIZATION PIPELINE")
    print("="*70)
    print(f"Groups: {groups_to_process}")
    print(f"Top-N: {args.top_n}")
    print(f"MATLAB cleaning required: {cleaning_required}")
    
    for group_name in groups_to_process:
        if not args.visualize_only:
            if not args.skip_generate:
                generate_tracts_for_group(group_name, args.dry_run)
            
            if not args.skip_sift:
                run_sift2_for_group(group_name, args.dry_run)
            
            # MATLAB cleaning - required step
            if not args.skip_clean:
                clean_tracts_with_matlab(group_name, args.dry_run)
            elif cleaning_required:
                print(f"\n⚠️  WARNING: MATLAB cleaning is REQUIRED but was skipped!")
                print(f"    Visualization will use uncleaned tracts.")
            
            if not args.skip_extract:
                extract_top_tracts_for_group(group_name, args.top_n, args.dry_run)
            
            if not args.skip_transform:
                transform_tracts_to_template(group_name, args.top_n, args.dry_run)
        
        # Open visualization
        open_group_visualization(group_name, args.dry_run)
    
    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
