#!/usr/bin/env python3
"""
Generic tract transformation script.

Transforms tracts to group template space using ANTs warps.
Supports all subjects, groups, and tract types.

Usage:
    python transform_tracts.py                                    # All tracts to templates
    python transform_tracts.py --group CTRL                       # CTRL group only
    python transform_tracts.py --subject sub-129-921-R1F          # Specific subject
    python transform_tracts.py --tract Isocortex_Striatum         # Specific tract
    python transform_tracts.py --dry-run                          # Preview commands
"""

import argparse
import subprocess
import numpy as np
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups, get_tract_definitions
)


def get_group_for_subject(subject_id):
    """Find which group a subject belongs to."""
    groups = get_subject_groups()
    for group_name, subjects in groups.items():
        for s in subjects:
            if s['id'] == subject_id:
                return group_name
    return None


def get_subjects(group=None, subject=None):
    """Get list of (subject_id, session, group) tuples."""
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
            raise ValueError(f"Group {group} not found")
        for s in groups[group]:
            subjects.append((s['id'], s['session'], group))
    else:
        for g, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session'], g))
    
    return subjects


def get_tracts(tract_name=None):
    """Get tract definitions, optionally filtered."""
    cfg = load_config()
    tract_defs = get_tract_definitions()
    top_n = cfg.get('top_tracts', {}).get('count', 200)
    
    tracts = []
    for t in tract_defs:
        if tract_name and t['name'] != tract_name:
            continue
        for hemi in ['left', 'right']:
            # Include both cleaned and raw variants
            tracts.append(f"tracks_{t['name']}_{hemi}_cleaned_top{top_n}")
            tracts.append(f"tracks_{t['name']}_{hemi}_top{top_n}")
    
    return tracts


def convert_warp_to_mrtrix(warp_file, affine_file, reference, output_deform):
    """Convert ANTs warp to MRtrix deformation field."""
    # Extract first 3 volumes (displacement field)
    disp_4d = output_deform.replace('.mif', '_disp4d.mif')
    cmd = ["mrconvert", str(warp_file), "-axes", "0,1,2,4", str(disp_4d), "-force"]
    subprocess.run(cmd, check=True)
    
    # Convert displacement to deformation
    cmd = ["warpconvert", str(disp_4d), "displacement2deformation", str(output_deform), "-force"]
    subprocess.run(cmd, check=True)
    
    return True


def transform_tract(subject_id, session, group_name, tract_file, dry_run):
    """Transform a single tract to group template space."""
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    template_dir = deriv_dir / "group_templates" / group_name
    transforms_dir = template_dir / "transforms" / subject_id
    output_dir = template_dir / "tracts"
    
    input_tract = dwi_dir / f"{tract_file}.tck"
    output_tract = output_dir / f"{subject_id}_{tract_file}.tck"
    
    if not input_tract.exists():
        return None  # Skip silently if tract doesn't exist
    
    print(f"    {tract_file}")
    
    if dry_run:
        print(f"      [DRY RUN] Would transform to {output_tract.name}")
        return True
    
    # Check for warp files
    warp_file = transforms_dir / "subject_to_template_Warp.nii.gz"
    affine_file = transforms_dir / "subject_to_template_0GenericAffine.mat"
    template_file = template_dir / f"{group_name}_template.nii.gz"
    
    if not warp_file.exists():
        print(f"      SKIP: Warp file not found")
        return False
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert ANTs warp to MRtrix format
    deform_file = transforms_dir / "deformation_field.mif"
    if not deform_file.exists():
        convert_warp_to_mrtrix(warp_file, affine_file, template_file, str(deform_file))
    
    # Transform tract
    cmd = [
        "tcktransform",
        str(input_tract),
        str(deform_file),
        str(output_tract),
        "-force"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"      ERROR: {result.stderr[:200]}")
        return False
    
    print(f"      ✓ -> {output_tract.name}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Transform tracts to group template space',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--group', type=str, help='Process specific group')
    parser.add_argument('--subject', type=str, help='Process specific subject')
    parser.add_argument('--tract', type=str, help='Process specific tract type')
    parser.add_argument('--dry-run', action='store_true', help='Preview without executing')
    
    args = parser.parse_args()
    
    subjects = get_subjects(group=args.group, subject=args.subject)
    tract_files = get_tracts(tract_name=args.tract)
    
    print(f"Transforming tracts to group templates")
    print(f"Subjects: {len(subjects)}, Tract variants: {len(tract_files)}")
    
    results = []
    for subject_id, session, group_name in subjects:
        print(f"\n{'='*60}")
        print(f"Subject: {subject_id} -> {group_name} template")
        print(f"{'='*60}")
        
        for tract_file in tract_files:
            result = transform_tract(
                subject_id, session, group_name, tract_file, args.dry_run
            )
            if result is not None:
                results.append((subject_id, tract_file, result))
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success_count = sum(1 for r in results if r[2])
    print(f"  Transformed: {success_count}/{len(results)}")


if __name__ == "__main__":
    main()
