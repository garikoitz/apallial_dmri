#!/usr/bin/env python3
"""
Generic tract visualization script.

Opens mrview with tracts and overlays for any subject/group.

Usage:
    python view_tracts.py --subject sub-129-921-R1F             # View specific subject
    python view_tracts.py --subject sub-129-921-R1F --tract Isocortex_Striatum
    python view_tracts.py --group CTRL --template               # View on group template
    python view_tracts.py --list-tracts                         # List available tracts
"""

import argparse
import subprocess
from pathlib import Path

from config_loader import (
    load_config, get_paths, get_subject_groups, get_tract_definitions
)


def get_subject_info(subject_id):
    """Get session and group for a subject."""
    groups = get_subject_groups()
    for group_name, subjects in groups.items():
        for s in subjects:
            if s['id'] == subject_id:
                return s['session'], group_name
    return None, None


def list_available_tracts(subject_id=None):
    """List available tract definitions."""
    tract_defs = get_tract_definitions()
    cfg = load_config()
    top_n = cfg.get('top_tracts', {}).get('count', 200)
    
    print("Available tract definitions:")
    for t in tract_defs:
        print(f"  - {t['name']}: {t['seed']} -> {t['target']}")
    
    print(f"\nFile naming convention:")
    print(f"  Raw:     tracks_<name>_<hemisphere>.tck")
    print(f"  Cleaned: tracks_<name>_<hemisphere>_cleaned.tck")
    print(f"  Top-N:   tracks_<name>_<hemisphere>_cleaned_top{top_n}.tck")
    
    if subject_id:
        paths = get_paths()
        session, _ = get_subject_info(subject_id)
        if session:
            dwi_dir = Path(paths['derivatives_dir']) / subject_id / session / "dwi"
            print(f"\nTract files for {subject_id}:")
            for f in sorted(dwi_dir.glob("tracks_*.tck")):
                print(f"  - {f.name}")


def view_subject_tracts(subject_id, tract_name=None, hemisphere='both', cleaned=True, top_n=None):
    """Open mrview for a specific subject."""
    paths = get_paths()
    cfg = load_config()
    
    session, group_name = get_subject_info(subject_id)
    if not session:
        print(f"ERROR: Subject {subject_id} not found")
        return
    
    deriv_dir = Path(paths['derivatives_dir'])
    dwi_dir = deriv_dir / subject_id / session / "dwi"
    
    # Background image
    bg_image = dwi_dir / "t2w_highres_flip.mif"
    if not bg_image.exists():
        bg_image = dwi_dir / "b0_final.mif"
    
    if not bg_image.exists():
        print(f"ERROR: No background image found in {dwi_dir}")
        return
    
    # Build mrview command
    cmd = ["mrview", str(bg_image)]
    
    # Get tract files to load
    tract_defs = get_tract_definitions()
    if tract_name:
        tract_defs = [t for t in tract_defs if t['name'] == tract_name]
    
    if top_n is None:
        top_n = cfg.get('top_tracts', {}).get('count', 200)
    
    hemispheres = ['left', 'right'] if hemisphere == 'both' else [hemisphere]
    
    for t in tract_defs:
        for hemi in hemispheres:
            # Try different file variants
            variants = []
            if cleaned:
                variants.append(f"tracks_{t['name']}_{hemi}_cleaned_top{top_n}.tck")
                variants.append(f"tracks_{t['name']}_{hemi}_cleaned.tck")
            variants.append(f"tracks_{t['name']}_{hemi}_top{top_n}.tck")
            variants.append(f"tracks_{t['name']}_{hemi}.tck")
            
            for v in variants:
                tract_file = dwi_dir / v
                if tract_file.exists():
                    cmd.extend(["-tractography.load", str(tract_file)])
                    print(f"  Loading: {v}")
                    break
    
    # Add atlas overlay if available
    atlas_file = dwi_dir / "atlas" / "atlas_template_warped.nii.gz"
    if atlas_file.exists():
        cmd.extend(["-overlay.load", str(atlas_file), "-overlay.opacity", "0.3"])
    
    print(f"\nOpening mrview for {subject_id}...")
    subprocess.Popen(cmd)


def view_template_tracts(group_name, tract_name=None):
    """Open mrview for group template with transformed tracts."""
    paths = get_paths()
    cfg = load_config()
    
    deriv_dir = Path(paths['derivatives_dir'])
    template_dir = deriv_dir / "group_templates" / group_name
    tracts_dir = template_dir / "tracts"
    
    template_file = template_dir / f"{group_name}_template.nii.gz"
    if not template_file.exists():
        print(f"ERROR: Template not found: {template_file}")
        return
    
    cmd = ["mrview", str(template_file)]
    
    # Find all tract files in template space
    if tracts_dir.exists():
        for tract_file in sorted(tracts_dir.glob("*.tck")):
            if tract_name and tract_name not in tract_file.name:
                continue
            cmd.extend(["-tractography.load", str(tract_file)])
            print(f"  Loading: {tract_file.name}")
    
    print(f"\nOpening mrview for {group_name} template...")
    subprocess.Popen(cmd)


def main():
    parser = argparse.ArgumentParser(
        description='Visualize tracts in mrview',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--subject', type=str, help='Subject to view')
    parser.add_argument('--group', type=str, help='Group template to view')
    parser.add_argument('--tract', type=str, help='Specific tract type')
    parser.add_argument('--hemisphere', choices=['left', 'right', 'both'],
                        default='both', help='Hemisphere')
    parser.add_argument('--raw', action='store_true', help='Use raw (uncleaned) tracts')
    parser.add_argument('--top-n', type=int, help='Top N tracks to view')
    parser.add_argument('--template', action='store_true', help='View on group template')
    parser.add_argument('--list-tracts', action='store_true', help='List available tracts')
    
    args = parser.parse_args()
    
    if args.list_tracts:
        list_available_tracts(args.subject)
        return
    
    if args.template and args.group:
        view_template_tracts(args.group, args.tract)
    elif args.subject:
        view_subject_tracts(
            args.subject,
            tract_name=args.tract,
            hemisphere=args.hemisphere,
            cleaned=not args.raw,
            top_n=args.top_n
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
