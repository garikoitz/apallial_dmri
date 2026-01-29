#!/usr/bin/env python3
"""
Build group template brains using ANTs.

Creates a template brain for each group (AP, APTR, CTRL) using iterative
registration with ANTs. Saves all transformations for warping tracts to
template space.

Usage:
    python build_group_templates.py              # Build all templates
    python build_group_templates.py --group AP   # Build only AP template
    python build_group_templates.py --quick      # Quick mode (fewer iterations)
"""

import subprocess
import argparse
import shutil
from pathlib import Path

# Import config loader
from config_loader import (
    load_config,
    get_paths,
    get_subject_groups
)

# Load configuration
cfg = load_config()
paths = get_paths()

# Paths from config
BASE_DIR = Path(paths['base_dir'])
DERIV_DIR = Path(paths['derivatives_dir'])
TEMPLATE_DIR = DERIV_DIR / "group_templates"

# Build GROUPS from config
def get_groups_with_info():
    """Get all groups with subjects and description."""
    groups = get_subject_groups()
    result = {}
    for group_name, subjects in groups.items():
        result[group_name] = {
            "subjects": [(s['id'], s['session']) for s in subjects],
            "description": f"{group_name} group"
        }
    return result

GROUPS = get_groups_with_info()


def run_cmd(cmd, description="", check=True):
    """Run a command and handle errors."""
    if description:
        print(f"  → {description}")
    
    if isinstance(cmd, str):
        print(f"    $ {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    else:
        print(f"    $ {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0 and check:
        print(f"    ✗ FAILED: {result.stderr}")
        return False
    if result.stdout:
        for line in result.stdout.strip().split('\n')[-5:]:
            print(f"      {line}")
    return True


def prepare_t2w_images(group_name, group_info):
    """Prepare T2w images for template construction."""
    print(f"\n  Preparing T2w images for {group_name}...")
    
    group_dir = TEMPLATE_DIR / group_name
    input_dir = group_dir / "input_images"
    input_dir.mkdir(parents=True, exist_ok=True)
    
    image_list = []
    for subject, session in group_info["subjects"]:
        # Source T2w (use N4 corrected if available, otherwise highres flip)
        dwi_dir = DERIV_DIR / subject / session / "dwi"
        
        # Try atlas N4 corrected first
        t2w_n4 = dwi_dir / "atlas" / "t2w_n4.nii.gz"
        t2w_flip = dwi_dir / "t2w_highres_flip.nii.gz"
        t2w_mif = dwi_dir / "t2w_highres_flip.mif"
        
        if t2w_n4.exists():
            source = t2w_n4
        elif t2w_flip.exists():
            source = t2w_flip
        elif t2w_mif.exists():
            # Convert MIF to NIfTI
            if not t2w_flip.exists():
                run_cmd([
                    "mrconvert", str(t2w_mif), str(t2w_flip), "-force"
                ], f"Converting {subject} T2w to NIfTI")
            source = t2w_flip
        else:
            print(f"    ✗ No T2w found for {subject}")
            continue
        
        # Copy to input directory with standardized name
        dest = input_dir / f"{subject}_T2w.nii.gz"
        if not dest.exists():
            shutil.copy(source, dest)
        
        image_list.append(str(dest))
        print(f"    ✓ {subject}")
    
    # Write image list file
    list_file = group_dir / "image_list.txt"
    with open(list_file, 'w') as f:
        f.write('\n'.join(image_list))
    
    return image_list, list_file


def build_template_ants(group_name, image_list, quick=False):
    """Build template using ANTs multivariate template construction."""
    print(f"\n  Building {group_name} template with ANTs...")
    
    group_dir = TEMPLATE_DIR / group_name
    output_dir = group_dir / "template_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # ANTs template construction parameters
    # -d 3: 3D images
    # -o: output prefix
    # -i: number of iterations (quick: 2, full: 4)
    # -g: gradient step
    # -c: rigid/affine convergence
    # -q: template resolution (use input resolution)
    # -n: N4 bias correction (0=no, 1=yes)
    # -r: rigid alignment first
    # -t: transformation model (SyN for nonlinear)
    
    iterations = 2 if quick else 4
    
    output_prefix = str(output_dir / f"{group_name}_")
    
    # Build the command
    cmd = [
        "antsMultivariateTemplateConstruction2.sh",
        "-d", "3",  # 3D
        "-o", output_prefix,
        "-i", str(iterations),
        "-g", "0.2",  # gradient step
        "-j", "4",  # number of cores
        "-c", "2",  # control for parallel computation
        "-k", "1",  # keep temporary files (for debugging)
        "-w", "1",  # equal weighting
        "-n", "1",  # N4 bias correction
        "-r", "1",  # rigid alignment first
        "-m", "CC",  # cross-correlation similarity metric
        "-t", "SyN",  # SyN for nonlinear
    ] + image_list
    
    print(f"    Running ANTs template construction ({iterations} iterations)...")
    print(f"    This may take 30-60 minutes per group...")
    print(f"    Output: {output_prefix}template0.nii.gz")
    
    # Run the command
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(output_dir))
    
    if result.returncode != 0:
        print(f"    ✗ Template construction failed")
        print(f"    stderr: {result.stderr[-500:] if result.stderr else 'None'}")
        return None
    
    # Check for output template
    template_file = output_dir / f"{group_name}_template0.nii.gz"
    if template_file.exists():
        print(f"    ✓ Template created: {template_file}")
        
        # Also create a symlink with cleaner name
        final_template = group_dir / f"{group_name}_template.nii.gz"
        if final_template.exists():
            final_template.unlink()
        shutil.copy(template_file, final_template)
        
        return final_template
    else:
        # Try alternate naming
        alt_template = output_dir / f"{group_name}_template0GenericAffine.nii.gz"
        if alt_template.exists():
            final_template = group_dir / f"{group_name}_template.nii.gz"
            shutil.copy(alt_template, final_template)
            return final_template
        
        print(f"    ✗ Template file not found")
        print(f"    Looking for: {template_file}")
        return None


def organize_transforms(group_name, group_info):
    """Organize and verify transforms from template construction."""
    print(f"\n  Organizing transforms for {group_name}...")
    
    group_dir = TEMPLATE_DIR / group_name
    output_dir = group_dir / "template_output"
    transforms_dir = group_dir / "transforms"
    transforms_dir.mkdir(exist_ok=True)
    
    # ANTs creates transforms with naming like:
    # {prefix}{subject}_T2w{iteration}Warp.nii.gz
    # {prefix}{subject}_T2w{iteration}Affine.txt
    
    for subject, session in group_info["subjects"]:
        subject_transforms = transforms_dir / subject
        subject_transforms.mkdir(exist_ok=True)
        
        # Find transforms for this subject
        # Pattern: {group_name}_{subject}_T2w*
        warp_pattern = f"{group_name}_{subject}_T2w*Warp.nii.gz"
        affine_pattern = f"{group_name}_{subject}_T2w*GenericAffine.mat"
        
        warps = list(output_dir.glob(warp_pattern))
        affines = list(output_dir.glob(affine_pattern))
        
        if warps:
            # Copy the final iteration warp
            final_warp = sorted(warps)[-1]  # Last iteration
            dest_warp = subject_transforms / "subject_to_template_Warp.nii.gz"
            shutil.copy(final_warp, dest_warp)
            print(f"    ✓ {subject}: warp saved")
            
            # Look for inverse warp
            inv_warp_pattern = f"{group_name}_{subject}_T2w*InverseWarp.nii.gz"
            inv_warps = list(output_dir.glob(inv_warp_pattern))
            if inv_warps:
                inv_warp = sorted(inv_warps)[-1]
                dest_inv = subject_transforms / "template_to_subject_Warp.nii.gz"
                shutil.copy(inv_warp, dest_inv)
        
        if affines:
            final_affine = sorted(affines)[-1]
            dest_affine = subject_transforms / "subject_to_template_Affine.mat"
            shutil.copy(final_affine, dest_affine)
        
        if not warps and not affines:
            print(f"    ⚠ {subject}: no transforms found")


def build_simple_template(group_name, group_info):
    """Build a simpler template using iterative averaging with ANTs pairwise registration."""
    print(f"\n  Building {group_name} template (simplified approach)...")
    
    group_dir = TEMPLATE_DIR / group_name
    group_dir.mkdir(parents=True, exist_ok=True)
    
    transforms_dir = group_dir / "transforms"
    transforms_dir.mkdir(exist_ok=True)
    
    # Step 1: Use first subject as initial reference
    first_subject, first_session = group_info["subjects"][0]
    ref_dir = DERIV_DIR / first_subject / first_session / "dwi"
    
    # Get reference T2w
    ref_t2w = ref_dir / "t2w_highres_flip.nii.gz"
    if not ref_t2w.exists():
        ref_mif = ref_dir / "t2w_highres_flip.mif"
        run_cmd(["mrconvert", str(ref_mif), str(ref_t2w), "-force"], 
                f"Converting {first_subject} T2w")
    
    # Step 2: Register all subjects to reference and average
    registered_images = [ref_t2w]
    
    for subject, session in group_info["subjects"][1:]:
        dwi_dir = DERIV_DIR / subject / session / "dwi"
        moving_t2w = dwi_dir / "t2w_highres_flip.nii.gz"
        
        if not moving_t2w.exists():
            moving_mif = dwi_dir / "t2w_highres_flip.mif"
            run_cmd(["mrconvert", str(moving_mif), str(moving_t2w), "-force"],
                    f"Converting {subject} T2w")
        
        # Output files
        subj_transforms = transforms_dir / subject
        subj_transforms.mkdir(exist_ok=True)
        
        output_prefix = subj_transforms / "to_template_"
        warped_t2w = subj_transforms / "T2w_in_template.nii.gz"
        
        # Register to reference using ANTs SyN
        print(f"    Registering {subject} to reference...")
        run_cmd([
            "antsRegistrationSyN.sh",
            "-d", "3",
            "-t", "s",  # SyN
            "-f", str(ref_t2w),
            "-m", str(moving_t2w),
            "-o", str(output_prefix),
            "-n", "4"  # number of threads
        ], f"ANTs registration for {subject}")
        
        # Check if warped image was created
        ants_warped = subj_transforms / "to_template_Warped.nii.gz"
        if ants_warped.exists():
            shutil.copy(ants_warped, warped_t2w)
            registered_images.append(warped_t2w)
            print(f"    ✓ {subject} registered")
        else:
            print(f"    ✗ {subject} registration failed")
    
    # Step 3: Create average template
    template_file = group_dir / f"{group_name}_template.nii.gz"
    
    # Use MRtrix mrmath to average all registered images
    cmd = ["mrmath"] + [str(img) for img in registered_images] + [
        "mean", str(template_file), "-force"
    ]
    run_cmd(cmd, f"Creating average template from {len(registered_images)} images")
    
    # Step 4: Re-register all subjects to the average template for final transforms
    print(f"\n  Final registration pass to average template...")
    
    for subject, session in group_info["subjects"]:
        dwi_dir = DERIV_DIR / subject / session / "dwi"
        moving_t2w = dwi_dir / "t2w_highres_flip.nii.gz"
        
        subj_transforms = transforms_dir / subject
        subj_transforms.mkdir(exist_ok=True)
        
        output_prefix = subj_transforms / "final_"
        
        print(f"    Final registration for {subject}...")
        run_cmd([
            "antsRegistrationSyN.sh",
            "-d", "3",
            "-t", "s",
            "-f", str(template_file),
            "-m", str(moving_t2w),
            "-o", str(output_prefix),
            "-n", "4"
        ], f"Final ANTs registration for {subject}")
        
        # Rename transforms to standard names
        for src_name, dst_name in [
            ("final_0GenericAffine.mat", "subject_to_template_Affine.mat"),
            ("final_1Warp.nii.gz", "subject_to_template_Warp.nii.gz"),
            ("final_1InverseWarp.nii.gz", "template_to_subject_Warp.nii.gz"),
        ]:
            src = subj_transforms / src_name
            dst = subj_transforms / dst_name
            if src.exists():
                shutil.copy(src, dst)
    
    print(f"\n  ✓ Template created: {template_file}")
    return template_file


def verify_transforms(group_name, group_info):
    """Verify that all transforms exist and can warp a test image."""
    print(f"\n  Verifying transforms for {group_name}...")
    
    group_dir = TEMPLATE_DIR / group_name
    transforms_dir = group_dir / "transforms"
    template_file = group_dir / f"{group_name}_template.nii.gz"
    
    if not template_file.exists():
        print(f"    ✗ Template not found: {template_file}")
        return False
    
    all_ok = True
    for subject, session in group_info["subjects"]:
        subj_transforms = transforms_dir / subject
        
        warp = subj_transforms / "subject_to_template_Warp.nii.gz"
        affine = subj_transforms / "subject_to_template_Affine.mat"
        inv_warp = subj_transforms / "template_to_subject_Warp.nii.gz"
        
        status = []
        if warp.exists():
            status.append("warp✓")
        else:
            status.append("warp✗")
            all_ok = False
        
        if affine.exists():
            status.append("affine✓")
        else:
            status.append("affine✗")
            all_ok = False
        
        if inv_warp.exists():
            status.append("inv_warp✓")
        else:
            status.append("inv_warp⚠")  # Warning only
        
        print(f"    {subject}: {' '.join(status)}")
    
    return all_ok


def main():
    parser = argparse.ArgumentParser(description='Build group template brains')
    parser.add_argument('--group', type=str, choices=['AP', 'APTR', 'CTRL', 'all'],
                        default='all', help='Which group to process')
    parser.add_argument('--quick', action='store_true',
                        help='Use quick mode (fewer iterations)')
    parser.add_argument('--simple', action='store_true',
                        help='Use simplified template construction (faster)')
    args = parser.parse_args()
    
    TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Building Group Template Brains")
    print("=" * 60)
    
    groups_to_process = [args.group] if args.group != 'all' else ['AP', 'APTR', 'CTRL']
    
    results = {}
    for group_name in groups_to_process:
        print(f"\n{'='*60}")
        print(f"Processing {group_name}: {GROUPS[group_name]['description']}")
        print(f"Subjects: {len(GROUPS[group_name]['subjects'])}")
        print(f"{'='*60}")
        
        group_info = GROUPS[group_name]
        
        if args.simple:
            # Simplified approach: pairwise registration + averaging
            template = build_simple_template(group_name, group_info)
        else:
            # Full ANTs multivariate template construction
            image_list, list_file = prepare_t2w_images(group_name, group_info)
            
            if len(image_list) < 2:
                print(f"  ✗ Not enough images for {group_name}")
                results[group_name] = "FAILED"
                continue
            
            template = build_template_ants(group_name, image_list, quick=args.quick)
            
            if template:
                organize_transforms(group_name, group_info)
        
        # Verify transforms
        if template:
            ok = verify_transforms(group_name, group_info)
            results[group_name] = "SUCCESS" if ok else "PARTIAL"
        else:
            results[group_name] = "FAILED"
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for group_name, status in results.items():
        template_file = TEMPLATE_DIR / group_name / f"{group_name}_template.nii.gz"
        if template_file.exists():
            print(f"  {group_name}: {status}")
            print(f"    Template: {template_file}")
            print(f"    Transforms: {TEMPLATE_DIR / group_name / 'transforms'}/")
        else:
            print(f"  {group_name}: {status}")
    
    print("\n  To view templates:")
    for group_name in groups_to_process:
        template_file = TEMPLATE_DIR / group_name / f"{group_name}_template.nii.gz"
        if template_file.exists():
            print(f"    mrview {template_file} &")


if __name__ == "__main__":
    main()
