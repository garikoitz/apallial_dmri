#!/usr/bin/env python3
"""
Regenerate tract files that were created with old minlength parameter.
Uses correct parameters from config.yaml.
"""

import subprocess
import sys
from pathlib import Path
import yaml

def load_config():
    config_path = Path(__file__).parent / "config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def get_tracts_to_regenerate():
    """Find all tracts with minlength 0.5 that need regeneration."""
    cfg = load_config()
    deriv_dir = Path(cfg['paths']['derivatives_dir'])
    
    tracts = []
    for tck in deriv_dir.rglob("tracks_Isocortex_*.tck"):
        # Skip derived files
        if any(x in tck.name for x in ['_cleaned', '_sift', '_top', '_SF', '_warped']):
            continue
        # Skip tracts subfolder (old format)
        if '/tracts/' in str(tck):
            continue
        
        # Check if this tract has minlength 0.5
        result = subprocess.run(['tckinfo', str(tck)], capture_output=True, text=True)
        if 'minlength 0.5' in result.stdout:
            tracts.append(tck)
    
    return sorted(tracts)

def regenerate_tract(tck_file, cfg, dry_run=False):
    """Regenerate a single tract file with correct parameters."""
    dwi_dir = tck_file.parent
    
    # Extract tract info from filename: tracks_Isocortex_Striatum_left.tck
    name = tck_file.stem.replace('tracks_', '')  # Isocortex_Striatum_left
    parts = name.split('_')  # ['Isocortex', 'Striatum', 'left']
    
    if len(parts) < 3:
        return False, f"Cannot parse tract name: {name}"
    
    seed_region = parts[0]  # Isocortex
    target_region = parts[1]  # Striatum or Thalamus
    side = parts[2]  # left or right
    
    # Find required files
    fod_file = dwi_dir / "fod_wm.mif"
    seed_roi = dwi_dir / "atlas" / "rois" / f"roi_{seed_region.lower()}_{side}.nii.gz"
    target_roi = dwi_dir / "atlas" / "rois" / f"roi_{target_region.lower()}_{side}.nii.gz"
    
    # Check files exist
    if not fod_file.exists():
        return False, f"FOD file not found: {fod_file}"
    if not seed_roi.exists():
        return False, f"Seed ROI not found: {seed_roi}"
    if not target_roi.exists():
        return False, f"Target ROI not found: {target_roi}"
    
    # Get tractography parameters from config
    tract_cfg = cfg.get('tractography', {})
    algorithm = tract_cfg.get('algorithm', 'iFOD2')
    select = tract_cfg.get('select', 5000)
    step_size = tract_cfg.get('step_size', 0.02)
    angle = tract_cfg.get('angle', 60)
    min_length = tract_cfg.get('min_length', 1.5)
    max_length = tract_cfg.get('max_length', 15)
    cutoff = tract_cfg.get('cutoff', 0.05)
    
    # Build tckgen command
    cmd = [
        'tckgen', str(fod_file), str(tck_file),
        '-algorithm', algorithm,
        '-seed_image', str(seed_roi),
        '-include', str(target_roi),
        '-select', str(select),
        '-step', str(step_size),
        '-angle', str(angle),
        '-minlength', str(min_length),
        '-maxlength', str(max_length),
        '-cutoff', str(cutoff),
        '-force'
    ]
    
    if dry_run:
        return True, f"Would run: tckgen ... -minlength {min_length} -force"
    
    # Backup old file
    backup = tck_file.with_suffix('.tck.old')
    if tck_file.exists():
        tck_file.rename(backup)
    
    # Run tckgen
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        # Restore backup
        if backup.exists():
            backup.rename(tck_file)
        return False, f"tckgen failed: {result.stderr[:200]}"
    
    # Remove backup
    if backup.exists():
        backup.unlink()
    
    # Also remove old cleaned files that are now invalid
    for suffix in ['_cleaned.tck', '_SF.tck', '_cleaned_top200.tck', '_cleaned_sift2_weights.txt']:
        old_file = dwi_dir / (tck_file.stem + suffix)
        if old_file.exists():
            old_file.unlink()
    
    return True, f"Regenerated with minlength={min_length}"

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Regenerate tracts with correct minlength')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done')
    args = parser.parse_args()
    
    cfg = load_config()
    min_length = cfg.get('tractography', {}).get('min_length', 1.5)
    
    print("=" * 60)
    print("TRACT REGENERATION - Fix minlength parameter")
    print("=" * 60)
    print(f"Target minlength: {min_length} mm")
    print()
    
    tracts = get_tracts_to_regenerate()
    print(f"Found {len(tracts)} tracts to regenerate")
    print()
    
    if args.dry_run:
        print("[DRY RUN MODE]")
        print()
    
    success = 0
    for i, tck in enumerate(tracts, 1):
        # Extract subject from path
        subject = "unknown"
        for part in tck.parts:
            if part.startswith('sub-'):
                subject = part
                break
        
        tract_name = tck.stem.replace('tracks_', '')
        print(f"[{i}/{len(tracts)}] {subject}: {tract_name}...", end=" ", flush=True)
        
        ok, msg = regenerate_tract(tck, cfg, dry_run=args.dry_run)
        
        if ok:
            success += 1
            print(f"✓ {msg}")
        else:
            print(f"✗ {msg}")
    
    print()
    print("=" * 60)
    print(f"Completed: {success}/{len(tracts)} tracts regenerated")
    print("=" * 60)
    
    if success < len(tracts):
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
