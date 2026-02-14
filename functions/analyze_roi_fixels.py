#!/usr/bin/env python3
"""
ROI-Based Fixel Analysis

Extract fixels from group-calibrated FODs and compute directional metrics within ROIs.

This script:
1. Runs fod2fixel on each subject's wm_fod_norm_group.mif
2. Loads fixel data (FD, directions) and ROI masks
3. Computes directional FD metrics: FD_total, FD_ap, FD_dv, FD_lm
4. Generates CSV reports and plots

Usage:
    python analyze_roi_fixels.py --subject sub-Apallial17NmF --roi graft_left
    python analyze_roi_fixels.py --group APTR --roi graft_left
    python analyze_roi_fixels.py --group AP --roi graft_left --output reports/fixel_analysis.csv
"""

import os
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import numpy as np
import nibabel as nib

# Set pipeline
os.environ['MBRAIN_PIPELINE'] = 'flip1'

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logger import setup_logger
from functions.config_loader import load_config

# Global logger
logger = None


def run_cmd(cmd: str, description: str = None, check: bool = True) -> bool:
    """Run shell command with logging."""
    if description:
        logger.info(f"  {description}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        if check:
            logger.error(f"  ✗ FAILED: {result.stderr[:500]}")
            return False
        else:
            logger.warning(f"  ⚠ Warning: {result.stderr[:200]}")
    
    logger.info(f"  ✓ Success")
    return True


def generate_fixels(fod_file: Path, output_dir: Path, force: bool = False) -> Optional[Path]:
    """
    Generate fixels from FOD using fod2fixel.
    
    Returns path to fixel directory if successful.
    """
    if output_dir.exists() and not force:
        logger.info(f"  Fixels already exist: {output_dir}")
        return output_dir
    
    # fod2fixel command - output AFD and peak amplitude
    cmd = f"fod2fixel '{fod_file}' '{output_dir}' -mask '{fod_file.parent / 'mask_anat.mif'}' -afd afd.mif -peak peak.mif -force -quiet"
    
    if not run_cmd(cmd, f"Generating fixels from {fod_file.name}..."):
        return None
    
    return output_dir


def load_mif_as_nifti(mif_path: Path) -> np.ndarray:
    """Convert MIF to NIfTI temporarily and load data."""
    with tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False) as tmp:
        tmp_path = tmp.name
    
    subprocess.run(['mrconvert', str(mif_path), tmp_path, '-force', '-quiet'], check=True)
    img = nib.load(tmp_path)
    data = img.get_fdata()
    os.remove(tmp_path)
    
    return data


def load_fixel_data(fixel_dir: Path) -> Dict[str, np.ndarray]:
    """
    Load fixel data from fod2fixel output directory.
    
    Returns dict with:
        'index': voxel-to-fixel mapping (X, Y, Z, 2) where [x,y,z,0]=count, [x,y,z,1]=offset
        'directions': unit direction vectors (N x 3 x 1)
        'fd': fiber density values (N,) - uses afd.mif if available, else fd.mif
    """
    index_file = fixel_dir / 'index.mif'
    directions_file = fixel_dir / 'directions.mif'
    
    # Prefer afd.mif (apparent fiber density) over fd.mif
    afd_file = fixel_dir / 'afd.mif'
    fd_file = fixel_dir / 'fd.mif'
    density_file = afd_file if afd_file.exists() else fd_file
    
    # Load index (maps voxels to fixels)
    index_data = load_mif_as_nifti(index_file)
    
    # Load directions (Nx3 matrix)
    directions_data = load_mif_as_nifti(directions_file)
    
    # Load FD values (N vector)
    fd_data = load_mif_as_nifti(density_file).flatten()
    
    return {
        'index': index_data,
        'directions': directions_data,
        'fd': fd_data
    }


def load_roi_mask(roi_path: Path) -> np.ndarray:
    """Load ROI mask (handles both .mif and .nii.gz)."""
    if roi_path.suffix == '.mif':
        return load_mif_as_nifti(roi_path) > 0.5
    else:
        img = nib.load(roi_path)
        return img.get_fdata() > 0.5


def compute_roi_fixel_metrics(fixel_data: Dict, roi_mask: np.ndarray) -> Dict:
    """
    Compute directional FD metrics for fixels within ROI.
    
    Uses the sparse fixel index format from fod2fixel:
        index[x,y,z,0] = number of fixels in this voxel
        index[x,y,z,1] = offset into the flat fixel arrays
    
    Axis convention (verified from FOD transform, near-identity):
        dir[0] = AP (Anterior-Posterior, axis 0)
        dir[1] = DV (Dorsal-Ventral, axis 1)
        dir[2] = LR (Left-Right, axis 2)
    
    Returns:
        fd_total: Total fiber density (sum of AFD across all fixels in ROI)
        fd_ap: FD weighted by |dir[0]| (AP component)
        fd_dv: FD weighted by |dir[1]| (DV component)
        fd_lm: FD weighted by |dir[2]| (LR/LM component)
        n_fixels: Number of fixels in ROI
        n_voxels: Number of ROI voxels with at least one fixel
    """
    index = fixel_data['index']       # shape: (X, Y, Z, 2)
    directions = fixel_data['directions']  # shape: (N, 3, 1)
    fd = fixel_data['fd']             # shape: (N,) flat
    
    # Squeeze directions to (N, 3)
    dirs = directions[:, :, 0] if directions.ndim == 3 else directions
    
    # Find voxels in ROI
    roi_coords = np.where(roi_mask)
    
    fd_total = 0.0
    fd_ap = 0.0
    fd_dv = 0.0
    fd_lm = 0.0
    n_fixels = 0
    n_voxels_with_fixels = 0
    
    for x, y, z in zip(roi_coords[0], roi_coords[1], roi_coords[2]):
        count = int(index[x, y, z, 0])
        offset = int(index[x, y, z, 1])
        
        if count == 0:
            continue
        
        n_voxels_with_fixels += 1
        
        for i in range(count):
            fi = offset + i
            if fi >= len(fd):
                break
            
            fd_value = fd[fi]
            direction = dirs[fi]
            
            if fd_value < 0.001:
                continue
            
            fd_total += fd_value
            fd_ap += fd_value * np.abs(direction[0])  # axis 0 = AP
            fd_dv += fd_value * np.abs(direction[1])  # axis 1 = DV
            fd_lm += fd_value * np.abs(direction[2])  # axis 2 = LR
            n_fixels += 1
    
    if logger:
        logger.info(f"  ROI voxels: {len(roi_coords[0])}, with fixels: {n_voxels_with_fixels}, total fixels: {n_fixels}")
    
    return {
        'fd_total': fd_total,
        'fd_ap': fd_ap,
        'fd_dv': fd_dv,
        'fd_lm': fd_lm,
        'n_fixels': n_fixels,
        'n_voxels': n_voxels_with_fixels
    }


def process_subject(subject: str, session: str, group: str, 
                   roi_name: str, cfg: dict, force: bool = False) -> Optional[Dict]:
    """
    Process one subject: generate fixels and compute ROI metrics.
    
    Returns: Dict with metrics or None if failed
    """
    logger.info(f"\n--- {subject} ({group}) ---")
    
    deriv_dir = Path(cfg['paths']['derivatives_flip1'])
    subj_dir = deriv_dir / subject / session
    
    # Input files
    fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'
    mask_file = subj_dir / 'dwi' / 'mask_anat.mif'
    
    # Check if group FOD exists
    if not fod_file.exists():
        logger.warning(f"  ⚠ Skipping: Group FOD not found")
        return None
    
    if not mask_file.exists():
        logger.warning(f"  ⚠ Skipping: Mask not found")
        return None
    
    # Generate fixels
    fixel_dir = subj_dir / 'fixels_group'
    fixel_dir = generate_fixels(fod_file, fixel_dir, force)
    
    if not fixel_dir:
        logger.error(f"  ✗ Failed to generate fixels")
        return None
    
    # Load fixel data
    logger.info("  Loading fixel data...")
    fixel_data = load_fixel_data(fixel_dir)
    
    # Load ROI mask
    roi_path = find_roi_path(subj_dir, roi_name, group)
    if not roi_path:
        logger.warning(f"  ⚠ ROI {roi_name} not found")
        return None
    
    logger.info(f"  Loading ROI: {roi_path.name}")
    roi_mask = load_roi_mask(roi_path)
    
    # Compute metrics
    logger.info("  Computing directional metrics...")
    metrics = compute_roi_fixel_metrics(fixel_data, roi_mask)
    
    # Add subject info
    metrics['subject'] = subject
    metrics['session'] = session
    metrics['group'] = group
    metrics['roi'] = roi_name
    
    logger.info(f"  ✓ FD_total={metrics['fd_total']:.2f}, n_fixels={metrics['n_fixels']}")
    
    return metrics


def find_roi_path(subj_dir: Path, roi_name: str, group: str) -> Optional[Path]:
    """
    Find ROI file path based on naming convention.
    
    roi_name examples:
        - 'graft_left' -> manual_roi/fod_space/roi_Graft_left_flip_fod.mif
        - 'isocortex_left' -> atlas/rois/fod_space/roi_isocortex_left.nii.gz
    """
    if 'graft' in roi_name.lower():
        # Manual ROI
        hemi = roi_name.split('_')[-1]  # left or right
        roi_path = subj_dir / 'manual_roi' / 'fod_space' / f'roi_Graft_{hemi}_flip_fod.mif'
    elif 'isocortex' in roi_name.lower():
        # Atlas ROI
        hemi = roi_name.split('_')[-1]
        roi_path = subj_dir / 'atlas' / 'rois' / 'fod_space' / f'roi_isocortex_{hemi}.nii.gz'
    else:
        logger.warning(f"Unknown ROI naming pattern: {roi_name}")
        return None
    
    if roi_path.exists():
        return roi_path
    
    return None


def get_subjects(cfg: dict, group: Optional[str] = None, subject: Optional[str] = None) -> List[Tuple[str, str, str]]:
    """Get list of (subject_id, session, group) tuples."""
    groups = cfg.get('groups', {})
    subjects = []
    
    if subject:
        for g_name, subj_list in groups.items():
            for s in subj_list:
                if s['id'] == subject:
                    subjects.append((s['id'], s['session'], g_name))
                    return subjects
        raise ValueError(f"Subject {subject} not found")
    
    if group:
        if group not in groups:
            raise ValueError(f"Group {group} not found")
        for s in groups[group]:
            subjects.append((s['id'], s['session'], group))
    else:
        # All groups
        for g_name, subj_list in groups.items():
            for s in subj_list:
                subjects.append((s['id'], s['session'], g_name))
    
    return subjects


def save_results_csv(results: List[Dict], output_file: Path):
    """Save results to CSV."""
    import csv
    
    if not results:
        logger.warning("No results to save")
        return
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = ['subject', 'session', 'group', 'roi', 'n_fixels', 
                  'fd_total', 'fd_ap', 'fd_dv', 'fd_lm',
                  'fd_ap_frac', 'fd_dv_frac', 'fd_lm_frac']
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for r in results:
            # Compute fractions
            total = r['fd_ap'] + r['fd_dv'] + r['fd_lm']
            row = {
                'subject': r['subject'],
                'session': r['session'],
                'group': r['group'],
                'roi': r['roi'],
                'n_fixels': r['n_fixels'],
                'fd_total': f"{r['fd_total']:.4f}",
                'fd_ap': f"{r['fd_ap']:.4f}",
                'fd_dv': f"{r['fd_dv']:.4f}",
                'fd_lm': f"{r['fd_lm']:.4f}",
                'fd_ap_frac': f"{r['fd_ap']/total:.4f}" if total > 0 else '0',
                'fd_dv_frac': f"{r['fd_dv']/total:.4f}" if total > 0 else '0',
                'fd_lm_frac': f"{r['fd_lm']/total:.4f}" if total > 0 else '0',
            }
            writer.writerow(row)
    
    logger.info(f"\n✓ Results saved to: {output_file}")


def analyze_probe_fixels(cfg: dict, output_dir: str, force: bool = False) -> Optional[str]:
    """
    Fixel-based probe comparison: same 3 subjects and single-slice probe ROIs
    as fod_direction_analysis.py --probe, but using proper fixel analysis.
    
    Generates radar chart: reports/fixel_probe_comparison.png
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    deriv = cfg['paths']['derivatives_flip1']
    
    # Same probe subjects/ROIs as FOD version
    probe_subjects = {
        'CTRL': {
            'subject': 'sub-129-921-R1F',
            'session': 'ses-260115121231',
            'roi_path': os.path.join(deriv, 'sub-129-921-R1F', 'ses-260115121231',
                                      'atlas', 'rois', 'fod_space',
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        },
        'AP': {
            'subject': 'sub-769-C-O5-R1M',
            'session': 'ses-240911115733',
            'roi_path': os.path.join(deriv, 'sub-769-C-O5-R1M', 'ses-240911115733',
                                      'manual_roi', 'fod_space',
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        },
        'APTR': {
            'subject': 'sub-AP72-862-901-NmF',
            'session': 'ses-251210093342',
            'roi_path': os.path.join(deriv, 'sub-AP72-862-901-NmF', 'ses-251210093342',
                                      'manual_roi', 'fod_space',
                                      'ROI_probe_square_left_1_coronal_slice_fod.mif')
        }
    }
    
    group_order = ['CTRL', 'AP', 'APTR']
    group_colors = {'CTRL': '#2196F3', 'AP': '#FF9800', 'APTR': '#4CAF50'}
    results = {}
    
    for group in group_order:
        info = probe_subjects[group]
        logger.info(f"\n--- {group}: {info['subject']} ---")
        
        # Check ROI exists
        roi_path = Path(info['roi_path'])
        if not roi_path.exists():
            logger.error(f"  Probe ROI not found: {roi_path}")
            continue
        
        # Paths
        subj_dir = Path(deriv) / info['subject'] / info['session']
        fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'
        
        if not fod_file.exists():
            logger.error(f"  FOD not found: {fod_file}")
            continue
        
        # Generate fixels
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        fixel_dir = generate_fixels(fod_file, fixel_dir, force)
        if not fixel_dir:
            continue
        
        # Load data
        logger.info("  Loading fixel data...")
        fixel_data = load_fixel_data(fixel_dir)
        
        logger.info(f"  Loading probe ROI: {roi_path.name}")
        roi_mask = load_roi_mask(roi_path)
        
        # Compute metrics
        metrics = compute_roi_fixel_metrics(fixel_data, roi_mask)
        
        if metrics['n_fixels'] == 0:
            logger.warning(f"  No fixels found in ROI!")
            continue
        
        results[group] = {
            'data': metrics,
            'subject': info['subject']
        }
        
        total_dir = metrics['fd_ap'] + metrics['fd_dv'] + metrics['fd_lm']
        logger.info(f"  FD_total={metrics['fd_total']:.3f}, fixels={metrics['n_fixels']}, "
                    f"AP={metrics['fd_ap']/total_dir*100:.1f}%, "
                    f"DV={metrics['fd_dv']/total_dir*100:.1f}%, "
                    f"LR={metrics['fd_lm']/total_dir*100:.1f}%")
    
    if not results:
        logger.error("No results collected!")
        return None
    
    # --- Generate radar comparison figure ---
    return plot_fixel_probe_radar(results, group_order, group_colors, output_dir)


def plot_fixel_probe_radar(results: Dict, group_order: list, group_colors: dict, 
                           output_dir: str) -> str:
    """
    Generate fixel probe comparison radar figure.
    
    Layout: one radar per group + one overlaid summary.
    Each radar: 3 axes (AP, DV, LR) showing AFD magnitude.
    Circle size = total FD (proxy for overall fiber content).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    
    categories = ['DV', 'AP', 'ML']
    n_cats = len(categories)
    angles = np.linspace(0, 2 * np.pi, n_cats, endpoint=False).tolist()
    angles += angles[:1]  # close polygon
    
    n_groups = len(results)
    n_cols = n_groups + 1  # individual + overlay
    
    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 5),
                             subplot_kw=dict(projection='polar'))
    
    # Find global max for normalization
    max_amplitude = 0
    max_fd_total = 0
    for group in group_order:
        if group not in results:
            continue
        d = results[group]['data']
        max_amplitude = max(max_amplitude, d['fd_ap'], d['fd_dv'], d['fd_lm'])
        max_fd_total = max(max_fd_total, d['fd_total'])
    
    # Individual radars
    for idx, group in enumerate(group_order):
        if group not in results:
            continue
        
        ax = axes[idx]
        d = results[group]['data']
        color = group_colors[group]
        
        # Normalize to global max — order: DV, AP, ML
        values = [d['fd_dv'] / max_amplitude,
                  d['fd_ap'] / max_amplitude,
                  d['fd_lm'] / max_amplitude]
        values += values[:1]
        
        ax.plot(angles, values, 'o-', linewidth=2.5, color=color, markersize=8)
        ax.fill(angles, values, alpha=0.25, color=color)
        
        # WM content circle (proportional to total FD)
        fd_ratio = d['fd_total'] / max_fd_total if max_fd_total > 0 else 0
        circle = Circle((0, 0), fd_ratio * 0.15, transform=ax.transData,
                        color=color, alpha=0.4, zorder=5)
        ax.add_patch(circle)
        
        ax.set_theta_zero_location("S")  # DV at bottom (South)
        ax.set_theta_direction(1)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, fontsize=11, fontweight='bold')
        ax.set_ylim(0, 1.1)
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['25%', '50%', '75%', '100%'], fontsize=7)
        ax.grid(True, alpha=0.3)
        
        subj_short = results[group]['subject'].replace('sub-', '')[:18]
        total_dir = d['fd_ap'] + d['fd_dv'] + d['fd_lm']
        ax.set_title(f"{group}\n{subj_short}\n"
                     f"FD={d['fd_total']:.2f}, fixels={d['n_fixels']}\n"
                     f"AP={d['fd_ap']/total_dir*100:.0f}% DV={d['fd_dv']/total_dir*100:.0f}% ML={d['fd_lm']/total_dir*100:.0f}%",
                     fontsize=10, fontweight='bold', pad=20)
    
    # Overlay plot
    ax_summary = axes[n_groups]
    for group in group_order:
        if group not in results:
            continue
        d = results[group]['data']
        color = group_colors[group]
        
        values = [d['fd_dv'] / max_amplitude,
                  d['fd_ap'] / max_amplitude,
                  d['fd_lm'] / max_amplitude]
        values += values[:1]
        
        ax_summary.plot(angles, values, 'o-', linewidth=3, color=color,
                        markersize=10, label=group)
        ax_summary.fill(angles, values, alpha=0.15, color=color)
    
    ax_summary.set_theta_zero_location("S")  # DV at bottom (South)
    ax_summary.set_theta_direction(1)
    ax_summary.set_xticks(angles[:-1])
    ax_summary.set_xticklabels(categories, fontsize=11, fontweight='bold')
    ax_summary.set_ylim(0, 1.1)
    ax_summary.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax_summary.set_yticklabels(['25%', '50%', '75%', '100%'], fontsize=7)
    ax_summary.grid(True, alpha=0.3)
    ax_summary.set_title('All Groups\nOverlaid', fontsize=11, fontweight='bold', pad=20)
    ax_summary.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0), fontsize=10)
    
    fig.suptitle('Fixel-Based Directional Analysis - Single-Slice Probe ROIs\n'
                 '(Radar: directional AFD, normalized to max across groups; '
                 'center circle = total FD)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'fixel_probe_comparison.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    logger.info(f"\nSaved: {output_path}")
    
    # Print summary table
    logger.info("\n" + "=" * 70)
    logger.info("FIXEL PROBE ANALYSIS SUMMARY")
    logger.info("=" * 70)
    logger.info(f"{'Group':>8} {'Subject':>20} {'Fixels':>8} {'FD_tot':>8} {'AP%':>7} {'DV%':>7} {'ML%':>7}")
    logger.info("-" * 70)
    
    for group in group_order:
        if group not in results:
            continue
        d = results[group]['data']
        subj = results[group]['subject'].replace('sub-', '')[:18]
        total_dir = d['fd_ap'] + d['fd_dv'] + d['fd_lm']
        logger.info(f"{group:>8} {subj:>20} {d['n_fixels']:>8} {d['fd_total']:>8.2f} "
                    f"{d['fd_ap']/total_dir*100:>6.1f}% {d['fd_dv']/total_dir*100:>6.1f}% "
                    f"{d['fd_lm']/total_dir*100:>6.1f}%")
    
    logger.info("=" * 70)
    logger.info("FD = Fiber Density (integral of FOD lobes from fod2fixel)")
    logger.info("AP/DV/ML = directional components weighted by fixel orientation")
    
    return output_path


def compute_slicewise_fixel_metrics(fixel_data: Dict, roi_mask: np.ndarray, 
                                     slice_axis: int = 0) -> Dict[int, Dict]:
    """
    Compute fixel metrics per slice along given axis within an ROI.
    
    Returns dict mapping slice_index -> metrics dict (same keys as compute_roi_fixel_metrics).
    Only slices with nonzero ROI voxels are included.
    """
    index = fixel_data['index']
    directions = fixel_data['directions']
    fd = fixel_data['fd']
    dirs = directions[:, :, 0] if directions.ndim == 3 else directions
    
    n_slices = roi_mask.shape[slice_axis]
    slice_data = {}
    
    for s in range(n_slices):
        # Extract single slice of ROI
        slicer = [slice(None)] * 3
        slicer[slice_axis] = s
        roi_slice = roi_mask[tuple(slicer)]
        
        if roi_slice.sum() == 0:
            continue
        
        # Build full 3D mask for this single slice
        mask_3d = np.zeros_like(roi_mask, dtype=bool)
        mask_3d[tuple(slicer)] = roi_slice
        
        # Get fixel metrics for this slice
        roi_coords = np.where(mask_3d)
        fd_total = fd_ap = fd_dv = fd_lm = 0.0
        n_fixels = 0
        n_voxels_with_fixels = 0
        
        for x, y, z in zip(roi_coords[0], roi_coords[1], roi_coords[2]):
            count = int(index[x, y, z, 0])
            offset = int(index[x, y, z, 1])
            if count == 0:
                continue
            n_voxels_with_fixels += 1
            for i in range(count):
                fi = offset + i
                if fi >= len(fd):
                    break
                fd_value = fd[fi]
                direction = dirs[fi]
                if fd_value < 0.001:
                    continue
                fd_total += fd_value
                fd_ap += fd_value * np.abs(direction[0])
                fd_dv += fd_value * np.abs(direction[1])
                fd_lm += fd_value * np.abs(direction[2])
                n_fixels += 1
        
        if n_fixels > 0:
            slice_data[s] = {
                'fd_total': fd_total,
                'fd_ap': fd_ap,
                'fd_dv': fd_dv,
                'fd_lm': fd_lm,
                'n_fixels': n_fixels,
                'n_voxels': n_voxels_with_fixels,
                'n_roi_voxels': int(roi_slice.sum())
            }
    
    return slice_data


def analyze_fixel_3group_comparison(hemisphere: str, cfg: dict, output_dir: str,
                                     use_shrink: bool = False, force: bool = False) -> Optional[str]:
    """
    Fixel-based 3-group comparison figure (CTRL vs AP vs APTR).
    
    Each row = one group. Columns = coronal slices (from anterior).
    Each cell = radar plot showing directional AFD (DV/AP/ML).
    Center circle = total FD (proxy for fiber content).
    Extra columns: per-subject average, group average, info panel.
    
    For CTRL: uses isocortex ROI from atlas/rois/fod_space/
    For AP/APTR: uses graft ROI from manual_roi/fod_space/ (optionally shrunk)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    shrink_label = " (SHRUNK ROIs)" if use_shrink else ""
    logger.info(f"\n{'='*60}")
    logger.info(f"FIXEL 3-GROUP COMPARISON - {hemisphere.upper()} hemisphere{shrink_label}")
    logger.info(f"{'='*60}")
    
    deriv = Path(cfg['paths']['derivatives_flip1'])
    
    # Subject definitions
    all_subjects = {
        'CTRL': [
            ('sub-129-921-R1F', 'ses-260115121231'),
            ('sub-129-921-R1L1F', 'ses-260116085433'),
            ('sub-129-1014-L1M', 'ses-260115113814'),
            ('sub-138-1010-L1F', 'ses-260115101441'),
        ],
        'AP': [
            ('sub-596-C-OO-NmM', 'ses-240911103349'),
            ('sub-769-C-O5-R1M', 'ses-240911115733'),
            ('sub-769-C-OO-NmM', 'ses-240911111853'),
            ('sub-875-C-OO-NM', 'ses-240911100253'),
        ],
        'APTR': [
            ('sub-AP72-862-901-NmF', 'ses-251210093342'),
            ('sub-AP72-862-901-R1M', 'ses-251210100043'),
            ('sub-Apallial17NmF', 'ses-240816091633'),
            ('sub-Apallial20R1F', 'ses-240909094232'),
        ],
    }
    
    colors = ['#e74c3c', '#3498db', '#27ae60', '#9b59b6']
    
    # Collect slicewise fixel data for all subjects
    all_group_data = {}
    graft_slice_counts = []
    
    for group in ['AP', 'APTR']:
        logger.info(f"\n--- Collecting {group} data{shrink_label} ---")
        all_group_data[group] = {}
        
        for (subj, ses), color in zip(all_subjects[group], colors):
            logger.info(f"\nProcessing {subj}...")
            
            subj_dir = deriv / subj / ses
            fixel_dir = subj_dir / 'dwi' / 'fixels_group'
            
            # Generate fixels if needed
            fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'
            fixel_dir = generate_fixels(fod_file, fixel_dir, force)
            if not fixel_dir:
                continue
            
            # Load fixel data
            fixel_data = load_fixel_data(fixel_dir)
            
            # Load ROI
            if use_shrink:
                roi_path = subj_dir / 'manual_roi' / 'fod_space' / 'shrink' / f'roi_Graft_{hemisphere}_flip_fod_shrink.mif'
            else:
                roi_path = subj_dir / 'manual_roi' / 'fod_space' / f'roi_Graft_{hemisphere}_flip_fod.mif'
            
            if not roi_path.exists():
                logger.warning(f"  ROI not found: {roi_path}")
                continue
            
            roi_mask = load_roi_mask(roi_path)
            
            # Compute slicewise metrics (axis 0 = AP axis)
            slice_data = compute_slicewise_fixel_metrics(fixel_data, roi_mask, slice_axis=0)
            
            if slice_data:
                all_group_data[group][subj] = {
                    'data': slice_data,
                    'color': color,
                    'session': ses
                }
                graft_slice_counts.append(len(slice_data))
                logger.info(f"  Found {len(slice_data)} slices with fixel data")
            else:
                logger.warning(f"  No fixel data in ROI!")
    
    if not graft_slice_counts:
        logger.error("No graft fixel data found!")
        return None
    
    target_n_slices = min(graft_slice_counts)
    logger.info(f"\nTarget slice count (smallest graft): {target_n_slices}")
    
    # Re-index graft data from anterior (sorted slices, take first target_n_slices)
    for group in ['AP', 'APTR']:
        for subj, subj_info in all_group_data[group].items():
            old_data = subj_info['data']
            sorted_slices = sorted(old_data.keys())[:target_n_slices]
            new_data = {i: old_data[s] for i, s in enumerate(sorted_slices)}
            subj_info['data'] = new_data
    
    # Collect CTRL data (isocortex ROI)
    logger.info(f"\n--- Collecting CTRL data (isocortex ROI) ---")
    all_group_data['CTRL'] = {}
    
    for (subj, ses), color in zip(all_subjects['CTRL'], colors):
        logger.info(f"\nProcessing {subj}...")
        
        subj_dir = deriv / subj / ses
        fixel_dir = subj_dir / 'dwi' / 'fixels_group'
        
        fod_file = subj_dir / 'dwi' / 'wm_fod_norm_group.mif'
        fixel_dir = generate_fixels(fod_file, fixel_dir, force)
        if not fixel_dir:
            continue
        
        fixel_data = load_fixel_data(fixel_dir)
        
        roi_path = subj_dir / 'atlas' / 'rois' / 'fod_space' / f'roi_isocortex_{hemisphere}.nii.gz'
        if not roi_path.exists():
            logger.warning(f"  ROI not found: {roi_path}")
            continue
        
        roi_mask = load_roi_mask(roi_path)
        slice_data = compute_slicewise_fixel_metrics(fixel_data, roi_mask, slice_axis=0)
        
        if slice_data:
            all_group_data['CTRL'][subj] = {
                'data': slice_data,
                'color': color,
                'session': ses
            }
            logger.info(f"  Found {len(slice_data)} slices with fixel data")
    
    # Re-index CTRL to match target slice count
    for subj, subj_info in all_group_data['CTRL'].items():
        old_data = subj_info['data']
        sorted_slices = sorted(old_data.keys())[:target_n_slices]
        new_data = {i: old_data[s] for i, s in enumerate(sorted_slices)}
        subj_info['data'] = new_data
    
    # Find global max FD and directional amplitude for normalization
    global_max_fd = 0
    global_max_amplitude = 0
    for group_data in all_group_data.values():
        for subj_info in group_data.values():
            for sd in subj_info['data'].values():
                global_max_fd = max(global_max_fd, sd['fd_total'])
                max_comp = max(sd['fd_ap'], sd['fd_dv'], sd['fd_lm'])
                global_max_amplitude = max(global_max_amplitude, max_comp)
    
    logger.info(f"\nGlobal max FD: {global_max_fd:.2f}")
    logger.info(f"Global max directional AFD component: {global_max_amplitude:.2f}")
    
    # Calculate average values per slice per group
    avg_fd_per_group = {}
    avg_dir_per_group = {}
    
    for group in ['CTRL', 'AP', 'APTR']:
        group_data = all_group_data[group]
        avg_fd = []
        avg_ap, avg_dv, avg_lm = [], [], []
        
        for slice_idx in range(target_n_slices):
            fd_sum, ap_sum, dv_sum, lm_sum = 0, 0, 0, 0
            count = 0
            for subj_info in group_data.values():
                if slice_idx in subj_info['data']:
                    sd = subj_info['data'][slice_idx]
                    fd_sum += sd['fd_total']
                    ap_sum += sd['fd_ap']
                    dv_sum += sd['fd_dv']
                    lm_sum += sd['fd_lm']
                    count += 1
            
            avg_fd.append(fd_sum / count if count > 0 else 0)
            avg_ap.append(ap_sum / count if count > 0 else 0)
            avg_dv.append(dv_sum / count if count > 0 else 0)
            avg_lm.append(lm_sum / count if count > 0 else 0)
        
        avg_fd_per_group[group] = avg_fd
        avg_dir_per_group[group] = {'ap': avg_ap, 'dv': avg_dv, 'lm': avg_lm}
    
    # Create figure
    n_cols = target_n_slices + 3  # +1 per-subject avg, +1 group avg, +1 info
    n_rows = 3
    
    fig = plt.figure(figsize=(2.2 * n_cols, 3.5 * n_rows))
    
    categories = ['DV', 'AP', 'ML']
    angles = [n / 3.0 * 2 * np.pi for n in range(3)]
    angles += angles[:1]
    
    row_groups = ['CTRL', 'AP', 'APTR']
    group_labels = ['CTRL\n(isocortex)', 'AP\n(graft, no Tx)', 'APTR\n(graft + Tx)']
    group_colors_bg = ['#e8f5e9', '#ffebee', '#e3f2fd']
    
    for row_idx, (group, label, bg_color) in enumerate(zip(row_groups, group_labels, group_colors_bg)):
        group_data = all_group_data[group]
        
        # Individual slices
        for col_idx in range(target_n_slices):
            ax_idx = row_idx * n_cols + col_idx + 1
            ax = fig.add_subplot(n_rows, n_cols, ax_idx, projection='polar')
            ax.set_facecolor(bg_color)
            
            for subj, subj_info in group_data.items():
                data = subj_info['data']
                color = subj_info['color']
                
                if col_idx in data:
                    sd = data[col_idx]
                    values = [sd['fd_dv'] / global_max_amplitude,
                              sd['fd_ap'] / global_max_amplitude,
                              sd['fd_lm'] / global_max_amplitude]
                    values += values[:1]
                    
                    ax.plot(angles, values, 'o-', linewidth=1.5, color=color,
                            markersize=3, alpha=0.7)
                    ax.fill(angles, values, alpha=0.1, color=color)
            
            # Center circle for average FD (scaled to global max)
            avg_fd = avg_fd_per_group[group][col_idx]
            circle_radius = 0.25 * np.sqrt(avg_fd / global_max_fd) if global_max_fd > 0 else 0
            
            if circle_radius > 0.005:
                circle_theta = np.linspace(0, 2 * np.pi, 50)
                circle_r = np.ones_like(circle_theta) * circle_radius
                ax.plot(circle_theta, circle_r, color='black',
                        linewidth=1.5, alpha=0.9, zorder=10)
            
            ax.set_theta_zero_location("S")  # DV at bottom
            ax.set_theta_direction(1)
            ax.set_ylim(0, 1.0)
            ax.set_yticks([0.25, 0.5, 0.75])
            if col_idx == 0:
                ax.set_yticklabels(['25%', '50%', '75%'], fontsize=5)
            else:
                ax.set_yticklabels(['', '', ''], fontsize=5)
            ax.grid(True, alpha=0.3)
            
            if row_idx == 0:
                ax.set_title(f'Slice {col_idx}', fontsize=9, fontweight='bold', pad=8)
            
            if col_idx == 0:
                ax.set_xticks(angles[:-1])
                ax.set_xticklabels(categories, fontsize=7, fontweight='bold')
            else:
                ax.set_xticks([])
        
        # AVERAGE column (per-subject average across slices)
        ax_avg = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + target_n_slices + 1, projection='polar')
        ax_avg.set_facecolor('#fffde7')
        
        for subj, subj_info in group_data.items():
            data = subj_info['data']
            color = subj_info['color']
            
            subj_dv = np.mean([data[i]['fd_dv'] for i in range(target_n_slices) if i in data])
            subj_ap = np.mean([data[i]['fd_ap'] for i in range(target_n_slices) if i in data])
            subj_lm = np.mean([data[i]['fd_lm'] for i in range(target_n_slices) if i in data])
            
            values = [subj_dv / global_max_amplitude,
                      subj_ap / global_max_amplitude,
                      subj_lm / global_max_amplitude]
            values += values[:1]
            
            ax_avg.plot(angles, values, 'o-', linewidth=1.5, color=color,
                        markersize=3, alpha=0.7)
            ax_avg.fill(angles, values, alpha=0.1, color=color)
        
        avg_fd_all = np.mean(avg_fd_per_group[group])
        circle_radius_avg = 0.25 * np.sqrt(avg_fd_all / global_max_fd) if global_max_fd > 0 else 0
        if circle_radius_avg > 0.005:
            circle_theta = np.linspace(0, 2 * np.pi, 50)
            circle_r = np.ones_like(circle_theta) * circle_radius_avg
            ax_avg.plot(circle_theta, circle_r, color='black',
                        linewidth=1.5, alpha=0.9, zorder=10)
        
        ax_avg.set_theta_zero_location("S")
        ax_avg.set_theta_direction(1)
        ax_avg.set_ylim(0, 1.0)
        ax_avg.set_yticks([0.25, 0.5, 0.75])
        ax_avg.set_yticklabels(['', '', ''], fontsize=5)
        ax_avg.grid(True, alpha=0.3)
        ax_avg.set_xticks([])
        if row_idx == 0:
            ax_avg.set_title('AVG/subj', fontsize=8, fontweight='bold', pad=8, color='darkred')
        
        # GROUP AVERAGE column
        ax_grp = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + target_n_slices + 2, projection='polar')
        ax_grp.set_facecolor('#ffecb3')
        
        grp_dv = np.mean(avg_dir_per_group[group]['dv'])
        grp_ap = np.mean(avg_dir_per_group[group]['ap'])
        grp_lm = np.mean(avg_dir_per_group[group]['lm'])
        
        grp_values = [grp_dv / global_max_amplitude,
                      grp_ap / global_max_amplitude,
                      grp_lm / global_max_amplitude]
        grp_values += grp_values[:1]
        
        ax_grp.plot(angles, grp_values, 'o-', linewidth=3, color='black',
                    markersize=6, alpha=0.9)
        ax_grp.fill(angles, grp_values, alpha=0.3, color=bg_color)
        
        grp_fd = np.mean(avg_fd_per_group[group])
        circle_radius_grp = 0.25 * np.sqrt(grp_fd / global_max_fd) if global_max_fd > 0 else 0
        if circle_radius_grp > 0.005:
            circle_theta = np.linspace(0, 2 * np.pi, 50)
            circle_r = np.ones_like(circle_theta) * circle_radius_grp
            ax_grp.plot(circle_theta, circle_r, color='darkred',
                        linewidth=2, alpha=0.9, zorder=10)
        
        ax_grp.set_theta_zero_location("S")
        ax_grp.set_theta_direction(1)
        ax_grp.set_ylim(0, 1.0)
        ax_grp.set_yticks([0.25, 0.5, 0.75])
        ax_grp.set_yticklabels(['', '', ''], fontsize=5)
        ax_grp.grid(True, alpha=0.3)
        ax_grp.set_xticks([])
        if row_idx == 0:
            ax_grp.set_title('GROUP\nAVG', fontsize=8, fontweight='bold', pad=8, color='darkred')
        
        # Info panel
        ax_info = fig.add_subplot(n_rows, n_cols, row_idx * n_cols + n_cols)
        ax_info.set_facecolor(bg_color)
        ax_info.axis('off')
        
        total_fd = np.mean(avg_fd_per_group[group])
        total_ap = np.mean(avg_dir_per_group[group]['ap'])
        total_dv = np.mean(avg_dir_per_group[group]['dv'])
        total_lm = np.mean(avg_dir_per_group[group]['lm'])
        total_dir = total_ap + total_dv + total_lm
        
        info_text = f"{label}\n\n"
        info_text += f"FD total: {total_fd:.1f}\n({100*total_fd/global_max_fd:.0f}% of max)\n\n"
        info_text += f"Dir AFD: {total_dir:.1f}\n(Σ|dir|×AFD)\n\n"
        info_text += f"Direction %:\n"
        if total_dir > 0:
            info_text += f"  DV: {100*total_dv/total_dir:.0f}%\n"
            info_text += f"  AP: {100*total_ap/total_dir:.0f}%\n"
            info_text += f"  ML: {100*total_lm/total_dir:.0f}%"
        
        ax_info.text(0.5, 0.5, info_text, transform=ax_info.transAxes,
                     fontsize=9, ha='center', va='center',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Scale legend at bottom
    ax_legend = fig.add_axes([0.1, 0.02, 0.8, 0.04])
    ax_legend.axis('off')
    scale_text = (f"SCALE: Radar = directional AFD (max={global_max_amplitude:.1f}). "
                  f"Circle radius ∝ √(FD/max), max FD={global_max_fd:.1f}. "
                  f"Units: AFD = Apparent Fiber Density from fod2fixel.")
    ax_legend.text(0.5, 0.5, scale_text, transform=ax_legend.transAxes,
                   fontsize=8, ha='center', va='center',
                   bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.9))
    
    # Color legend
    legend_elements = [plt.Line2D([0], [0], color=c, linewidth=2, marker='o',
                                   label=f'Mouse {i+1}')
                       for i, c in enumerate(colors)]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=8,
               bbox_to_anchor=(0.99, 0.98), title='Subjects')
    
    shrink_title = " (Shrunk ROI -1 voxel)" if use_shrink else ""
    fig.suptitle(f'Fixel Directional AFD & Fiber Density - {hemisphere.upper()} Hemisphere{shrink_title}\n'
                 f'Radar: directional AFD (DV/AP/ML), Circle: total FD',
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.06, 1, 0.95])
    
    suffix = '_shrink' if use_shrink else ''
    output_path = os.path.join(output_dir, f'fixel_3group_comparison_{hemisphere}{suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    logger.info(f"\nSaved: {output_path}")
    plt.close()
    
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='ROI-based fixel analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--subject', help='Process specific subject')
    parser.add_argument('--group', help='Process all subjects in group (CTRL, AP, APTR)')
    parser.add_argument('--roi',
                        help='ROI name (e.g., graft_left, isocortex_left)')
    parser.add_argument('--probe', action='store_true',
                        help='Generate fixel probe comparison (3 groups, single-slice ROIs)')
    parser.add_argument('--compare', action='store_true',
                        help='Generate fixel-based 3-group comparison (CTRL vs AP vs APTR)')
    parser.add_argument('--shrink', action='store_true',
                        help='Use shrunk graft ROIs (eroded by 1 voxel) for AP/APTR')
    parser.add_argument('--output', default='reports/fixel_analysis.csv',
                        help='Output CSV file')
    parser.add_argument('--force', action='store_true',
                        help='Force regeneration of fixels')
    parser.add_argument('--log-file', help='Path to log file')
    
    args = parser.parse_args()
    
    # Setup logger
    global logger
    logger = setup_logger(__name__, args.log_file)
    
    # Load config
    cfg = load_config()
    
    logger.info("="*60)
    logger.info("ROI-Based Fixel Analysis")
    logger.info("="*60)
    
    if args.probe:
        logger.info("PROBE MODE: Fixel-based comparison using single-slice probe ROIs")
        output_dir = 'reports'
        result = analyze_probe_fixels(cfg, output_dir, force=args.force)
        if result:
            logger.info(f"\nDone! Figure saved to: {result}")
        return 0
    
    if args.compare:
        logger.info("COMPARE MODE: Fixel-based 3-group slicewise comparison")
        output_dir = 'reports'
        for hemisphere in ['left', 'right']:
            analyze_fixel_3group_comparison(hemisphere, cfg, output_dir,
                                            use_shrink=args.shrink, force=args.force)
        return 0
    
    # Standard ROI mode
    if not args.roi:
        parser.error("Must specify --roi, --probe, or --compare")
    
    logger.info(f"ROI: {args.roi}")
    
    if not args.group and not args.subject:
        parser.error("Must specify --group or --subject")
    
    # Get subjects
    subjects = get_subjects(cfg, group=args.group, subject=args.subject)
    logger.info(f"Processing {len(subjects)} subjects...")
    
    # Process each subject
    results = []
    for subject, session, group in subjects:
        metrics = process_subject(subject, session, group, args.roi, cfg, args.force)
        if metrics:
            results.append(metrics)
    
    # Save results
    output_path = Path(args.output)
    save_results_csv(results, output_path)
    
    logger.info(f"\n✓ Processed {len(results)}/{len(subjects)} subjects successfully")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
