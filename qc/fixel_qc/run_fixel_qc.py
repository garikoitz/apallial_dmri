#!/usr/bin/env python3
"""
Fixel QC: Generate voxel-level QC maps and statistics for all subjects.

Based on Raffelt et al. 2017 (FBA framework) and Dhollander et al. 2021.

QC checks:
1. Fixel count map (voxels should have 1-3 fixels in WM)
2. AFD map (total apparent fiber density per voxel)
3. Peak amplitude map
4. Fixel complexity (fraction of multi-fixel voxels)
5. Mask coverage check (fixels should fill the mask)
6. Cross-subject consistency (outlier detection)

Outputs:
    qc/fixel_qc/fixel_qc_report.txt  - summary stats
    qc/fixel_qc/<subject>/           - per-subject voxel maps
"""

import os
import sys
import subprocess
import numpy as np
from pathlib import Path

# Add parent for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
os.environ['MBRAIN_PIPELINE'] = 'flip1'
from functions.config_loader import load_config

cfg = load_config()
DERIV = Path(cfg['paths']['derivatives_flip1'])
QC_DIR = Path(__file__).parent

SUBJECTS = {
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


def run(cmd):
    """Run shell command, return True on success."""
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  WARNING: {r.stderr[:200]}")
        return False
    return True


def generate_voxel_maps(subj, ses, fixel_dir, out_dir):
    """Generate fixel2voxel QC maps."""
    os.makedirs(out_dir, exist_ok=True)
    
    afd_file = fixel_dir / 'afd.mif'
    peak_file = fixel_dir / 'peak.mif'
    
    maps = {}
    
    # 1. Fixel count per voxel
    count_map = out_dir / 'fixel_count.mif'
    if run(f"fixel2voxel '{afd_file}' count '{count_map}' -force -quiet"):
        maps['count'] = count_map
    
    # 2. Total AFD per voxel (sum across fixels)
    if afd_file.exists():
        afd_sum = out_dir / 'afd_total.mif'
        if run(f"fixel2voxel '{afd_file}' sum '{afd_sum}' -force -quiet"):
            maps['afd_sum'] = afd_sum
        
        # 3. Max AFD per voxel
        afd_max = out_dir / 'afd_max.mif'
        if run(f"fixel2voxel '{afd_file}' max '{afd_max}' -force -quiet"):
            maps['afd_max'] = afd_max
    
    # 4. Peak amplitude per voxel (sum)
    if peak_file.exists():
        peak_sum = out_dir / 'peak_total.mif'
        if run(f"fixel2voxel '{peak_file}' sum '{peak_sum}' -force -quiet"):
            maps['peak_sum'] = peak_sum
    
    # 5. Complexity (fraction of second fixel relative to first)
    complexity = out_dir / 'complexity.mif'
    if run(f"fixel2voxel '{afd_file}' complexity '{complexity}' -force -quiet"):
        maps['complexity'] = complexity
    
    return maps


def get_mrstats(filepath):
    """Get mean, median, std, min, max from mrstats."""
    r = subprocess.run(
        f"mrstats '{filepath}' -output mean -output median -output std -output min -output max -ignorezero",
        shell=True, capture_output=True, text=True
    )
    if r.returncode != 0:
        return None
    vals = r.stdout.strip().split()
    if len(vals) >= 5:
        return {
            'mean': float(vals[0]),
            'median': float(vals[1]),
            'std': float(vals[2]),
            'min': float(vals[3]),
            'max': float(vals[4]),
        }
    return None


def count_voxels_in_mask(mask_path):
    """Count non-zero voxels in mask."""
    r = subprocess.run(
        f"mrstats '{mask_path}' -output count -mask '{mask_path}'",
        shell=True, capture_output=True, text=True
    )
    if r.returncode == 0:
        return int(float(r.stdout.strip()))
    return 0


def main():
    report_path = QC_DIR / 'fixel_qc_report.txt'
    
    all_stats = []
    
    with open(report_path, 'w') as report:
        report.write("=" * 100 + "\n")
        report.write("FIXEL QC REPORT\n")
        report.write("Based on Raffelt et al. 2017, Dhollander et al. 2021\n")
        report.write("=" * 100 + "\n\n")
        
        report.write("QC CHECKS (what to look for):\n")
        report.write("-" * 60 + "\n")
        report.write("1. FIXEL COUNT: WM voxels should have 1-3 fixels (mean ~1.2-2.0)\n")
        report.write("   - If max > 5: possible noisy FODs or bad masking\n")
        report.write("   - If mean < 1.0: poor FOD estimation\n")
        report.write("2. AFD (Apparent Fiber Density): should be positive in WM\n")
        report.write("   - Outlier subjects (>2 SD from group mean) suggest issues\n")
        report.write("   - Very low AFD: signal dropout or bad normalization\n")
        report.write("3. PEAK AMPLITUDE: FOD peak values at fixel directions\n")
        report.write("   - Should correlate with AFD\n")
        report.write("4. COMPLEXITY: fraction of voxels with crossing fibers\n")
        report.write("   - Expected range 0.1-0.4 in mouse brain\n")
        report.write("   - Very high (>0.5): noisy FODs\n")
        report.write("5. MASK COVERAGE: fixels should fill the brain mask\n")
        report.write("   - Voxels with fixels / mask voxels should be >80%\n")
        report.write("6. VISUAL CHECK in mrview:\n")
        report.write("   - fixel plot: sticks should align with known WM tracts\n")
        report.write("   - No fixels in CSF/outside brain\n")
        report.write("   - Crossing fiber regions should show 2-3 fixels\n")
        report.write("\n")
        
        # Header for stats table
        header = (f"{'Group':<6} {'Subject':<28} {'MaskVox':>8} {'FixVox':>8} "
                  f"{'Cover%':>7} {'CountM':>7} {'CountMax':>8} "
                  f"{'AFD_mean':>9} {'AFD_std':>8} "
                  f"{'Peak_m':>8} {'Cmplx_m':>8}")
        report.write(header + "\n")
        report.write("-" * len(header) + "\n")
        
        for group in ['CTRL', 'AP', 'APTR']:
            for subj, ses in SUBJECTS[group]:
                subj_short = subj.replace('sub-', '')
                subj_dir = DERIV / subj / ses
                fixel_dir = subj_dir / 'dwi' / 'fixels_group'
                mask_path = subj_dir / 'dwi' / 'mask_anat.mif'
                
                if not fixel_dir.exists():
                    report.write(f"{group:<6} {subj_short:<28} FIXELS MISSING\n")
                    continue
                
                # Generate QC maps
                out_dir = QC_DIR / subj
                maps = generate_voxel_maps(subj, ses, fixel_dir, out_dir)
                
                # Collect stats
                mask_vox = count_voxels_in_mask(mask_path)
                
                count_stats = get_mrstats(maps.get('count', '')) if 'count' in maps else None
                afd_stats = get_mrstats(maps.get('afd_sum', '')) if 'afd_sum' in maps else None
                peak_stats = get_mrstats(maps.get('peak_sum', '')) if 'peak_sum' in maps else None
                cmplx_stats = get_mrstats(maps.get('complexity', '')) if 'complexity' in maps else None
                
                # Count voxels with at least 1 fixel
                if 'count' in maps:
                    r = subprocess.run(
                        f"mrstats '{maps['count']}' -output count -ignorezero",
                        shell=True, capture_output=True, text=True
                    )
                    fix_vox = int(float(r.stdout.strip())) if r.returncode == 0 else 0
                else:
                    fix_vox = 0
                
                coverage = 100 * fix_vox / mask_vox if mask_vox > 0 else 0
                
                row = (f"{group:<6} {subj_short:<28} {mask_vox:>8} {fix_vox:>8} "
                       f"{coverage:>6.1f}% "
                       f"{count_stats['mean']:>7.2f} {count_stats['max']:>8.0f} " if count_stats else "     -        - ")
                
                if afd_stats:
                    row += f"{afd_stats['mean']:>9.3f} {afd_stats['std']:>8.3f} "
                else:
                    row += "        -        - "
                
                if peak_stats:
                    row += f"{peak_stats['mean']:>8.3f} "
                else:
                    row += "       - "
                
                if cmplx_stats:
                    row += f"{cmplx_stats['mean']:>8.3f}"
                else:
                    row += "       -"
                
                report.write(row + "\n")
                
                stats_entry = {
                    'group': group, 'subject': subj_short,
                    'mask_vox': mask_vox, 'fix_vox': fix_vox,
                    'coverage': coverage,
                    'count_mean': count_stats['mean'] if count_stats else 0,
                    'count_max': count_stats['max'] if count_stats else 0,
                    'afd_mean': afd_stats['mean'] if afd_stats else 0,
                    'afd_std': afd_stats['std'] if afd_stats else 0,
                }
                all_stats.append(stats_entry)
                
                print(f"  {group} {subj_short}: coverage={coverage:.1f}%, "
                      f"count_mean={count_stats['mean']:.2f}, "
                      f"afd_mean={afd_stats['mean']:.3f}" if count_stats and afd_stats else f"  {group} {subj_short}: incomplete")
        
        # Outlier detection
        report.write("\n" + "=" * 80 + "\n")
        report.write("OUTLIER DETECTION (>2 SD from group mean)\n")
        report.write("=" * 80 + "\n")
        
        for group in ['CTRL', 'AP', 'APTR']:
            group_stats = [s for s in all_stats if s['group'] == group]
            if len(group_stats) < 2:
                continue
            
            afd_vals = [s['afd_mean'] for s in group_stats]
            cov_vals = [s['coverage'] for s in group_stats]
            
            afd_mean = np.mean(afd_vals)
            afd_std = np.std(afd_vals)
            cov_mean = np.mean(cov_vals)
            cov_std = np.std(cov_vals)
            
            report.write(f"\n{group}: AFD mean={afd_mean:.3f} ± {afd_std:.3f}, "
                        f"Coverage={cov_mean:.1f}% ± {cov_std:.1f}%\n")
            
            for s in group_stats:
                flags = []
                if afd_std > 0 and abs(s['afd_mean'] - afd_mean) > 2 * afd_std:
                    flags.append(f"AFD={s['afd_mean']:.3f} (>{2*afd_std:.3f} from mean)")
                if cov_std > 0 and abs(s['coverage'] - cov_mean) > 2 * cov_std:
                    flags.append(f"Coverage={s['coverage']:.1f}% (>{2*cov_std:.1f} from mean)")
                if s['count_max'] > 5:
                    flags.append(f"Max fixels/voxel={s['count_max']:.0f} (>5)")
                
                if flags:
                    report.write(f"  ⚠ {s['subject']}: {'; '.join(flags)}\n")
            
            if not any(True for s in group_stats 
                      if (afd_std > 0 and abs(s['afd_mean'] - afd_mean) > 2 * afd_std) or
                         (cov_std > 0 and abs(s['coverage'] - cov_mean) > 2 * cov_std)):
                report.write(f"  ✓ No outliers detected\n")
        
        # mrview commands for visual inspection
        report.write("\n" + "=" * 80 + "\n")
        report.write("VISUAL QC COMMANDS (copy-paste to terminal)\n")
        report.write("=" * 80 + "\n\n")
        
        report.write("# --- FIXEL VISUALIZATION (sticks colored by AFD) ---\n")
        report.write("# Look for: sticks aligned with WM tracts, no fixels in CSF\n\n")
        
        # One representative per group
        reps = {
            'CTRL': ('sub-129-921-R1F', 'ses-260115121231'),
            'AP': ('sub-769-C-O5-R1M', 'ses-240911115733'),
            'APTR': ('sub-AP72-862-901-NmF', 'ses-251210093342'),
        }
        
        for group, (subj, ses) in reps.items():
            base = DERIV / subj / ses
            t2w = base / 't2w_highres_flip.mif'
            fixdir = base / 'dwi' / 'fixels_group'
            report.write(f"# {group}: {subj}\n")
            report.write(f"mrview {t2w} -fixel.load {fixdir}/afd.mif &\n\n")
        
        report.write("# --- VOXEL MAP QC (overlay on T2w) ---\n")
        report.write("# Look for: uniform AFD in WM, count=1-3, no holes in mask\n\n")
        
        for group, (subj, ses) in reps.items():
            base = DERIV / subj / ses
            t2w = base / 't2w_highres_flip.mif'
            qc_subj = QC_DIR / subj
            report.write(f"# {group}: {subj} - AFD map\n")
            report.write(f"mrview {t2w} -overlay.load {qc_subj}/afd_total.mif -overlay.opacity 0.5 &\n\n")
            report.write(f"# {group}: {subj} - Fixel count map\n")
            report.write(f"mrview {t2w} -overlay.load {qc_subj}/fixel_count.mif -overlay.opacity 0.5 &\n\n")
        
        report.write("# --- FOD CHECK (view FODs directly) ---\n")
        report.write("# Look for: smooth FODs, no spurious lobes, FODs match anatomy\n\n")
        
        for group, (subj, ses) in reps.items():
            base = DERIV / subj / ses
            t2w = base / 't2w_highres_flip.mif'
            fod = base / 'dwi' / 'wm_fod_norm_group.mif'
            report.write(f"# {group}: {subj} - FOD glyphs\n")
            report.write(f"mrview {t2w} -odf.load_sh {fod} &\n\n")
    
    print(f"\nReport saved: {report_path}")
    print(f"Per-subject QC maps saved in: {QC_DIR}/sub-*/")


if __name__ == '__main__':
    main()
