#!/usr/bin/env python3
"""
Generate HTML report with tract statistics for all subjects.
Calculates FA and FOD values along each tract.
"""

import subprocess
import os
import json
from pathlib import Path
from datetime import datetime
import numpy as np

from config_loader import load_config, get_paths, get_subject_groups


def run_cmd(cmd, capture=True):
    """Run shell command and return output."""
    result = subprocess.run(cmd, shell=True, capture_output=capture, text=True)
    if capture:
        return result.stdout.strip()
    return result.returncode


def generate_fa_map(dwi_path, fa_path, mask_path):
    """Generate FA map from DWI if it doesn't exist."""
    if not os.path.exists(fa_path):
        print(f"    Generating FA map...")
        # First create tensor
        tensor_path = fa_path.replace('fa_flip2.mif', 'tensor_flip2.mif')
        run_cmd(f"dwi2tensor {dwi_path} {tensor_path} -mask {mask_path} -force", capture=False)
        # Then extract FA
        run_cmd(f"tensor2metric {tensor_path} -fa {fa_path} -force", capture=False)
    return fa_path


def sample_tract_values(tract_path, image_path, stat='mean'):
    """Sample values from an image along a tract and return statistics."""
    if not os.path.exists(tract_path) or not os.path.exists(image_path):
        return None
    
    # Use tcksample to get values along tract (write to temp file)
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp_path = tmp.name
    
    try:
        run_cmd(f"tcksample {tract_path} {image_path} {tmp_path} -stat_tck {stat} -force 2>/dev/null", capture=False)
        
        if not os.path.exists(tmp_path):
            return None
        
        with open(tmp_path, 'r') as f:
            output = f.read()
        
        # Parse values - skip comment lines and handle both space and newline separated
        values = []
        for line in output.strip().split('\n'):
            if line.startswith('#'):
                continue
            # Values can be space-separated on one line
            for val_str in line.strip().split():
                try:
                    val = float(val_str)
                    if not np.isnan(val):
                        values.append(val)
                except ValueError:
                    continue
        
        if len(values) == 0:
            return None
        
        return {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'n_streamlines': len(values)
        }
    except Exception as e:
        return None
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def get_tract_count(tract_path):
    """Get number of streamlines in a tract file."""
    if not os.path.exists(tract_path):
        return 0
    output = run_cmd(f"tckinfo {tract_path} 2>/dev/null | grep 'count:' | head -1")
    try:
        return int(output.split(':')[1].strip())
    except:
        return 0


def get_subject_info():
    """Get all subjects organized by group."""
    groups = get_subject_groups()
    return groups


def main():
    config = load_config()
    paths = get_paths()
    deriv_dir = Path(paths['derivatives_dir'])
    
    print("=" * 60)
    print("GENERATING TRACT STATISTICS REPORT")
    print("=" * 60)
    
    groups = get_subject_info()
    
    # Define tract types
    tracts = [
        ('thalamus', 'left'),
        ('thalamus', 'right'),
        ('striatum', 'left'),
        ('striatum', 'right')
    ]
    
    # Collect all statistics
    all_stats = {}
    
    for group_name, subjects in groups.items():
        print(f"\n### {group_name} GROUP ###")
        all_stats[group_name] = {}
        
        for subj in subjects:
            subj_id = subj['id']
            session = subj['session']
            
            print(f"\n  Subject: {subj_id}")
            
            base_path = deriv_dir / subj_id / session / "dwi"
            
            # Check if subject has tracts
            tract_dir = base_path / "segment_tracts"
            if not tract_dir.exists():
                print(f"    No tract directory found, skipping...")
                continue
            
            # Generate FA map if needed
            dwi_path = base_path / "dwi_flip2.mif"
            fa_path = base_path / "fa_flip2.mif"
            mask_path = base_path / "mask_flip2.mif"
            fod_path = base_path / "wm_fod_flip2.mif"
            
            if dwi_path.exists() and mask_path.exists():
                generate_fa_map(str(dwi_path), str(fa_path), str(mask_path))
            
            all_stats[group_name][subj_id] = {
                'session': session,
                'tracts': {}
            }
            
            # Determine if AP/APTR (graft) or CTRL (isocortex)
            is_graft = group_name in ['AP', 'APTR']
            seed_prefix = 'graft' if is_graft else 'isocortex'
            
            for target, hemi in tracts:
                tract_name = f"{seed_prefix}_{hemi}_anterior_dorsal_to_{target}_{hemi}"
                tract_file = tract_dir / f"tracks_{tract_name}_cleaned.tck"
                
                if not tract_file.exists():
                    print(f"    {hemi} {target}: NOT FOUND")
                    all_stats[group_name][subj_id]['tracts'][f"{target}_{hemi}"] = None
                    continue
                
                count = get_tract_count(str(tract_file))
                
                # Sample FA values
                fa_stats = None
                if fa_path.exists():
                    fa_stats = sample_tract_values(str(tract_file), str(fa_path))
                
                # Sample FOD amplitude (L0 component approximation using first volume)
                fod_stats = None
                if fod_path.exists():
                    # Extract first SH coefficient (L=0, approximates overall FOD amplitude)
                    fod_l0_path = base_path / "fod_l0_temp.mif"
                    run_cmd(f"mrconvert {fod_path} {fod_l0_path} -coord 3 0 -force", capture=False)
                    fod_stats = sample_tract_values(str(tract_file), str(fod_l0_path))
                    if fod_l0_path.exists():
                        os.remove(fod_l0_path)
                
                all_stats[group_name][subj_id]['tracts'][f"{target}_{hemi}"] = {
                    'count': count,
                    'fa': fa_stats,
                    'fod': fod_stats
                }
                
                fa_str = f"FA={fa_stats['mean']:.3f}±{fa_stats['std']:.3f}" if fa_stats else "FA=N/A"
                fod_str = f"FOD={fod_stats['mean']:.3f}±{fod_stats['std']:.3f}" if fod_stats else "FOD=N/A"
                print(f"    {hemi} {target}: {count} fibers, {fa_str}, {fod_str}")
    
    # Generate HTML report
    generate_html_report(all_stats, deriv_dir)
    
    print("\n" + "=" * 60)
    print("REPORT GENERATED!")
    print("=" * 60)


def generate_html_report(all_stats, deriv_dir):
    """Generate HTML report with all statistics."""
    
    report_path = deriv_dir / "tract_statistics_report.html"
    
    html = """<!DOCTYPE html>
<html>
<head>
    <title>Mouse Brain Tractography Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
        }
        h3 {
            color: #7f8c8d;
        }
        .group-section {
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        .subject-card {
            background: #ecf0f1;
            padding: 15px;
            margin: 10px 0;
            border-radius: 5px;
            border-left: 4px solid #3498db;
        }
        .ctrl-card { border-left-color: #27ae60; }
        .ap-card { border-left-color: #e74c3c; }
        .aptr-card { border-left-color: #9b59b6; }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 10px 0;
        }
        th, td {
            padding: 8px 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background: #3498db;
            color: white;
        }
        tr:hover {
            background: #f1f1f1;
        }
        .stat-value {
            font-family: monospace;
            font-size: 0.9em;
        }
        .na-value {
            color: #95a5a6;
            font-style: italic;
        }
        .summary-table {
            margin-top: 30px;
        }
        .summary-table th {
            background: #2c3e50;
        }
        .good { color: #27ae60; }
        .warning { color: #f39c12; }
        .bad { color: #e74c3c; }
        .timestamp {
            color: #7f8c8d;
            font-size: 0.9em;
            margin-bottom: 20px;
        }
        .legend {
            background: #fff;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }
        .legend-item {
            display: inline-block;
            margin-right: 20px;
        }
        .legend-color {
            display: inline-block;
            width: 20px;
            height: 20px;
            margin-right: 5px;
            vertical-align: middle;
            border-radius: 3px;
        }
    </style>
</head>
<body>
    <h1>🧠 Mouse Brain Tractography Report</h1>
    <p class="timestamp">Generated: """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
    
    <div class="legend">
        <strong>Groups:</strong>
        <span class="legend-item"><span class="legend-color" style="background:#27ae60"></span> CTRL (Allen Atlas)</span>
        <span class="legend-item"><span class="legend-color" style="background:#e74c3c"></span> AP (Graft)</span>
        <span class="legend-item"><span class="legend-color" style="background:#9b59b6"></span> APTR (Graft + Treatment)</span>
    </div>
    
    <div class="legend">
        <strong>Tracts:</strong> anterior_dorsal segment → thalamus (left/right) and striatum (left/right)
    </div>
"""
    
    # Group summary statistics
    html += """
    <h2>📊 Group Summary Statistics</h2>
    <div class="group-section">
"""
    
    for group_name in ['CTRL', 'AP', 'APTR']:
        if group_name not in all_stats:
            continue
            
        group_data = all_stats[group_name]
        
        # Calculate group averages
        tract_types = ['thalamus_left', 'thalamus_right', 'striatum_left', 'striatum_right']
        
        html += f"""
        <h3>{group_name} Group (n={len(group_data)})</h3>
        <table>
            <tr>
                <th>Tract</th>
                <th>Avg Count</th>
                <th>Avg FA (mean±std)</th>
                <th>Avg FOD (mean±std)</th>
            </tr>
"""
        
        for tract_type in tract_types:
            counts = []
            fa_means = []
            fod_means = []
            
            for subj_id, subj_data in group_data.items():
                if 'tracts' in subj_data and tract_type in subj_data['tracts']:
                    tract_data = subj_data['tracts'][tract_type]
                    if tract_data:
                        counts.append(tract_data['count'])
                        if tract_data.get('fa'):
                            fa_means.append(tract_data['fa']['mean'])
                        if tract_data.get('fod'):
                            fod_means.append(tract_data['fod']['mean'])
            
            avg_count = f"{np.mean(counts):.0f}±{np.std(counts):.0f}" if counts else "N/A"
            avg_fa = f"{np.mean(fa_means):.3f}±{np.std(fa_means):.3f}" if fa_means else "N/A"
            avg_fod = f"{np.mean(fod_means):.3f}±{np.std(fod_means):.3f}" if fod_means else "N/A"
            
            html += f"""
            <tr>
                <td>{tract_type.replace('_', ' ').title()}</td>
                <td class="stat-value">{avg_count}</td>
                <td class="stat-value">{avg_fa}</td>
                <td class="stat-value">{avg_fod}</td>
            </tr>
"""
        
        html += """
        </table>
"""
    
    html += """
    </div>
"""
    
    # Individual subject details
    html += """
    <h2>📋 Individual Subject Details</h2>
"""
    
    for group_name in ['CTRL', 'AP', 'APTR']:
        if group_name not in all_stats:
            continue
            
        group_data = all_stats[group_name]
        card_class = {'CTRL': 'ctrl-card', 'AP': 'ap-card', 'APTR': 'aptr-card'}[group_name]
        
        html += f"""
    <div class="group-section">
        <h2>{group_name} Subjects</h2>
"""
        
        for subj_id, subj_data in group_data.items():
            session = subj_data.get('session', 'unknown')
            
            html += f"""
        <div class="subject-card {card_class}">
            <h3>{subj_id}</h3>
            <p>Session: {session}</p>
            <table>
                <tr>
                    <th>Tract</th>
                    <th>Streamlines</th>
                    <th>FA Mean</th>
                    <th>FA Std</th>
                    <th>FOD Mean</th>
                    <th>FOD Std</th>
                </tr>
"""
            
            for tract_type in ['thalamus_left', 'thalamus_right', 'striatum_left', 'striatum_right']:
                tract_data = subj_data.get('tracts', {}).get(tract_type)
                
                if tract_data:
                    count = tract_data['count']
                    fa = tract_data.get('fa', {})
                    fod = tract_data.get('fod', {})
                    
                    fa_mean = f"{fa['mean']:.4f}" if fa else '<span class="na-value">N/A</span>'
                    fa_std = f"{fa['std']:.4f}" if fa else '<span class="na-value">N/A</span>'
                    fod_mean = f"{fod['mean']:.4f}" if fod else '<span class="na-value">N/A</span>'
                    fod_std = f"{fod['std']:.4f}" if fod else '<span class="na-value">N/A</span>'
                    
                    count_class = 'good' if count > 500 else ('warning' if count > 100 else 'bad')
                    
                    html += f"""
                <tr>
                    <td>{tract_type.replace('_', ' ').title()}</td>
                    <td class="stat-value {count_class}">{count}</td>
                    <td class="stat-value">{fa_mean}</td>
                    <td class="stat-value">{fa_std}</td>
                    <td class="stat-value">{fod_mean}</td>
                    <td class="stat-value">{fod_std}</td>
                </tr>
"""
                else:
                    html += f"""
                <tr>
                    <td>{tract_type.replace('_', ' ').title()}</td>
                    <td class="stat-value na-value">0</td>
                    <td class="stat-value na-value">N/A</td>
                    <td class="stat-value na-value">N/A</td>
                    <td class="stat-value na-value">N/A</td>
                    <td class="stat-value na-value">N/A</td>
                </tr>
"""
            
            html += """
            </table>
        </div>
"""
        
        html += """
    </div>
"""
    
    html += """
    <hr>
    <p style="color: #7f8c8d; font-size: 0.8em;">
        Generated by MBRAIN Pipeline | FA = Fractional Anisotropy | FOD = Fiber Orientation Distribution (L0 amplitude)
    </p>
</body>
</html>
"""
    
    with open(report_path, 'w') as f:
        f.write(html)
    
    print(f"\n✓ HTML report saved to: {report_path}")
    return report_path


if __name__ == "__main__":
    main()
