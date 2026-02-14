#!/usr/bin/env python3
"""
Generate v12 HTML Report for SHRINK Passthrough Tractography.

This report presents passthrough connectivity results from SHRUNK graft segments
(eroded by 1 voxel) for AP and APTR groups only.

Usage:
    python functions/generate_shrink_report_v12.py

Output:
    reports/wholebrain_tractography_report_v12.html
"""

import os
import sys
import subprocess
import json
from pathlib import Path
from datetime import datetime
import numpy as np

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config

DERIV = "/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/Shared drives/DTI Collab/mbrain_data/data/derivatives/mrtrix_preproc_flip1"

# Only AP and APTR subjects for shrink analysis
subjects = {
    'AP': ['sub-596-C-OO-NmM', 'sub-769-C-O5-R1M', 'sub-769-C-OO-NmM', 'sub-875-C-OO-NM'],
    'APTR': ['sub-AP72-862-901-NmF', 'sub-AP72-862-901-R1M', 'sub-Apallial17NmF', 'sub-Apallial20R1F']
}

# Graft segments (3 dorsal segments)
graft_segments = ['anterior_dorsal', 'central_dorsal', 'posterior_dorsal']

# Targets for AP/APTR
targets = ['Thalamus', 'Striatum', 'cerebellum', 'Palecortex']


def get_session(subj_path):
    """Get the session folder for a subject."""
    try:
        sessions = [d for d in os.listdir(subj_path) if d.startswith('ses-')]
        return sessions[0] if sessions else None
    except:
        return None


def get_tck_count(tck_path):
    """Get streamline count from a .tck file."""
    if not os.path.exists(tck_path):
        return 0
    try:
        result = subprocess.run(['tckinfo', tck_path], capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if 'count:' in line and 'total_count' not in line:
                return int(line.split(':')[1].strip())
    except:
        pass
    return 0


def collect_shrink_data(hemisphere='left'):
    """Collect tract data from passthrough_tracts_shrink folder."""
    hemi = hemisphere
    data = {'AP': {}, 'APTR': {}}
    
    for group, subj_list in subjects.items():
        for subj in subj_list:
            subj_path = os.path.join(DERIV, subj)
            ses = get_session(subj_path)
            if not ses:
                print(f"  WARNING: No session found for {subj}")
                continue
            
            # Use passthrough_tracts_shrink folder
            tract_dir = os.path.join(subj_path, ses, 'dwi', 'wholebrain', 'passthrough_tracts_shrink')
            
            subj_data = {}
            
            for seg in graft_segments:
                seg_data = {}
                for tgt in targets:
                    # Construct tract filename
                    tck_name = f"seg9_graft_{hemi}_{seg}_to_{tgt}_{hemi}_cleaned.tck"
                    tck_path = os.path.join(tract_dir, tck_name)
                    
                    # Also check raw (uncleaned) if cleaned doesn't exist
                    if not os.path.exists(tck_path):
                        tck_name_raw = f"seg9_graft_{hemi}_{seg}_to_{tgt}_{hemi}.tck"
                        tck_path_raw = os.path.join(tract_dir, tck_name_raw)
                        if os.path.exists(tck_path_raw):
                            count = get_tck_count(tck_path_raw)
                            seg_data[tgt.lower()] = {'cleaned': 0, 'raw': count}
                            continue
                    
                    count = get_tck_count(tck_path)
                    
                    # Also get raw count for comparison
                    tck_name_raw = f"seg9_graft_{hemi}_{seg}_to_{tgt}_{hemi}.tck"
                    tck_path_raw = os.path.join(tract_dir, tck_name_raw)
                    raw_count = get_tck_count(tck_path_raw)
                    
                    seg_data[tgt.lower()] = {'cleaned': count, 'raw': raw_count}
                
                subj_data[seg] = seg_data
            
            # Calculate totals
            total_cleaned = sum(
                sum(sd.get(t, {}).get('cleaned', 0) for t in ['thalamus', 'striatum', 'cerebellum', 'palecortex'])
                for sd in subj_data.values()
            )
            total_raw = sum(
                sum(sd.get(t, {}).get('raw', 0) for t in ['thalamus', 'striatum', 'cerebellum', 'palecortex'])
                for sd in subj_data.values()
            )
            subj_data['total_cleaned'] = total_cleaned
            subj_data['total_raw'] = total_raw
            
            data[group][subj.replace('sub-', '')] = subj_data
    
    return data


def generate_report_html(left_data, right_data):
    """Generate the v12 HTML report."""
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Calculate group statistics
    def calc_group_stats(data, group):
        subjects = data[group]
        totals_cleaned = [s['total_cleaned'] for s in subjects.values()]
        totals_raw = [s['total_raw'] for s in subjects.values()]
        return {
            'n': len(subjects),
            'mean_cleaned': np.mean(totals_cleaned) if totals_cleaned else 0,
            'std_cleaned': np.std(totals_cleaned) if totals_cleaned else 0,
            'mean_raw': np.mean(totals_raw) if totals_raw else 0,
            'std_raw': np.std(totals_raw) if totals_raw else 0,
            'subjects': subjects
        }
    
    ap_left = calc_group_stats(left_data, 'AP')
    aptr_left = calc_group_stats(left_data, 'APTR')
    ap_right = calc_group_stats(right_data, 'AP')
    aptr_right = calc_group_stats(right_data, 'APTR')
    
    # Generate subject rows
    def gen_subject_rows(subjects, group_name):
        rows = []
        for subj, data in subjects.items():
            rows.append(f"""
                <tr>
                    <td>{subj}</td>
                    <td>{group_name}</td>
                    <td>{data['total_raw']:,}</td>
                    <td>{data['total_cleaned']:,}</td>
                    <td>{data['total_cleaned']/data['total_raw']*100:.1f}%</td>
                </tr>
            """)
        return '\n'.join(rows)
    
    # Generate detailed segment/target breakdown
    def gen_detail_rows(data, group):
        rows = []
        for subj, subj_data in data[group].items():
            for seg in graft_segments:
                if seg not in subj_data:
                    continue
                seg_data = subj_data[seg]
                for tgt in ['thalamus', 'striatum', 'cerebellum', 'palecortex']:
                    if tgt not in seg_data:
                        continue
                    cleaned = seg_data[tgt].get('cleaned', 0)
                    raw = seg_data[tgt].get('raw', 0)
                    rows.append(f"""
                        <tr>
                            <td>{subj}</td>
                            <td>{seg}</td>
                            <td>{tgt}</td>
                            <td>{raw:,}</td>
                            <td>{cleaned:,}</td>
                        </tr>
                    """)
        return '\n'.join(rows)
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MBRAIN Whole-Brain Tractography Report v12 - SHRUNK SEGMENTS</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --ap-color: #27ae60;
            --aptr-color: #e74c3c;
            --left-color: #9b59b6;
            --right-color: #e67e22;
            --bg-color: #f5f6fa;
            --card-bg: #ffffff;
            --text-color: #2c3e50;
            --border-color: #dcdde1;
        }}
        
        * {{ box-sizing: border-box; }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
            margin: 0;
            padding: 20px;
        }}
        
        .container {{ max-width: 1400px; margin: 0 auto; }}
        
        header {{
            background: linear-gradient(135deg, #8e44ad, #9b59b6);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        
        header h1 {{ margin: 0 0 10px 0; font-size: 2.5em; }}
        header p {{ margin: 0; opacity: 0.9; font-size: 1.1em; }}
        
        .card {{
            background: var(--card-bg);
            border-radius: 10px;
            padding: 25px;
            margin-bottom: 25px;
            box-shadow: 0 2px 15px rgba(0,0,0,0.1);
        }}
        
        .card h2 {{
            color: var(--primary-color);
            border-bottom: 2px solid var(--secondary-color);
            padding-bottom: 10px;
            margin-top: 0;
        }}
        
        .methodology-note {{
            background: #f3e5f5;
            border-left: 4px solid #8e44ad;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        
        .key-insight {{
            background: #e8f6ff;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 15px 0;
            border-radius: 0 8px 8px 0;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        
        th {{
            background: var(--primary-color);
            color: white;
            font-weight: 600;
        }}
        
        tr:hover {{ background: #f5f6fa; }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .stat-box {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
        }}
        
        .stat-box.ap {{ background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); }}
        .stat-box.aptr {{ background: linear-gradient(135deg, #eb3349 0%, #f45c43 100%); }}
        
        .stat-value {{ font-size: 2.5em; font-weight: bold; }}
        .stat-label {{ opacity: 0.9; margin-top: 5px; }}
        
        .hemisphere-section {{
            border: 3px solid;
            border-radius: 12px;
            padding: 25px;
            margin-bottom: 30px;
        }}
        
        .hemisphere-section.left {{ border-color: var(--left-color); }}
        .hemisphere-section.right {{ border-color: var(--right-color); }}
        
        .hemisphere-header {{
            padding: 15px 25px;
            margin: -25px -25px 20px -25px;
            border-radius: 8px 8px 0 0;
            color: white;
        }}
        
        .hemisphere-header.left {{ background: var(--left-color); }}
        .hemisphere-header.right {{ background: var(--right-color); }}
        
        .hemisphere-header h2 {{ color: white; border: none; margin: 0; padding: 0; }}
        
        .chart-container {{ height: 400px; margin: 20px 0; }}
        
        .group-label {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        
        .group-label.ap {{ background: #d5f5e3; color: #1e8449; }}
        .group-label.aptr {{ background: #fadbd8; color: #c0392b; }}
        
        footer {{
            text-align: center;
            padding: 20px;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧠 MBRAIN Tractography Report v12</h1>
            <p>SHRUNK SEGMENTS Analysis - AP and APTR Groups</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: {timestamp}</p>
        </header>
        
        <div class="card">
            <h2>📋 Report Overview</h2>
            <div class="methodology-note">
                <strong>🔬 SHRUNK SEGMENTS Method:</strong>
                <p>This report presents passthrough tractography results using <strong>eroded graft segments</strong> (shrunk by 1 voxel). 
                This erosion removes edge voxels to reduce spurious connections from partial volume effects at segment boundaries.</p>
                <ul>
                    <li><strong>Erosion:</strong> Binary erosion with 3D ball structuring element (6-connectivity)</li>
                    <li><strong>Segments used:</strong> 3 dorsal segments (anterior_dorsal, central_dorsal, posterior_dorsal)</li>
                    <li><strong>Targets:</strong> Thalamus, Striatum, Cerebellum, Paleocortex</li>
                    <li><strong>Groups:</strong> AP (n=4) and APTR (n=4) - CTRL not applicable (no graft)</li>
                </ul>
            </div>
            
            <div class="key-insight">
                <strong>💡 Key Changes from v11:</strong>
                <ul>
                    <li>Graft segments eroded by 1 voxel to reduce edge effects</li>
                    <li>Output folder: <code>passthrough_tracts_shrink/</code></li>
                    <li>Expected reduction in streamline counts due to smaller ROIs</li>
                </ul>
            </div>
        </div>
        
        <!-- Group Summary Statistics -->
        <div class="card">
            <h2>📊 Group Summary Statistics</h2>
            
            <div class="stats-grid">
                <div class="stat-box ap">
                    <div class="stat-value">{ap_left['mean_cleaned'] + ap_right['mean_cleaned']:.0f}</div>
                    <div class="stat-label">AP Mean Total (Both Hemispheres)</div>
                </div>
                <div class="stat-box aptr">
                    <div class="stat-value">{aptr_left['mean_cleaned'] + aptr_right['mean_cleaned']:.0f}</div>
                    <div class="stat-label">APTR Mean Total (Both Hemispheres)</div>
                </div>
            </div>
            
            <table>
                <tr>
                    <th>Group</th>
                    <th>N</th>
                    <th>Left Mean ± SD</th>
                    <th>Right Mean ± SD</th>
                    <th>Total Mean</th>
                </tr>
                <tr>
                    <td><span class="group-label ap">AP</span></td>
                    <td>{ap_left['n']}</td>
                    <td>{ap_left['mean_cleaned']:.0f} ± {ap_left['std_cleaned']:.0f}</td>
                    <td>{ap_right['mean_cleaned']:.0f} ± {ap_right['std_cleaned']:.0f}</td>
                    <td>{ap_left['mean_cleaned'] + ap_right['mean_cleaned']:.0f}</td>
                </tr>
                <tr>
                    <td><span class="group-label aptr">APTR</span></td>
                    <td>{aptr_left['n']}</td>
                    <td>{aptr_left['mean_cleaned']:.0f} ± {aptr_left['std_cleaned']:.0f}</td>
                    <td>{aptr_right['mean_cleaned']:.0f} ± {aptr_right['std_cleaned']:.0f}</td>
                    <td>{aptr_left['mean_cleaned'] + aptr_right['mean_cleaned']:.0f}</td>
                </tr>
            </table>
        </div>
        
        <!-- LEFT HEMISPHERE -->
        <div class="card hemisphere-section left">
            <div class="hemisphere-header left">
                <h2>🧠 LEFT Hemisphere Results</h2>
            </div>
            
            <h3>Subject Totals (Cleaned Tracts)</h3>
            <table>
                <tr>
                    <th>Subject</th>
                    <th>Group</th>
                    <th>Raw Streamlines</th>
                    <th>Cleaned Streamlines</th>
                    <th>Retention %</th>
                </tr>
                {gen_subject_rows(ap_left['subjects'], 'AP')}
                {gen_subject_rows(aptr_left['subjects'], 'APTR')}
            </table>
            
            <div id="chart-left" class="chart-container"></div>
        </div>
        
        <!-- RIGHT HEMISPHERE -->
        <div class="card hemisphere-section right">
            <div class="hemisphere-header right">
                <h2>🧠 RIGHT Hemisphere Results</h2>
            </div>
            
            <h3>Subject Totals (Cleaned Tracts)</h3>
            <table>
                <tr>
                    <th>Subject</th>
                    <th>Group</th>
                    <th>Raw Streamlines</th>
                    <th>Cleaned Streamlines</th>
                    <th>Retention %</th>
                </tr>
                {gen_subject_rows(ap_right['subjects'], 'AP')}
                {gen_subject_rows(aptr_right['subjects'], 'APTR')}
            </table>
            
            <div id="chart-right" class="chart-container"></div>
        </div>
        
        <!-- Detailed Breakdown -->
        <div class="card">
            <h2>📋 Detailed Breakdown - Left Hemisphere</h2>
            <table>
                <tr>
                    <th>Subject</th>
                    <th>Segment</th>
                    <th>Target</th>
                    <th>Raw</th>
                    <th>Cleaned</th>
                </tr>
                {gen_detail_rows(left_data, 'AP')}
                {gen_detail_rows(left_data, 'APTR')}
            </table>
        </div>
        
        <div class="card">
            <h2>📋 Detailed Breakdown - Right Hemisphere</h2>
            <table>
                <tr>
                    <th>Subject</th>
                    <th>Segment</th>
                    <th>Target</th>
                    <th>Raw</th>
                    <th>Cleaned</th>
                </tr>
                {gen_detail_rows(right_data, 'AP')}
                {gen_detail_rows(right_data, 'APTR')}
            </table>
        </div>
        
        <footer>
            <p>MBRAIN Project | Passthrough Tractography Analysis with SHRUNK Graft Segments</p>
            <p>Pipeline: flip1 | MATLAB AFQ Cleaning | Report Version: v12</p>
        </footer>
    </div>
    
    <!-- Charts -->
    <script>
        // Left hemisphere chart
        var ap_left_subjects = {json.dumps([s for s in ap_left['subjects'].keys()])};
        var ap_left_cleaned = {json.dumps([s['total_cleaned'] for s in ap_left['subjects'].values()])};
        var aptr_left_subjects = {json.dumps([s for s in aptr_left['subjects'].keys()])};
        var aptr_left_cleaned = {json.dumps([s['total_cleaned'] for s in aptr_left['subjects'].values()])};
        
        var trace_ap_left = {{
            x: ap_left_subjects,
            y: ap_left_cleaned,
            name: 'AP',
            type: 'bar',
            marker: {{ color: '#27ae60' }}
        }};
        
        var trace_aptr_left = {{
            x: aptr_left_subjects,
            y: aptr_left_cleaned,
            name: 'APTR',
            type: 'bar',
            marker: {{ color: '#e74c3c' }}
        }};
        
        Plotly.newPlot('chart-left', [trace_ap_left, trace_aptr_left], {{
            title: 'Left Hemisphere - Cleaned Streamlines per Subject',
            barmode: 'group',
            xaxis: {{ title: 'Subject' }},
            yaxis: {{ title: 'Cleaned Streamlines' }}
        }});
        
        // Right hemisphere chart
        var ap_right_subjects = {json.dumps([s for s in ap_right['subjects'].keys()])};
        var ap_right_cleaned = {json.dumps([s['total_cleaned'] for s in ap_right['subjects'].values()])};
        var aptr_right_subjects = {json.dumps([s for s in aptr_right['subjects'].keys()])};
        var aptr_right_cleaned = {json.dumps([s['total_cleaned'] for s in aptr_right['subjects'].values()])};
        
        var trace_ap_right = {{
            x: ap_right_subjects,
            y: ap_right_cleaned,
            name: 'AP',
            type: 'bar',
            marker: {{ color: '#27ae60' }}
        }};
        
        var trace_aptr_right = {{
            x: aptr_right_subjects,
            y: aptr_right_cleaned,
            name: 'APTR',
            type: 'bar',
            marker: {{ color: '#e74c3c' }}
        }};
        
        Plotly.newPlot('chart-right', [trace_ap_right, trace_aptr_right], {{
            title: 'Right Hemisphere - Cleaned Streamlines per Subject',
            barmode: 'group',
            xaxis: {{ title: 'Subject' }},
            yaxis: {{ title: 'Cleaned Streamlines' }}
        }});
    </script>
</body>
</html>
"""
    return html


def main():
    print("=" * 60)
    print("GENERATING v12 SHRINK PASSTHROUGH REPORT")
    print("=" * 60)
    
    # Collect data for both hemispheres
    print("\nCollecting LEFT hemisphere data...")
    left_data = collect_shrink_data('left')
    
    print("\nCollecting RIGHT hemisphere data...")
    right_data = collect_shrink_data('right')
    
    # Print summary
    print("\n=== LEFT HEMISPHERE SUMMARY ===")
    for group in ['AP', 'APTR']:
        totals = [v['total_cleaned'] for v in left_data[group].values()]
        print(f"{group}: {totals} (mean={np.mean(totals):.0f})")
    
    print("\n=== RIGHT HEMISPHERE SUMMARY ===")
    for group in ['AP', 'APTR']:
        totals = [v['total_cleaned'] for v in right_data[group].values()]
        print(f"{group}: {totals} (mean={np.mean(totals):.0f})")
    
    # Generate HTML report
    print("\nGenerating HTML report...")
    html = generate_report_html(left_data, right_data)
    
    # Save report
    reports_dir = Path(__file__).parent.parent / 'reports'
    reports_dir.mkdir(parents=True, exist_ok=True)
    output_path = reports_dir / 'wholebrain_tractography_report_v12.html'
    
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"\n✅ Report saved to: {output_path}")
    print("=" * 60)
    
    return output_path


if __name__ == '__main__':
    main()
