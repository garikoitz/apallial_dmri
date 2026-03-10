#!/usr/bin/env python3
"""
Generate Per-Tract Metrics Report (v2).

Computes for every passthrough-cleaned tract:
  1. Streamline count  (number of streamlines in cleaned .tck)
  2. SIFT2 weighted     (sum of wholebrain SIFT2 weights for matched raw streamlines)
  3. Density            (streamline_count / target_ROI_volume_mm³)

Source ROIs:
  - CTRL:     isocortex {anterior, central, posterior} per hemisphere
  - AP/APTR:  seg9 graft {anterior_dorsal, central_dorsal, posterior_dorsal} (shrunk) per hemisphere

Target ROIs (all groups):
  - {Striatum, Thalamus, Paleocortex/olfactory_cortex, Cerebellum} per hemisphere

Outputs:
  - tract_metrics_{version}.csv             — long-format with all metrics
  - tract_metrics_count_{version}.csv       — companion: streamline counts
  - tract_metrics_sift2_{version}.csv       — companion: SIFT2 weighted values
  - tract_metrics_density_{version}.csv     — companion: density values
  - tract_metrics_report_{version}.html     — tables + bar charts with error bars + individual points

Usage:
    python functions/generate_tract_metrics_report_v2.py --version v20250601
    python functions/generate_tract_metrics_report_v2.py --version v20250601 --derivatives-dir /path/to/derivs

Author: MBRAIN Project
"""

import os
import sys
import argparse
import subprocess
import csv
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np
import math

os.environ['MBRAIN_PIPELINE'] = 'flip1'

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config


# ============================================================
# Configuration
# ============================================================

EXCLUDED_SUBJECTS = {'sub-APtestmouseNmM'}

CTRL_SEGMENTS  = ["anterior", "central", "posterior"]
CTRL_TARGETS   = ["thalamus", "striatum", "cerebellum", "olfactory_cortex"]

APAPTR_SEGMENTS = ["anterior_dorsal", "central_dorsal", "posterior_dorsal"]
APAPTR_TARGETS  = ["Thalamus", "Striatum", "cerebellum", "Palecortex"]

SEGMENT_NORMALIZE = {
    "anterior": "anterior", "anterior_dorsal": "anterior",
    "central": "central",   "central_dorsal": "central",
    "posterior": "posterior", "posterior_dorsal": "posterior",
}

TARGET_NORMALIZE = {
    "thalamus": "thalamus",     "Thalamus": "thalamus",
    "striatum": "striatum",     "Striatum": "striatum",
    "cerebellum": "cerebellum",
    "olfactory_cortex": "paleocortex", "Palecortex": "paleocortex",
    "paleocortex": "paleocortex",
}

TARGET_DISPLAY = {
    "thalamus": "Thalamus",
    "striatum": "Striatum",
    "cerebellum": "Cerebellum",
    "paleocortex": "Paleocortex",
}


def get_paths():
    cfg = load_config()
    base = Path(cfg['paths']['derivatives_dir'])
    return cfg, base


def get_subjects(cfg, group=None, subject=None):
    subjects = []
    for group_name in ['CTRL', 'AP', 'APTR']:
        if group and group != group_name:
            continue
        for subj in cfg.get('groups', {}).get(group_name, []):
            sid = subj['id']
            if sid in EXCLUDED_SUBJECTS:
                continue
            if subject and subject != sid:
                continue
            subjects.append((sid, subj['session'], group_name))
    return subjects


def get_tract_info(group):
    if group == 'CTRL':
        return CTRL_SEGMENTS, CTRL_TARGETS
    return APAPTR_SEGMENTS, APAPTR_TARGETS


def get_tract_filename(group, hemi, segment, target):
    if group == 'CTRL':
        return f"isocortex_{hemi}_{segment}_to_{target}_{hemi}.tck"
    return f"seg9_graft_{hemi}_{segment}_to_{target}_{hemi}.tck"


def get_target_roi_path(group, base_dir, hemi, target):
    if group == 'CTRL':
        return base_dir / "atlas" / "rois" / "fod_space" / f"roi_{target}_{hemi}.nii.gz"
    return base_dir / "manual_roi" / "fod_space" / f"roi_{target}_{hemi}_flip_fod.mif"


# ============================================================
# Shell / MRtrix helpers
# ============================================================

def run_cmd(cmd, timeout=300):
    """Run a command, return (success, stdout, stderr)."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode == 0, r.stdout, r.stderr
    except Exception as e:
        return False, "", str(e)


def compute_roi_volume_mm3(roi_file):
    """Compute ROI volume in mm³ = n_nonzero_voxels × voxel_volume."""
    if not roi_file.exists():
        return 0.0
    # Get voxel spacing
    ok, stdout, _ = run_cmd(["mrinfo", "-spacing", str(roi_file)])
    if not ok:
        return 0.0
    try:
        spacings = [float(x) for x in stdout.strip().split()]
        voxel_vol = 1.0
        for s in spacings[:3]:
            voxel_vol *= s
    except Exception:
        return 0.0
    # Get non-zero voxel count
    ok, stdout, _ = run_cmd(["mrstats", "-mask", str(roi_file),
                              "-output", "count", str(roi_file)])
    if not ok:
        return 0.0
    try:
        n_voxels = int(float(stdout.strip().split()[0]))
    except Exception:
        return 0.0
    return n_voxels * voxel_vol


# ============================================================
# Fast .tck binary parsing for SIFT2 weight extraction
# ============================================================

def parse_tck_header(filepath):
    hdr = {}
    with open(filepath, 'rb') as f:
        while True:
            line = f.readline()
            line_str = line.decode('latin-1').strip()
            if line_str == 'END':
                break
            if ':' in line_str:
                key, val = line_str.split(':', 1)
                hdr[key.strip()] = val.strip()
        data_offset = f.tell()
    if 'file' in hdr:
        parts = hdr['file'].split()
        if len(parts) >= 2:
            data_offset = int(parts[-1])
    return data_offset, hdr


def count_tck_streamlines(tck_path):
    """Count streamlines in a .tck file via NaN separators."""
    data_offset, _ = parse_tck_header(str(tck_path))
    with open(str(tck_path), 'rb') as f:
        f.seek(data_offset)
        raw = np.fromfile(f, dtype='<f4')
    data = raw.reshape(-1, 3)
    return int(np.sum(np.isnan(data[:, 0])))


def build_wholebrain_hash_index(wholebrain_tck, n_hash_verts=3):
    """Build hash map: {first_N_vertices_bytes → streamline_index}."""
    import time
    t0 = time.time()
    data_offset, _ = parse_tck_header(str(wholebrain_tck))
    with open(str(wholebrain_tck), 'rb') as f:
        f.seek(data_offset)
        raw = np.fromfile(f, dtype='<f4')
    data = raw.reshape(-1, 3)
    nan_mask = np.isnan(data[:, 0])
    nan_indices = np.where(nan_mask)[0]
    inf_mask = np.isinf(data[:, 0])
    starts = np.concatenate([[0], nan_indices + 1])
    starts = starts[starts < len(data)]
    starts = starts[~inf_mask[starts]]
    end_idx = np.searchsorted(nan_indices, starts)
    end_pos = np.where(end_idx < len(nan_indices), nan_indices[end_idx], len(data))
    n_verts = end_pos - starts
    valid = starts[n_verts >= n_hash_verts]
    flat = data.ravel()
    n_floats = n_hash_verts * 3
    hash_map = {}
    for i, s in enumerate(valid):
        key = flat[s * 3: s * 3 + n_floats].tobytes()
        hash_map[key] = i
    elapsed = time.time() - t0
    print(f"      Hash index: {len(hash_map)} streamlines in {elapsed:.1f}s")
    return hash_map


def match_passthrough_weights(passthrough_tck, hash_map, sift2_weights, n_hash_verts=3):
    """Match passthrough streamlines to wholebrain and get SIFT2 weights.
    Returns (n_streamlines, sift2_sum)."""
    data_offset, _ = parse_tck_header(str(passthrough_tck))
    with open(str(passthrough_tck), 'rb') as f:
        f.seek(data_offset)
        raw = np.fromfile(f, dtype='<f4')
    data = raw.reshape(-1, 3)
    nan_mask = np.isnan(data[:, 0])
    nan_indices = np.where(nan_mask)[0]
    inf_mask = np.isinf(data[:, 0])
    starts = np.concatenate([[0], nan_indices + 1])
    starts = starts[starts < len(data)]
    starts = starts[~inf_mask[starts]]
    end_idx = np.searchsorted(nan_indices, starts)
    end_pos = np.where(end_idx < len(nan_indices), nan_indices[end_idx], len(data))
    n_verts = end_pos - starts
    valid = starts[n_verts >= n_hash_verts]
    flat = data.ravel()
    n_floats = n_hash_verts * 3
    matched_indices = []
    for s in valid:
        key = flat[s * 3: s * 3 + n_floats].tobytes()
        idx = hash_map.get(key)
        if idx is not None and idx < len(sift2_weights):
            matched_indices.append(idx)
    n_total = len(valid)
    if matched_indices:
        w = sift2_weights[np.array(matched_indices)]
        return n_total, float(w.sum())
    return n_total, 0.0


def load_sift2_weights(filepath):
    with open(str(filepath), 'r') as f:
        content = f.read()
    lines = [l for l in content.strip().split('\n') if not l.startswith('#')]
    return np.array([float(x) for x in ' '.join(lines).split()], dtype=np.float64)


# ============================================================
# Compute metrics for one subject
# ============================================================

def compute_subject_metrics(subject_id, session, group, deriv_dir):
    """Compute all tract metrics for one subject. Returns list of dicts."""
    base_dir = deriv_dir / subject_id / session
    dwi_dir = base_dir / "dwi"
    wb_dir = dwi_dir / "wholebrain"

    # AP/APTR use shrink tracts; CTRL uses non-shrink
    if group in ('AP', 'APTR'):
        tract_dir = wb_dir / "passthrough_tracts_shrink"
    else:
        tract_dir = wb_dir / "passthrough_tracts"

    raw_tract_dir = wb_dir / "passthrough_tracts"
    wholebrain_tck = wb_dir / "wholebrain_1000k.tck"
    sift2_file = wb_dir / "sift2_weights.txt"

    if not wholebrain_tck.exists():
        print(f"    WARNING: No wholebrain tractogram for {subject_id}")
        return []
    if not sift2_file.exists():
        print(f"    WARNING: No SIFT2 weights for {subject_id}")
        return []
    if not tract_dir.exists():
        print(f"    WARNING: No tract dir {tract_dir} for {subject_id}")
        return []

    segments, targets = get_tract_info(group)

    # Build hash index and load SIFT2 weights (once per subject)
    print(f"    Building wholebrain hash index for {subject_id}...")
    hash_map = build_wholebrain_hash_index(wholebrain_tck)
    print(f"    Loading SIFT2 weights...")
    sift2_weights = load_sift2_weights(sift2_file)
    print(f"      → {len(sift2_weights)} weights loaded")

    # Cache target ROI volumes
    roi_volumes = {}

    metrics_list = []
    for hemi in ['left', 'right']:
        for segment in segments:
            for target in targets:
                tract_base = get_tract_filename(group, hemi, segment, target).replace('.tck', '')
                tract_cleaned = tract_dir / f"{tract_base}_cleaned.tck"
                tract_raw = raw_tract_dir / f"{tract_base}.tck"

                # Streamline count from cleaned tract
                n_streamlines = 0
                if tract_cleaned.exists():
                    n_streamlines = count_tck_streamlines(tract_cleaned)
                elif tract_raw.exists():
                    n_streamlines = count_tck_streamlines(tract_raw)

                # SIFT2 from raw tract (cleaning resamples vertices, breaking hash)
                sift2_sum = 0.0
                sift2_source = tract_raw
                if sift2_source.exists():
                    _, sift2_sum = match_passthrough_weights(
                        sift2_source, hash_map, sift2_weights
                    )

                # Target ROI volume for density
                target_roi = get_target_roi_path(group, base_dir, hemi, target)
                roi_key = str(target_roi)
                if roi_key not in roi_volumes:
                    roi_volumes[roi_key] = compute_roi_volume_mm3(target_roi)
                target_vol = roi_volumes[roi_key]

                density = (n_streamlines / target_vol) if target_vol > 0 else 0.0

                seg_norm = SEGMENT_NORMALIZE.get(segment, segment)
                tgt_norm = TARGET_NORMALIZE.get(target, target.lower())

                metrics_list.append({
                    'subject': subject_id,
                    'group': group,
                    'hemisphere': hemi,
                    'segment': seg_norm,
                    'target': tgt_norm,
                    'streamline_count': n_streamlines,
                    'sift2_weighted': round(sift2_sum, 6),
                    'target_roi_volume_mm3': round(target_vol, 6),
                    'density': round(density, 6),
                    'tract_file': tract_base,
                })

    return metrics_list


# ============================================================
# CSV output
# ============================================================

FIELDNAMES = [
    'subject', 'group', 'hemisphere', 'segment', 'target',
    'streamline_count', 'sift2_weighted', 'target_roi_volume_mm3', 'density',
    'tract_file'
]


def write_main_csv(all_metrics, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        w.writerows(all_metrics)
    print(f"  CSV written: {path} ({len(all_metrics)} rows)")


def write_companion_csv(all_metrics, path, metric_key):
    """Write a wide-format companion CSV for one metric.
    Rows = subjects, columns = hemi_segment_target."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # Build pivot
    subjects_order = []
    seen = set()
    for m in all_metrics:
        key = (m['subject'], m['group'])
        if key not in seen:
            seen.add(key)
            subjects_order.append(key)

    # Column keys
    combos = []
    combos_seen = set()
    for m in all_metrics:
        ck = (m['hemisphere'], m['segment'], m['target'])
        if ck not in combos_seen:
            combos_seen.add(ck)
            combos.append(ck)
    combos.sort()

    col_names = [f"{h}_{s}_{t}" for h, s, t in combos]

    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['subject', 'group'] + col_names)
        for subj, grp in subjects_order:
            row = [subj, grp]
            for h, s, t in combos:
                match = [m for m in all_metrics
                         if m['subject'] == subj and m['hemisphere'] == h
                         and m['segment'] == s and m['target'] == t]
                if match:
                    row.append(match[0][metric_key])
                else:
                    row.append('')
            w.writerow(row)
    print(f"  Companion CSV: {path}")


# ============================================================
# HTML report generation
# ============================================================

def generate_html_report(all_metrics, output_path):
    """Generate HTML report with tables + bar charts + error bars + individual points.

    Layout: grouped by SEGMENT (anterior / central / posterior).
    X-axis of charts = target regions (Striatum, Thalamus, Paleocortex, Cerebellum).
    Bars = groups (CTRL, AP, APTR).
    """
    if not all_metrics:
        print("No metrics for report!")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Group data
    groups_data = defaultdict(list)
    for m in all_metrics:
        groups_data[m['group']].append(m)

    group_names = ['CTRL', 'AP', 'APTR']
    group_colors = {
        'CTRL': {'bg': 'rgba(46,204,113,0.5)', 'border': '#2ecc71'},
        'AP':   {'bg': 'rgba(231,76,60,0.5)',   'border': '#e74c3c'},
        'APTR': {'bg': 'rgba(243,156,18,0.5)',  'border': '#f39c12'},
    }

    segments_order = ['anterior', 'central', 'posterior']
    targets_order = ['striatum', 'thalamus', 'paleocortex', 'cerebellum']

    metrics_defs = [
        ('streamline_count', 'Streamline Count', 'Count'),
        ('sift2_weighted',   'SIFT2 Weighted',   'SIFT2 Sum'),
        ('density',          'Density',           'Streamlines / mm³'),
    ]

    # ── HTML header ──
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MBRAIN Tract Metrics Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       margin: 0; padding: 20px; background: #f5f5f5; color: #333; }}
.container {{ max-width: 1600px; margin: 0 auto; }}
h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
h2 {{ color: #34495e; margin-top: 30px; }}
h3 {{ color: #7f8c8d; margin-top: 20px; }}
.info-box {{ background: #fdf6e3; border-left: 4px solid #f39c12; padding: 12px 18px;
             margin: 10px 0; font-size: 13px; }}
.info-box code {{ background: #eee; padding: 1px 5px; border-radius: 3px; }}
.summary-box {{ background: white; border-radius: 8px; padding: 20px;
                margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 12px; }}
th, td {{ border: 1px solid #ddd; padding: 5px 8px; text-align: right; }}
th {{ background: #3498db; color: white; position: sticky; top: 0; font-size: 11px; }}
td:first-child, td:nth-child(2) {{ text-align: left; }}
tr:nth-child(even) {{ background: #f8f9fa; }}
tr:hover {{ background: #e8f4fd; }}
tr.group-mean {{ font-weight: bold; background: #eee !important; }}
.chart-grid {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 15px; margin: 15px 0; }}
.chart-box {{ background: white; border-radius: 8px; padding: 15px;
              box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
.tabs {{ display: flex; gap: 5px; margin: 20px 0 10px; }}
.tab {{ padding: 8px 16px; background: #ecf0f1; border: none; border-radius: 5px 5px 0 0;
        cursor: pointer; font-size: 14px; }}
.tab.active {{ background: #3498db; color: white; }}
.tab-content {{ display: none; }} .tab-content.active {{ display: block; }}
.group-ctrl {{ color: #2ecc71; }} .group-ap {{ color: #e74c3c; }} .group-aptr {{ color: #f39c12; }}
.zero {{ color: #ccc; }}
.footer {{ text-align: center; color: #95a5a6; margin-top: 40px; padding: 20px; font-size: 12px; }}
</style>
</head>
<body>
<div class="container">
<h1>MBRAIN Tract Metrics Report</h1>
<p>Generated: {timestamp}</p>

<div class="info-box">
<strong>Segments (source ROIs):</strong>
CTRL: <code>isocortex_{{anterior, central, posterior}}</code> &harr;
AP/APTR: <code>seg9_graft_{{anterior, central, posterior}}_dorsal</code> (shrunk ROIs)<br>
<strong>Targets (X-axis):</strong> Striatum, Thalamus, Paleocortex, Cerebellum (L/R)<br>
<strong>Tracts:</strong> AP/APTR from <code>passthrough_tracts_shrink/</code>;
CTRL from <code>passthrough_tracts/</code>. All MATLAB AFQ-cleaned.
</div>
<div class="info-box" style="border-left-color: #3498db;">
<strong>Metrics:</strong><br>
<b>Count</b>: Number of streamlines in cleaned passthrough tract.<br>
<b>SIFT2 Weighted</b>: Sum of SIFT2 weights for matched streamlines (whole-brain SIFT2).
Biologically-corrected connection strength.<br>
<b>Density</b>: Streamline count / target ROI volume (mm&sup3;).
Normalizes for target region size.
</div>
"""

    # ── Summary ──
    html += '<div class="summary-box"><h2>Summary</h2>\n'
    for gn in group_names:
        gd = groups_data.get(gn, [])
        if not gd:
            continue
        subjs = sorted(set(m['subject'] for m in gd))
        total_n = sum(m['streamline_count'] for m in gd)
        total_sift2 = sum(m['sift2_weighted'] for m in gd)
        html += f'<h3 class="group-{gn.lower()}">{gn} ({len(subjs)} subjects)</h3>\n'
        html += f'<p>Subjects: {", ".join(subjs)}<br>'
        html += f'Total streamlines: {total_n:,} | Total SIFT2: {total_sift2:.2f}</p>\n'
    html += '</div>\n'

    # ── Tabs: Left / Right ──
    html += '<div class="tabs">\n'
    html += '<button class="tab active" onclick="showTab(\'left\')">Left Hemisphere</button>\n'
    html += '<button class="tab" onclick="showTab(\'right\')">Right Hemisphere</button>\n'
    html += '</div>\n'

    chart_counter = 0

    for hemi in ['left', 'right']:
        active = ' active' if hemi == 'left' else ''
        html += f'<div id="tab-{hemi}" class="tab-content{active}">\n'
        html += f'<h2>{hemi.capitalize()} Hemisphere</h2>\n'

        # ──── Sections grouped by SEGMENT ────
        for segment in segments_order:
            seg_disp = segment.capitalize()
            html += f'<div class="summary-box"><h3>Source Segment: {seg_disp}</h3>\n'

            # ──── Table: columns = targets ────
            html += '<table><thead><tr>'
            html += '<th>Subject</th><th>Group</th>'
            for tgt in targets_order:
                tgt_disp = TARGET_DISPLAY.get(tgt, tgt.capitalize())
                html += f'<th colspan="3">{tgt_disp}</th>'
            html += '</tr><tr><th></th><th></th>'
            for _ in targets_order:
                html += '<th>Count</th><th>SIFT2</th><th>Density</th>'
            html += '</tr></thead><tbody>\n'

            subj_list = sorted(set(
                (m['subject'], m['group']) for m in all_metrics
                if m['hemisphere'] == hemi and m['segment'] == segment
            ), key=lambda x: ({'CTRL': 0, 'AP': 1, 'APTR': 2}.get(x[1], 9), x[0]))

            for subj_id, grp in subj_list:
                gc = f'group-{grp.lower()}'
                html += f'<tr><td>{subj_id}</td><td class="{gc}">{grp}</td>'
                for tgt in targets_order:
                    match = [m for m in all_metrics
                             if m['subject'] == subj_id and m['hemisphere'] == hemi
                             and m['segment'] == segment and m['target'] == tgt]
                    if match:
                        mm = match[0]
                        cnt = mm['streamline_count']
                        s2 = mm['sift2_weighted']
                        dn = mm['density']
                        zc = ' class="zero"' if cnt == 0 else ''
                        html += f'<td{zc}>{cnt:,}</td><td>{s2:.4f}</td><td>{dn:.2f}</td>'
                    else:
                        html += '<td>-</td><td>-</td><td>-</td>'
                html += '</tr>\n'

            # Group mean rows
            for gn in group_names:
                g_subjs = [s for s, g in subj_list if g == gn]
                if not g_subjs:
                    continue
                gc = f'group-{gn.lower()}'
                html += f'<tr class="group-mean"><td>{gn} Mean &plusmn; SEM</td><td class="{gc}">{gn}</td>'
                for tgt in targets_order:
                    g_m = [m for m in all_metrics
                           if m['group'] == gn and m['hemisphere'] == hemi
                           and m['segment'] == segment and m['target'] == tgt]
                    if g_m:
                        vals_c = [m['streamline_count'] for m in g_m]
                        vals_s = [m['sift2_weighted'] for m in g_m]
                        vals_d = [m['density'] for m in g_m]
                        n = len(vals_c)
                        sem_c = np.std(vals_c, ddof=1) / math.sqrt(n) if n > 1 else 0
                        sem_s = np.std(vals_s, ddof=1) / math.sqrt(n) if n > 1 else 0
                        sem_d = np.std(vals_d, ddof=1) / math.sqrt(n) if n > 1 else 0
                        html += (f'<td>{np.mean(vals_c):.0f}&plusmn;{sem_c:.0f}</td>'
                                 f'<td>{np.mean(vals_s):.4f}&plusmn;{sem_s:.4f}</td>'
                                 f'<td>{np.mean(vals_d):.2f}&plusmn;{sem_d:.2f}</td>')
                    else:
                        html += '<td>-</td><td>-</td><td>-</td>'
                html += '</tr>\n'

            html += '</tbody></table>\n'

            # ──── Charts: 3 bar charts (Count, SIFT2, Density) ────
            # X-axis = targets, bars = groups (CTRL, AP, APTR)
            html += '<div class="chart-grid">\n'
            for mk, mtitle, ylabel in metrics_defs:
                cid = f'chart_{chart_counter}'
                chart_counter += 1
                html += f'<div class="chart-box"><canvas id="{cid}"></canvas></div>\n'
            html += '</div>\n'  # chart-grid

            html += '<script>\n'

            chart_counter -= 3  # rewind to match canvas IDs
            for mk, mtitle, ylabel in metrics_defs:
                cid = f'chart_{chart_counter}'
                chart_counter += 1

                bar_datasets = []
                for gn in group_names:
                    gd = groups_data.get(gn, [])
                    if not gd:
                        continue
                    g_subjs = sorted(set(m['subject'] for m in gd))
                    means = []
                    sems = []
                    all_points = []

                    for tgt in targets_order:
                        vals = []
                        for subj in g_subjs:
                            match = [m for m in gd
                                     if m['subject'] == subj and m['hemisphere'] == hemi
                                     and m['segment'] == segment and m['target'] == tgt]
                            if match:
                                v = match[0][mk]
                                if v is not None:
                                    vals.append(float(v))
                        if vals:
                            means.append(round(float(np.mean(vals)), 6))
                            sem = float(np.std(vals, ddof=1) / math.sqrt(len(vals))) if len(vals) > 1 else 0.0
                            sems.append(round(sem, 6))
                            all_points.append([round(v, 6) for v in vals])
                        else:
                            means.append(0)
                            sems.append(0)
                            all_points.append([])

                    colors = group_colors[gn]
                    bar_datasets.append({
                        'label': gn,
                        'data': means,
                        'backgroundColor': colors['bg'],
                        'borderColor': colors['border'],
                        'borderWidth': 1.5,
                        '_errorBars': sems,
                        '_dataPoints': all_points,
                    })

                chart_title = f'{mtitle}: {seg_disp} segment ({hemi.capitalize()})'
                target_labels = json.dumps([TARGET_DISPLAY.get(t, t.capitalize()) for t in targets_order])

                html += f"""
(function() {{
  var ctx = document.getElementById('{cid}');
  new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels: {target_labels},
      datasets: {json.dumps(bar_datasets)}
    }},
    options: {{
      responsive: true,
      plugins: {{
        title: {{ display: true, text: '{chart_title}', font: {{size: 13}} }},
        legend: {{ position: 'top', labels: {{font: {{size: 11}}}} }}
      }},
      scales: {{
        y: {{
          beginAtZero: true,
          title: {{ display: true, text: '{ylabel}', font: {{size: 11}} }}
        }}
      }}
    }},
    plugins: [errorBarsPlugin, dataPointsPlugin]
  }});
}})();
"""

            html += '</script>\n'
            html += '</div>\n'  # summary-box

        html += '</div>\n'  # tab-content

    # ── Chart.js plugins ──
    plugin_js = """
<script>
// ── Custom Chart.js plugin: Error Bars ──
const errorBarsPlugin = {
  id: 'customErrorBars',
  afterDatasetsDraw(chart) {
    const ctx = chart.ctx;
    chart.data.datasets.forEach((ds, di) => {
      if (!ds._errorBars) return;
      const meta = chart.getDatasetMeta(di);
      if (meta.hidden) return;
      meta.data.forEach((bar, idx) => {
        const err = ds._errorBars[idx];
        if (err == null || err === 0 || isNaN(err)) return;
        const x = bar.x;
        const yScale = chart.scales.y;
        const mean = ds.data[idx];
        const yTop = yScale.getPixelForValue(mean + err);
        const yBottom = yScale.getPixelForValue(Math.max(0, mean - err));
        const capW = bar.width ? bar.width * 0.3 : 5;
        ctx.save();
        ctx.strokeStyle = 'rgba(0,0,0,0.7)';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(x, yTop);
        ctx.lineTo(x, yBottom);
        ctx.moveTo(x - capW, yTop);
        ctx.lineTo(x + capW, yTop);
        ctx.moveTo(x - capW, yBottom);
        ctx.lineTo(x + capW, yBottom);
        ctx.stroke();
        ctx.restore();
      });
    });
  }
};

// ── Custom Chart.js plugin: Individual Data Points ──
const dataPointsPlugin = {
  id: 'customDataPoints',
  afterDatasetsDraw(chart) {
    const ctx = chart.ctx;
    chart.data.datasets.forEach((ds, di) => {
      if (!ds._dataPoints) return;
      const meta = chart.getDatasetMeta(di);
      if (meta.hidden) return;
      meta.data.forEach((bar, idx) => {
        const pts = ds._dataPoints[idx];
        if (!pts || pts.length === 0) return;
        const x = bar.x;
        const yScale = chart.scales.y;
        const n = pts.length;
        const spread = Math.min(bar.width ? bar.width * 0.6 : 8, 16);
        const jitter = n > 1
          ? pts.map((_, i) => (i - (n-1)/2) * (spread / Math.max(n-1, 1)))
          : [0];
        ctx.save();
        ctx.fillStyle = ds.borderColor || '#333';
        ctx.globalAlpha = 0.7;
        pts.forEach((val, pi) => {
          const px = x + jitter[pi];
          const py = yScale.getPixelForValue(val);
          ctx.beginPath();
          ctx.arc(px, py, 3.5, 0, 2 * Math.PI);
          ctx.fill();
        });
        ctx.restore();
      });
    });
  }
};
</script>
"""

    tab_js = """
<script>
function showTab(hemi) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
  document.getElementById('tab-' + hemi).classList.add('active');
  event.target.classList.add('active');
}
</script>
"""

    footer = f"""
<div class="footer">
<p>MBRAIN Tract Metrics Report | Generated {timestamp}</p>
<p>Metrics: Streamline Count, SIFT2 Weighted, Density (Count / Target ROI Volume)</p>
</div>
</div>
</body>
</html>
"""

    final_html = html.replace(
        '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>',
        '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>\n' + plugin_js
    )
    final_html += tab_js + footer

    with open(output_path, 'w') as f:
        f.write(final_html)
    print(f"  HTML report: {output_path}")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate tract metrics report (v2): Count + SIFT2 + Density',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--group', choices=['CTRL', 'AP', 'APTR'],
                        help='Process specific group only')
    parser.add_argument('--subject', type=str,
                        help='Process specific subject only')
    parser.add_argument('--output-dir', type=str, default='reports',
                        help='Output directory (default: reports/)')
    parser.add_argument('--version', type=str, default='v01',
                        help='Report version suffix')
    parser.add_argument('--derivatives-dir', type=str,
                        help='Override derivatives directory')

    args = parser.parse_args()

    if args.derivatives_dir:
        os.environ['MBRAIN_DERIVATIVES_DIR'] = args.derivatives_dir

    print("=" * 70)
    print("MBRAIN TRACT METRICS REPORT (v2)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    cfg, deriv_dir = get_paths()
    subjects = get_subjects(cfg, group=args.group, subject=args.subject)

    if not subjects:
        print("ERROR: No subjects found!")
        sys.exit(1)

    print(f"\nProcessing {len(subjects)} subjects:")
    for sid, ses, grp in subjects:
        print(f"  {sid} ({grp})")

    # Process all subjects
    all_metrics = []
    for i, (sid, ses, grp) in enumerate(subjects, 1):
        print(f"\n{'=' * 60}")
        print(f"[{i}/{len(subjects)}] {sid} ({grp})")
        print(f"{'=' * 60}")
        metrics = compute_subject_metrics(sid, ses, grp, deriv_dir)
        all_metrics.extend(metrics)
        print(f"  → {len(metrics)} tract metrics computed")

    # Write outputs
    out_dir = Path(args.output_dir)
    ver = args.version

    # Main CSV
    csv_main = out_dir / f"tract_metrics_{ver}.csv"
    write_main_csv(all_metrics, csv_main)

    # Companion CSVs (wide format)
    write_companion_csv(all_metrics, out_dir / f"tract_metrics_count_{ver}.csv", 'streamline_count')
    write_companion_csv(all_metrics, out_dir / f"tract_metrics_sift2_{ver}.csv", 'sift2_weighted')
    write_companion_csv(all_metrics, out_dir / f"tract_metrics_density_{ver}.csv", 'density')

    # HTML report
    html_path = out_dir / f"tract_metrics_report_{ver}.html"
    generate_html_report(all_metrics, html_path)

    print(f"\n{'=' * 70}")
    print(f"REPORT COMPLETE")
    print(f"  Main CSV:     {csv_main}")
    print(f"  Count CSV:    {out_dir / f'tract_metrics_count_{ver}.csv'}")
    print(f"  SIFT2 CSV:    {out_dir / f'tract_metrics_sift2_{ver}.csv'}")
    print(f"  Density CSV:  {out_dir / f'tract_metrics_density_{ver}.csv'}")
    print(f"  HTML:         {html_path}")
    print(f"  Total metrics: {len(all_metrics)}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
