#!/usr/bin/env python3
"""
Generate Per-Tract Metrics Report: SIFT2 connectivity + Fiber Density (AFD).

Computes for every passthrough tract (segment → target, per hemisphere):
  1. Streamline count
  2. SIFT2 connectivity (sum of SIFT2 weights for tract streamlines)
  3. SIFT2 density (connectivity / seed ROI voxels × 1000)
  4. Mean Fiber Density (AFD sampled along streamlines via fixel→voxel→tcksample)

Outputs:
  - CSV with all per-tract metrics for all subjects
  - HTML report with group comparison tables and bar charts

Usage:
    python functions/generate_tract_metrics_report.py
    python functions/generate_tract_metrics_report.py --group CTRL
    python functions/generate_tract_metrics_report.py --subject sub-129-921-R1F
    python functions/generate_tract_metrics_report.py --skip-extraction  # reuse existing weight files

Author: MBRAIN Project
"""

import os
import sys
os.environ['MBRAIN_PIPELINE'] = 'flip1'
import argparse
import subprocess
import csv
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))
from functions.config_loader import load_config


# ============================================================
# Configuration
# ============================================================

def get_paths():
    """Get paths from config."""
    cfg = load_config()
    base = Path(cfg['paths']['derivatives_dir'])
    return cfg, base


def get_subjects(cfg, group=None, subject=None):
    """Get list of (subject_id, session, group) tuples."""
    # APtestmouseNmM is a test subject borrowed from 769-C-OO-NmM — always excluded
    EXCLUDED_SUBJECTS = {'sub-APtestmouseNmM'}
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


# ============================================================
# Tract definitions per group
# ============================================================
# 3-segment mapping for cross-group comparison:
#   CTRL anterior      <->  AP/APTR anterior_dorsal   -> "anterior"
#   CTRL central       <->  AP/APTR central_dorsal    -> "central"
#   CTRL posterior     <->  AP/APTR posterior_dorsal   -> "posterior"
# (CTRL anterior_central and posterior_central are excluded)
#
# Target mapping:
#   CTRL olfactory_cortex  <->  AP/APTR Palecortex  -> "olfactory_cortex"
#   CTRL thalamus          <->  AP/APTR Thalamus    -> "thalamus"
#   CTRL striatum          <->  AP/APTR Striatum    -> "striatum"
#   CTRL cerebellum        <->  AP/APTR cerebellum  -> "cerebellum"

CTRL_SEGMENTS = ["anterior", "central", "posterior"]
CTRL_TARGETS = ["thalamus", "striatum", "cerebellum", "olfactory_cortex"]

APAPTR_SEGMENTS = ["anterior_dorsal", "central_dorsal", "posterior_dorsal"]
APAPTR_TARGETS = ["Thalamus", "Striatum", "cerebellum", "Palecortex"]

# Normalized segment names for cross-group comparison
SEGMENT_NORMALIZE = {
    "anterior": "anterior",
    "anterior_dorsal": "anterior",
    "central": "central",
    "central_dorsal": "central",
    "posterior": "posterior",
    "posterior_dorsal": "posterior",
}

# Normalized target names for cross-group comparison
TARGET_NORMALIZE = {
    "thalamus": "thalamus",
    "Thalamus": "thalamus",
    "striatum": "striatum",
    "Striatum": "striatum",
    "cerebellum": "cerebellum",
    "olfactory_cortex": "olfactory_cortex",
    "Palecortex": "olfactory_cortex",
    "paleocortex": "olfactory_cortex",
}


def get_tract_info(group):
    """Return (segments, targets) for a group."""
    if group == 'CTRL':
        return CTRL_SEGMENTS, CTRL_TARGETS
    else:
        return APAPTR_SEGMENTS, APAPTR_TARGETS


def get_tract_filename(group, hemi, segment, target):
    """Get the passthrough tract .tck filename."""
    if group == 'CTRL':
        return f"isocortex_{hemi}_{segment}_to_{target}_{hemi}.tck"
    else:
        return f"seg9_graft_{hemi}_{segment}_to_{target}_{hemi}.tck"


def get_seed_roi_path(group, base_dir, hemi, segment, use_shrink=False):
    """Get the seed ROI path for a segment."""
    if group == 'CTRL':
        return base_dir / "dwi" / "isocortex_segments_v2" / f"isocortex_{hemi}_{segment}.nii.gz"
    else:
        if use_shrink:
            return base_dir / "dwi" / "graft_segments_fod" / "shrink" / f"seg9_graft_{hemi}_{segment}_fod_shrink.mif"
        return base_dir / "dwi" / "graft_segments_fod" / f"seg9_graft_{hemi}_{segment}_fod.mif"


def get_target_roi_path(group, base_dir, hemi, target):
    """Get the target ROI path."""
    if group == 'CTRL':
        return base_dir / "atlas" / "rois" / "fod_space" / f"roi_{target}_{hemi}.nii.gz"
    else:
        return base_dir / "manual_roi" / "fod_space" / f"roi_{target}_{hemi}_flip_fod.mif"


# ============================================================
# Fast .tck binary parsing for SIFT2 weight extraction
# ============================================================

def parse_tck_header(filepath):
    """Parse .tck file header, return data_offset and header dict."""
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
    """Count the number of streamlines in a .tck file by counting NaN separators."""
    data_offset, _ = parse_tck_header(str(tck_path))
    with open(str(tck_path), 'rb') as f:
        f.seek(data_offset)
        raw_data = np.fromfile(f, dtype='<f4')
    data = raw_data.reshape(-1, 3)
    nan_count = np.sum(np.isnan(data[:, 0]))
    # Subtract 1 for the Inf terminator if present, but NaN count = streamline count
    return int(nan_count)


def build_wholebrain_hash_index(wholebrain_tck, n_hash_verts=3):
    """
    Read the wholebrain .tck file once and build a hash map:
    {first_N_vertices_bytes: streamline_index}.
    ~4s for 1M streamlines (4GB file).
    """
    import time
    t0 = time.time()
    data_offset, _ = parse_tck_header(str(wholebrain_tck))

    with open(str(wholebrain_tck), 'rb') as f:
        f.seek(data_offset)
        raw_data = np.fromfile(f, dtype='<f4')

    data = raw_data.reshape(-1, 3)
    nan_mask = np.isnan(data[:, 0])
    nan_indices = np.where(nan_mask)[0]
    inf_mask = np.isinf(data[:, 0])

    starts = np.concatenate([[0], nan_indices + 1])
    starts = starts[starts < len(data)]
    starts = starts[~inf_mask[starts]]

    end_idx = np.searchsorted(nan_indices, starts)
    end_positions = np.where(end_idx < len(nan_indices), nan_indices[end_idx], len(data))
    n_verts = end_positions - starts
    valid_starts = starts[n_verts >= n_hash_verts]

    flat = data.ravel()
    n_floats = n_hash_verts * 3
    hash_map = {}
    for i, s in enumerate(valid_starts):
        key = flat[s * 3: s * 3 + n_floats].tobytes()
        hash_map[key] = i

    elapsed = time.time() - t0
    print(f"    Wholebrain hash index: {len(hash_map)} streamlines in {elapsed:.1f}s")
    return hash_map


def match_passthrough_weights(passthrough_tck, hash_map, sift2_weights, n_hash_verts=3):
    """
    Match passthrough tract streamlines to the wholebrain hash index
    and extract their SIFT2 weights.
    Returns (n_streamlines, sift2_sum, sift2_mean, matched_weights_array).
    """
    data_offset, _ = parse_tck_header(str(passthrough_tck))

    with open(str(passthrough_tck), 'rb') as f:
        f.seek(data_offset)
        raw_data = np.fromfile(f, dtype='<f4')

    data = raw_data.reshape(-1, 3)
    nan_mask = np.isnan(data[:, 0])
    nan_indices = np.where(nan_mask)[0]
    inf_mask = np.isinf(data[:, 0])

    starts = np.concatenate([[0], nan_indices + 1])
    starts = starts[starts < len(data)]
    starts = starts[~inf_mask[starts]]

    end_idx = np.searchsorted(nan_indices, starts)
    end_positions = np.where(end_idx < len(nan_indices), nan_indices[end_idx], len(data))
    n_verts = end_positions - starts
    valid_starts = starts[n_verts >= n_hash_verts]

    flat = data.ravel()
    n_floats = n_hash_verts * 3

    matched_indices = []
    for s in valid_starts:
        key = flat[s * 3: s * 3 + n_floats].tobytes()
        idx = hash_map.get(key)
        if idx is not None:
            matched_indices.append(idx)

    n_total = len(valid_starts)
    n_matched = len(matched_indices)

    if n_matched > 0:
        tract_weights = sift2_weights[matched_indices]
        return n_total, float(tract_weights.sum()), float(tract_weights.mean()), tract_weights
    return n_total, 0.0, 0.0, np.array([])


def load_sift2_weights(filepath):
    """Load SIFT2 weights from text file. Returns numpy array."""
    with open(str(filepath), 'r') as f:
        content = f.read()
    lines = [l for l in content.strip().split('\n') if not l.startswith('#')]
    return np.array([float(x) for x in ' '.join(lines).split()], dtype=np.float64)


# ============================================================
# Step 1: Prepare voxel AFD map (fixel2voxel)
# ============================================================

def prepare_afd_voxel_map(dwi_dir):
    """Convert fixel AFD to voxel-space sum. Returns path to afd_voxel_sum.mif."""
    fixels_dir = dwi_dir / "fixels_group"
    afd_file = fixels_dir / "afd.mif"
    output_file = dwi_dir / "afd_voxel_sum.mif"

    if not afd_file.exists():
        print(f"    WARNING: No fixel AFD found: {afd_file}")
        return None

    if output_file.exists():
        return output_file

    ok, _, stderr = run_cmd(
        ["fixel2voxel", str(afd_file), "sum", str(output_file), "-force"],
        "Creating voxel AFD map (fixel2voxel sum)"
    )
    if ok:
        return output_file
    print(f"    ERROR: fixel2voxel failed: {stderr[:200]}")
    return None


# ============================================================
# MRtrix command helpers
# ============================================================

def run_cmd(cmd, description=None, timeout=600):
    """Run a shell command, return (success, stdout, stderr)."""
    if description:
        print(f"    {description}...", end=" ", flush=True)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if result.returncode != 0:
            if description:
                print(f"FAILED")
            return False, result.stdout, result.stderr
        if description:
            print("OK")
        return True, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        if description:
            print("TIMEOUT")
        return False, "", "Timeout"
    except Exception as e:
        if description:
            print(f"ERROR: {e}")
        return False, "", str(e)


def count_roi_voxels(roi_file):
    """Count non-zero voxels in an ROI."""
    if not roi_file.exists():
        return 0
    suffix = ''.join(roi_file.suffixes).lower()
    if '.mif' in suffix:
        ok, stdout, _ = run_cmd(["mrstats", "-mask", str(roi_file), "-output", "count", str(roi_file)])
        if ok:
            return int(float(stdout.strip()))
    else:
        try:
            import nibabel as nib
            data = nib.load(str(roi_file)).get_fdata()
            return int(np.sum(data > 0))
        except Exception:
            pass
    return 0


def sample_tract_afd(tract_tck, afd_voxel_map, output_txt):
    """
    Sample mean AFD along each streamline using tcksample.
    Returns mean AFD across all streamlines, or None on failure.
    """
    if not tract_tck.exists() or afd_voxel_map is None:
        return None

    ok, _, _ = run_cmd([
        "tcksample",
        str(tract_tck),
        str(afd_voxel_map),
        str(output_txt),
        "-stat_tck", "mean",
        "-force"
    ], timeout=120)

    if not ok or not output_txt.exists():
        return None

    try:
        with open(output_txt, 'r') as f:
            values = []
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                values.extend([float(x) for x in line.split()])
            if values:
                return float(np.mean(values))
    except Exception:
        pass
    return None


# ============================================================
# Compute all metrics for one subject
# ============================================================

def filter_cerebellum_via_thalamus(tract_tck, thalamus_roi, output_tck):
    """Filter a cerebellum passthrough tract to only keep streamlines passing through thalamus.
    
    Uses tckedit -include to apply the thalamus ROI as additional constraint.
    Returns path to filtered .tck (or None if filtering fails).
    """
    import subprocess
    if not tract_tck.exists():
        return None
    if not thalamus_roi.exists():
        print(f"      WARNING: Thalamus ROI not found: {thalamus_roi}")
        return None
    
    cmd = [
        'tckedit', str(tract_tck),
        str(output_tck),
        '-include', str(thalamus_roi),
        '-force'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"      ERROR tckedit: {result.stderr.strip()}")
        return None
    if output_tck.exists():
        return output_tck
    return None


def compute_subject_metrics(subject_id, session, group, deriv_dir, skip_extraction=False,
                            use_shrink=False, use_cleaned=False, cerebellum_via_thalamus=False):
    """Compute all tract metrics for one subject. Returns list of metric dicts.
    
    Args:
        use_cleaned: If True, prefer _cleaned.tck files over raw .tck.
                     Falls back to raw .tck if _cleaned doesn't exist.
        cerebellum_via_thalamus: If True, for cerebellum target only, filter
                                tracts to keep only those passing through thalamus.
    """
    import time

    base_dir = deriv_dir / subject_id / session
    dwi_dir = base_dir / "dwi"
    wholebrain_dir = dwi_dir / "wholebrain"
    # AP/APTR ALWAYS use shrink tracts (primary analysis); CTRL uses non-shrink
    if group in ('AP', 'APTR'):
        tract_dir = wholebrain_dir / "passthrough_tracts_shrink"
        # Also keep raw tract dir for SIFT2 matching fallback
        raw_tract_dir = wholebrain_dir / "passthrough_tracts"
    else:
        tract_dir = wholebrain_dir / "passthrough_tracts"
        raw_tract_dir = tract_dir
    wholebrain_tck = wholebrain_dir / "wholebrain_1000k.tck"
    sift2_weights_file = wholebrain_dir / "sift2_weights.txt"

    if not wholebrain_tck.exists():
        print(f"    ERROR: No wholebrain tractogram")
        return []
    if not sift2_weights_file.exists():
        print(f"    ERROR: No SIFT2 weights")
        return []

    segments, targets = get_tract_info(group)
    metrics_list = []

    # Step 1: Prepare voxel AFD map
    afd_map = prepare_afd_voxel_map(dwi_dir)

    # Step 2: Build wholebrain hash index (ONE read of the 4GB file)
    print(f"    Building wholebrain hash index...")
    hash_map = build_wholebrain_hash_index(wholebrain_tck)

    # Step 3: Load SIFT2 weights
    print(f"    Loading SIFT2 weights...")
    sift2_weights = load_sift2_weights(sift2_weights_file)
    print(f"    → {len(sift2_weights)} weights loaded")

    # Step 4: Process all tracts
    for hemi in ['left', 'right']:
        afd_dir = tract_dir / "afd_samples"
        afd_dir.mkdir(parents=True, exist_ok=True)

        for segment in segments:
            for target in targets:
                tract_name = get_tract_filename(group, hemi, segment, target).replace('.tck', '')
                # Prefer _cleaned.tck if --cleaned and file exists
                tract_tck_cleaned = tract_dir / f"{tract_name}_cleaned.tck"
                tract_tck_raw = tract_dir / f"{tract_name}.tck"
                if use_cleaned and tract_tck_cleaned.exists():
                    tract_tck = tract_tck_cleaned
                    used_cleaned = True
                else:
                    tract_tck = tract_tck_raw
                    used_cleaned = False

                # Cerebellum-via-thalamus: filter cerebellum tracts through thalamus ROI
                target_normalized_check = TARGET_NORMALIZE.get(target, target.lower())
                if cerebellum_via_thalamus and target_normalized_check == 'cerebellum' and tract_tck.exists():
                    # Get the thalamus ROI for this group/hemi
                    thal_target = 'thalamus' if group == 'CTRL' else 'Thalamus'
                    thalamus_roi = get_target_roi_path(group, base_dir, hemi, thal_target)
                    suffix = '_cleaned_via_thalamus.tck' if used_cleaned else '_via_thalamus.tck'
                    via_thal_tck = tract_dir / f"{tract_name}{suffix}"
                    filtered = filter_cerebellum_via_thalamus(tract_tck, thalamus_roi, via_thal_tck)
                    if filtered:
                        tract_tck = filtered

                # SIFT2 via hash matching (fast!)
                # IMPORTANT: Always use RAW tract for SIFT2 matching, because
                # MATLAB AFQ cleaning resamples streamlines (changes vertex coords),
                # which breaks the hash match to wholebrain streamlines.
                # The streamline count reported is from the file actually used
                # (cleaned if available, raw otherwise).
                n_streamlines = 0
                sift2_sum = 0.0
                sift2_mean = 0.0
                sift2_file = tract_tck_raw  # Always match SIFT2 from raw tract
                if sift2_file.exists():
                    _, sift2_sum, sift2_mean, matched_weights = match_passthrough_weights(
                        sift2_file, hash_map, sift2_weights
                    )
                if tract_tck.exists():
                    # Get streamline count from the tract we're actually reporting
                    n_streamlines, _, _, _ = match_passthrough_weights(
                        tract_tck, hash_map, sift2_weights
                    )
                    # If cleaned file doesn't hash-match (resampled), count streamlines directly
                    if n_streamlines == 0 and used_cleaned:
                        n_streamlines = count_tck_streamlines(tract_tck)

                # Seed ROI voxels for density normalization
                seed_roi = get_seed_roi_path(group, base_dir, hemi, segment,
                                             use_shrink=use_shrink)
                roi_voxels = count_roi_voxels(seed_roi)
                sift2_density = (sift2_sum / roi_voxels * 1000) if roi_voxels > 0 else 0.0

                # AFD sampling (tcksample on the small passthrough file — fast)
                afd_txt = afd_dir / f"{tract_name}_afd.txt"
                mean_afd = None
                if afd_map and tract_tck.exists() and n_streamlines > 0:
                    mean_afd = sample_tract_afd(tract_tck, afd_map, afd_txt)

                # Normalize segment/target names for cross-group comparison
                segment_normalized = SEGMENT_NORMALIZE.get(segment, segment)
                target_normalized = TARGET_NORMALIZE.get(target, target.lower())

                metrics_list.append({
                    'subject': subject_id,
                    'group': group,
                    'hemisphere': hemi,
                    'segment': segment_normalized,
                    'target': target_normalized,
                    'streamline_count': n_streamlines,
                    'sift2_sum': round(sift2_sum, 4),
                    'sift2_mean': round(sift2_mean, 6),
                    'sift2_density': round(sift2_density, 4),
                    'roi_voxels': roi_voxels,
                    'mean_afd': round(mean_afd, 6) if mean_afd is not None else None,
                    'tract_name': tract_name,
                    'cleaned': used_cleaned if use_cleaned else None,
                })

    return metrics_list


# ============================================================
# CSV output
# ============================================================

def write_csv(all_metrics, output_path):
    """Write all metrics to CSV."""
    if not all_metrics:
        print("No metrics to write!")
        return

    fieldnames = [
        'subject', 'group', 'hemisphere', 'segment', 'target',
        'streamline_count', 'sift2_sum', 'sift2_mean', 'sift2_density',
        'roi_voxels', 'mean_afd', 'tract_name', 'cleaned'
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_metrics)

    print(f"\nCSV written: {output_path} ({len(all_metrics)} rows)")


# ============================================================
# HTML report
# ============================================================

def generate_html_report(all_metrics, output_path, cfg, use_shrink=False, use_cleaned=False,
                          cerebellum_via_thalamus=False):
    """Generate HTML report with tables and charts."""
    if not all_metrics:
        print("No metrics for report!")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Group data
    groups = defaultdict(list)
    for m in all_metrics:
        groups[m['group']].append(m)

    # Determine all unique segment-target-hemi combos
    all_combos = set()
    for m in all_metrics:
        all_combos.add((m['hemisphere'], m['segment'], m['target']))

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MBRAIN Per-Tract Metrics Report (3-Segment Comparison)</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
<style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
           margin: 0; padding: 20px; background: #f5f5f5; color: #333; }}
    .container {{ max-width: 1400px; margin: 0 auto; }}
    h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
    h2 {{ color: #34495e; margin-top: 30px; }}
    h3 {{ color: #7f8c8d; }}
    .summary-box {{ background: white; border-radius: 8px; padding: 20px;
                    margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
    .mapping-box {{ background: #fdf6e3; border-left: 4px solid #f39c12; padding: 12px 18px;
                    margin: 10px 0; font-size: 13px; }}
    .mapping-box code {{ background: #eee; padding: 1px 5px; border-radius: 3px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }}
    th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: right; }}
    th {{ background: #3498db; color: white; position: sticky; top: 0; }}
    td:first-child, td:nth-child(2) {{ text-align: left; }}
    tr:nth-child(even) {{ background: #f8f9fa; }}
    tr:hover {{ background: #e8f4fd; }}
    .chart-container {{ background: white; border-radius: 8px; padding: 20px;
                       margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
    .chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    .group-ctrl {{ color: #2ecc71; }} .group-ap {{ color: #e74c3c; }} .group-aptr {{ color: #f39c12; }}
    .metric-high {{ background: #d4edda; }} .metric-low {{ background: #f8d7da; }}
    .zero {{ color: #ccc; }}
    .tabs {{ display: flex; gap: 5px; margin: 20px 0 10px; }}
    .tab {{ padding: 8px 16px; background: #ecf0f1; border: none; border-radius: 5px 5px 0 0;
            cursor: pointer; font-size: 14px; }}
    .tab.active {{ background: #3498db; color: white; }}
    .tab-content {{ display: none; }} .tab-content.active {{ display: block; }}
    .footer {{ text-align: center; color: #95a5a6; margin-top: 40px; padding: 20px; font-size: 12px; }}
</style>
</head>
<body>
<div class="container">
<h1>MBRAIN Per-Tract Metrics Report</h1>
<p>Generated: {timestamp}</p>
{'<div class="mapping-box" style="border-left-color: #e74c3c;"><strong>&#x26A0; SHRINK MODE:</strong> AP/APTR passthrough tracts extracted using <b>eroded (shrink) segment ROIs</b> from <code>graft_segments_fod/shrink/</code>. Density normalized by shrink ROI voxel count.</div>' if use_shrink else ''}
{'<div class="mapping-box" style="border-left-color: #27ae60;"><strong>&#x2705; CLEANED TRACTS:</strong> Tracts processed with MATLAB AFQ_removeFiberOutliers (maxDist=3, maxLen=3, numNodes=100). Outlier streamlines removed before computing metrics. Falls back to raw .tck if _cleaned.tck not available or &lt;20 streamlines.</div>' if use_cleaned else ''}
{'<div class="mapping-box" style="border-left-color: #8e44ad;"><strong>&#x1F9E0; CEREBELLUM VIA THALAMUS:</strong> For the <b>cerebellum</b> target, tracts are additionally filtered to only keep streamlines that also pass through the <b>thalamus</b> ROI (tckedit -include). Other targets (striatum, olfactory_cortex, thalamus) are unchanged.</div>' if cerebellum_via_thalamus else ''}
<div class="mapping-box">
<strong>Segment alignment (3-segment comparison):</strong><br>
CTRL uses <code>isocortex_segments_v2/</code> with 5 AP segments &mdash; only 3 kept for comparison:<br>
<code>anterior</code> &harr; AP/APTR <code>seg9_graft_*_anterior_dorsal</code>
&emsp;|&emsp; <code>central</code> &harr; <code>seg9_graft_*_central_dorsal</code>
&emsp;|&emsp; <code>posterior</code> &harr; <code>seg9_graft_*_posterior_dorsal</code><br>
<strong>Target alignment:</strong> CTRL <code>olfactory_cortex</code> = AP/APTR <code>Palecortex</code> (piriform complex)
</div>
<div class="mapping-box" style="border-left-color: #3498db;">
<strong>Metric definitions:</strong><br>
<b>Count</b>: Raw number of streamlines in each passthrough tract. Least meaningful &mdash; driven by tractography parameters, not biology.<br>
<b>SIFT2</b>: Sum of SIFT2 weights for matched streamlines. Biologically-corrected total connection strength (SIFT2 adjusts weights to fit the FOD). Useful for absolute connectivity, but confounded by segment ROI size.<br>
<b>Density</b>: SIFT2 sum normalized by seed ROI voxel count (&times;1000). Most useful for cross-group connectivity comparison &mdash; controls for differences in segment ROI size.<br>
<b>AFD</b>: Mean Apparent Fiber Density along each tract (fixel2voxel &rarr; tcksample). Measures microstructural fiber integrity along the tract path, independent of streamline count.
</div>
"""

    # ---- Summary section ----
    html += '<div class="summary-box">\n<h2>Summary</h2>\n'
    for gname in ['CTRL', 'AP', 'APTR']:
        gdata = groups.get(gname, [])
        if not gdata:
            continue
        subjects = sorted(set(m['subject'] for m in gdata))
        total_streams = sum(m['streamline_count'] for m in gdata)
        total_sift2 = sum(m['sift2_sum'] for m in gdata)
        afd_vals = [m['mean_afd'] for m in gdata if m['mean_afd'] is not None]
        mean_afd = np.mean(afd_vals) if afd_vals else 0
        html += f'<h3 class="group-{gname.lower()}">{gname} ({len(subjects)} subjects)</h3>\n'
        html += f'<p>Subjects: {", ".join(subjects)}<br>'
        html += f'Total streamlines: {total_streams:,} | Total SIFT2: {total_sift2:.2f} | Mean AFD: {mean_afd:.4f}</p>\n'
    html += '</div>\n'

    # ---- Tabs: Left / Right ----
    html += '<div class="tabs">\n'
    html += '<button class="tab active" onclick="showTab(\'left\')">Left Hemisphere</button>\n'
    html += '<button class="tab" onclick="showTab(\'right\')">Right Hemisphere</button>\n'
    html += '</div>\n'

    for hemi in ['left', 'right']:
        active = ' active' if hemi == 'left' else ''
        html += f'<div id="tab-{hemi}" class="tab-content{active}">\n'
        html += f'<h2>{hemi.capitalize()} Hemisphere</h2>\n'

        # Get unique targets for this hemisphere
        hemi_targets = sorted(set(
            m['target'] for m in all_metrics if m['hemisphere'] == hemi
        ))

        for target in hemi_targets:
            html += f'<div class="summary-box">\n<h3>Target: {target.capitalize()}</h3>\n'

            # Table
            html += '<table>\n<thead><tr>'
            html += '<th>Subject</th><th>Group</th>'
            # Get all segments for this target
            seg_set = sorted(set(
                m['segment'] for m in all_metrics
                if m['hemisphere'] == hemi and m['target'] == target
            ))
            for seg in seg_set:
                html += f'<th colspan="4">{seg}</th>'
            html += '</tr>\n<tr><th></th><th></th>'
            for _ in seg_set:
                html += '<th>Count</th><th>SIFT2</th><th>Density</th><th>AFD</th>'
            html += '</tr></thead>\n<tbody>\n'

            # Rows per subject
            subject_list = sorted(set(
                (m['subject'], m['group']) for m in all_metrics
                if m['hemisphere'] == hemi and m['target'] == target
            ), key=lambda x: (x[1], x[0]))

            for subj_id, grp in subject_list:
                gclass = f'group-{grp.lower()}'
                html += f'<tr><td>{subj_id}</td><td class="{gclass}">{grp}</td>'

                for seg in seg_set:
                    match = [m for m in all_metrics
                             if m['subject'] == subj_id and m['hemisphere'] == hemi
                             and m['target'] == target and m['segment'] == seg]
                    if match:
                        m = match[0]
                        cnt = m['streamline_count']
                        s2 = m['sift2_sum']
                        dens = m['sift2_density']
                        afd = m['mean_afd']
                        cnt_cls = ' class="zero"' if cnt == 0 else ''
                        html += f'<td{cnt_cls}>{cnt:,}</td>'
                        html += f'<td>{s2:.2f}</td>'
                        html += f'<td>{dens:.2f}</td>'
                        html += f'<td>{afd:.4f}</td>' if afd is not None else '<td class="zero">-</td>'
                    else:
                        html += '<td>-</td>' * 4
                html += '</tr>\n'

            # Group means row
            for gname in ['CTRL', 'AP', 'APTR']:
                g_subjects = [s for s, g in subject_list if g == gname]
                if not g_subjects:
                    continue
                gclass = f'group-{gname.lower()}'
                html += f'<tr style="font-weight:bold;background:#eee"><td>{gname} Mean</td><td class="{gclass}">{gname}</td>'
                for seg in seg_set:
                    g_matches = [m for m in all_metrics
                                 if m['group'] == gname and m['hemisphere'] == hemi
                                 and m['target'] == target and m['segment'] == seg]
                    if g_matches:
                        cnt_mean = np.mean([m['streamline_count'] for m in g_matches])
                        s2_mean = np.mean([m['sift2_sum'] for m in g_matches])
                        dens_mean = np.mean([m['sift2_density'] for m in g_matches])
                        afd_vals = [m['mean_afd'] for m in g_matches if m['mean_afd'] is not None]
                        afd_mean = np.mean(afd_vals) if afd_vals else None
                        html += f'<td>{cnt_mean:.0f}</td><td>{s2_mean:.2f}</td><td>{dens_mean:.2f}</td>'
                        html += f'<td>{afd_mean:.4f}</td>' if afd_mean else '<td>-</td>'
                    else:
                        html += '<td>-</td>' * 4
                html += '</tr>\n'

            html += '</tbody></table>\n'

            # Charts - 4 dot+mean plots (Count, SIFT2, Density, AFD)
            chart_id = f'chart_{hemi}_{target}'.replace(' ', '_')
            metrics_defs = [
                ('count', 'streamline_count',
                 f'Streamline Count \u2192 {target.capitalize()} ({hemi})', 'Count'),
                ('sift2', 'sift2_sum',
                 f'SIFT2 Sum \u2192 {target.capitalize()} ({hemi})', 'SIFT2 Sum'),
                ('density', 'sift2_density',
                 f'SIFT2 Density \u2192 {target.capitalize()} ({hemi})', 'Density (\u00d71000)'),
                ('afd', 'mean_afd',
                 f'Mean AFD \u2192 {target.capitalize()} ({hemi})', 'Mean AFD'),
            ]

            html += f'<div class="chart-grid">\n'
            for msuffix, _, _, _ in metrics_defs:
                html += f'<div class="chart-container"><canvas id="{chart_id}_{msuffix}"></canvas></div>\n'
            html += f'</div>\n'

            html += f'<script>\n'

            group_colors = {'CTRL': '#2ecc71', 'AP': '#e74c3c', 'APTR': '#f39c12'}
            group_colors_alpha = {'CTRL': 'rgba(46,204,113,0.4)', 'AP': 'rgba(231,76,60,0.4)',
                                  'APTR': 'rgba(243,156,18,0.4)'}
            group_offsets = {'CTRL': -0.2, 'AP': 0.0, 'APTR': 0.2}

            for msuffix, mkey, mtitle, ylabel in metrics_defs:
                datasets = []
                for gname in ['CTRL', 'AP', 'APTR']:
                    g_data = groups.get(gname, [])
                    if not g_data:
                        continue
                    g_offset = group_offsets[gname]
                    g_subjects = sorted(set(m['subject'] for m in g_data))
                    n_subj = len(g_subjects)
                    sub_jitter = list(np.linspace(-0.06, 0.06, n_subj)) if n_subj > 1 else [0.0]

                    points = []
                    mean_points = []
                    for seg_idx, seg in enumerate(seg_set):
                        seg_vals = []
                        for si, subj in enumerate(g_subjects):
                            matches = [m for m in g_data
                                       if m['subject'] == subj and m['hemisphere'] == hemi
                                       and m['target'] == target and m['segment'] == seg]
                            if matches:
                                val = matches[0].get(mkey)
                                if val is not None:
                                    points.append({
                                        'x': round(seg_idx + g_offset + sub_jitter[si], 3),
                                        'y': round(float(val), 4)
                                    })
                                    seg_vals.append(float(val))
                        if seg_vals:
                            mean_points.append({
                                'x': round(seg_idx + g_offset, 2),
                                'y': round(float(np.mean(seg_vals)), 4)
                            })

                    # Individual dots
                    datasets.append({
                        'label': gname,
                        'data': points,
                        'backgroundColor': group_colors_alpha[gname],
                        'borderColor': group_colors[gname],
                        'borderWidth': 1,
                        'pointRadius': 5,
                        'pointStyle': 'circle',
                        'showLine': False,
                    })
                    # Mean markers (diamond)
                    datasets.append({
                        'label': f'{gname} mean',
                        'data': mean_points,
                        'backgroundColor': group_colors[gname],
                        'borderColor': '#333',
                        'borderWidth': 2,
                        'pointRadius': 9,
                        'pointStyle': 'rectRot',
                        'showLine': False,
                    })

                seg_labels_json = json.dumps(list(seg_set))
                n_segs = len(seg_set)

                html += f"""
new Chart(document.getElementById('{chart_id}_{msuffix}'), {{
    type: 'scatter',
    data: {{ datasets: {json.dumps(datasets)} }},
    options: {{
        responsive: true,
        plugins: {{
            title: {{ display: true, text: '{mtitle}' }},
            legend: {{
                labels: {{
                    filter: function(item) {{ return !item.text.includes('mean'); }}
                }}
            }},
            tooltip: {{
                callbacks: {{
                    label: function(ctx) {{
                        var lbl = ctx.dataset.label || '';
                        return lbl + ': ' + ctx.parsed.y.toFixed(4);
                    }}
                }}
            }}
        }},
        scales: {{
            x: {{
                type: 'linear',
                min: -0.5,
                max: {n_segs - 0.5},
                ticks: {{
                    stepSize: 1,
                    callback: function(value) {{
                        var labels = {seg_labels_json};
                        return labels[value] || '';
                    }}
                }},
                title: {{ display: false }}
            }},
            y: {{
                beginAtZero: true,
                title: {{ display: true, text: '{ylabel}' }}
            }}
        }}
    }}
}});
"""

            html += '</script>\n'
            html += '</div>\n'  # summary-box

        html += '</div>\n'  # tab-content

    # Tab switching JS
    html += """
<script>
function showTab(hemi) {
    document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
    document.getElementById('tab-' + hemi).classList.add('active');
    event.target.classList.add('active');
}
</script>
"""

    # Footer
    html += f"""
<div class="footer">
    <p>MBRAIN Per-Tract Metrics Report | Generated {timestamp}</p>
    <p>Pipeline: mrtrix_preproc_flip1 | SIFT2 + Fixel AFD Analysis</p>
</div>
</div>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)
    print(f"HTML report written: {output_path}")


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate per-tract metrics report (SIFT2 connectivity + AFD)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--group', choices=['CTRL', 'AP', 'APTR'],
                        help='Process specific group only')
    parser.add_argument('--subject', type=str,
                        help='Process specific subject only')
    parser.add_argument('--skip-extraction', action='store_true',
                        help='Skip SIFT2 weight extraction (reuse existing files)')
    parser.add_argument('--shrink', action='store_true',
                        help='Use shrunk segment ROIs for AP/APTR (passthrough_tracts_shrink/)')
    parser.add_argument('--cleaned', action='store_true',
                        help='Use MATLAB-cleaned tracts (_cleaned.tck) when available')
    parser.add_argument('--cerebellum-via-thalamus', action='store_true',
                        help='For cerebellum target, only keep tracts passing through thalamus')
    parser.add_argument('--output-dir', type=str, default='reports',
                        help='Output directory for report files')
    parser.add_argument('--version', type=str, default='v01',
                        help='Report version suffix')
    parser.add_argument('--derivatives-dir', type=str,
                        help='Override derivatives directory (default: from config1.yaml)')

    args = parser.parse_args()

    # Override derivatives directory if specified
    if args.derivatives_dir:
        os.environ['MBRAIN_DERIVATIVES_DIR'] = args.derivatives_dir

    print("=" * 70)
    print("MBRAIN PER-TRACT METRICS REPORT")
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
    for i, (subject_id, session, group) in enumerate(subjects, 1):
        print(f"\n{'='*60}")
        print(f"[{i}/{len(subjects)}] {subject_id} ({group})")
        print(f"{'='*60}")

        metrics = compute_subject_metrics(
            subject_id, session, group, deriv_dir,
            skip_extraction=args.skip_extraction,
            use_shrink=args.shrink,
            use_cleaned=args.cleaned,
            cerebellum_via_thalamus=args.cerebellum_via_thalamus
        )
        all_metrics.extend(metrics)
        print(f"  → {len(metrics)} tract metrics computed")

    # Output
    out_dir = Path(args.output_dir)
    csv_path = out_dir / f"tract_metrics_{args.version}.csv"
    html_path = out_dir / f"tract_metrics_report_{args.version}.html"

    write_csv(all_metrics, csv_path)
    generate_html_report(all_metrics, html_path, cfg, use_shrink=args.shrink, use_cleaned=args.cleaned,
                          cerebellum_via_thalamus=args.cerebellum_via_thalamus)

    print(f"\n{'='*70}")
    print(f"REPORT COMPLETE")
    print(f"  CSV:  {csv_path}")
    print(f"  HTML: {html_path}")
    print(f"  Total metrics: {len(all_metrics)}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
