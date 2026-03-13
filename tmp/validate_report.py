#!/usr/bin/env python3
"""Validate tract_metrics_report HTML against CSVs."""
import csv, re, ast, numpy as np, math, json

REPORTS = '/home/glerma/glerma/soft/mbrain/reports'
MAIN_CSV = f'{REPORTS}/tract_metrics_v20260306_newmask.csv'
HTML_FILE = f'{REPORTS}/tract_metrics_report_v20260306_newmask.html'

# Read main CSV
with open(MAIN_CSV) as f:
    data = list(csv.DictReader(f))
print(f"Main CSV: {len(data)} rows")

# Read HTML
with open(HTML_FILE) as f:
    html = f.read()

# ── Extract chart titles and dataset arrays ──
# Each chart block starts with (function(){...})();
# We need to extract the datasets JSON, but keys are unquoted JS.
# Strategy: extract the raw datasets string and fix it for JSON parsing.

# Find all chart instantiations
chart_re = re.compile(
    r"text:\s*'([^']+)'.*?datasets:\s*(\[.*?\])\s*\}\s*,\s*\n\s*options:",
    re.DOTALL
)
matches = chart_re.findall(html)
print(f"Charts found: {len(matches)}")

# Parse JS object notation into Python dict
def js_to_json(s):
    """Convert JS object notation to valid JSON by quoting unquoted keys."""
    # Quote unquoted keys: word followed by :
    s = re.sub(r'(?<=[{,\s])(\w+)\s*:', r'"\1":', s)
    # Remove trailing commas before } or ]
    s = re.sub(r',\s*([}\]])', r'\1', s)
    return s

groups = ['CTRL', 'AP', 'APTR']
targets_order = ['striatum', 'thalamus', 'paleocortex', 'cerebellum']
segments_order = ['anterior', 'central', 'posterior']
metrics_order = [
    ('streamline_count', 'Streamline Count'),
    ('sift2_weighted', 'SIFT2 Weighted'),
    ('density', 'Density'),
]

errors = []
chart_idx = 0

for hemi in ['left', 'right']:
    for seg in segments_order:
        for metric_key, metric_label in metrics_order:
            if chart_idx >= len(matches):
                errors.append(f"Missing chart {chart_idx}")
                chart_idx += 1
                continue

            title, ds_raw = matches[chart_idx]
            # Verify title matches expected
            expected_title = f"{metric_label}: {seg.capitalize()} segment ({hemi.capitalize()})"
            if title != expected_title:
                errors.append(f"Chart {chart_idx} title mismatch: '{title}' vs expected '{expected_title}'")

            # Parse datasets
            ds_json = js_to_json(ds_raw)
            try:
                datasets = json.loads(ds_json)
            except json.JSONDecodeError as e:
                errors.append(f"Chart {chart_idx} JSON parse error: {e}")
                chart_idx += 1
                continue

            for gi, gn in enumerate(groups):
                if gi >= len(datasets):
                    errors.append(f"Chart {chart_idx}: missing dataset for {gn}")
                    continue
                ds = datasets[gi]
                chart_means = ds['data']
                chart_points = ds['_dataPoints']
                chart_errbars = ds['_errorBars']

                # Verify label
                if ds['label'] != gn:
                    errors.append(f"Chart {chart_idx} ds[{gi}] label '{ds['label']}' != '{gn}'")

                for ti, tgt in enumerate(targets_order):
                    # Get CSV values for this group/hemi/seg/target
                    csv_vals = sorted([
                        float(r[metric_key])
                        for r in data
                        if r['group'] == gn and r['hemisphere'] == hemi
                        and r['segment'] == seg and r['target'] == tgt
                    ])

                    # Check mean
                    csv_mean = float(np.mean(csv_vals)) if csv_vals else 0.0
                    chart_mean = chart_means[ti]
                    if abs(csv_mean - chart_mean) > 0.01:
                        errors.append(
                            f"MEAN chart_{chart_idx} {gn} {hemi}/{seg}/{tgt}: "
                            f"CSV={csv_mean:.6f} Chart={chart_mean}"
                        )

                    # Check individual points
                    chart_pts_sorted = sorted(chart_points[ti])
                    if len(csv_vals) != len(chart_pts_sorted):
                        errors.append(
                            f"NPTS chart_{chart_idx} {gn} {hemi}/{seg}/{tgt}: "
                            f"CSV={len(csv_vals)} Chart={len(chart_pts_sorted)}"
                        )
                    else:
                        for ci, (cv, chv) in enumerate(zip(csv_vals, chart_pts_sorted)):
                            if abs(cv - chv) > 0.01:
                                errors.append(
                                    f"PT chart_{chart_idx} {gn} {hemi}/{seg}/{tgt}[{ci}]: "
                                    f"CSV={cv} Chart={chv}"
                                )

                    # Check SEM
                    if len(csv_vals) > 1:
                        csv_sem = float(np.std(csv_vals, ddof=1) / math.sqrt(len(csv_vals)))
                    else:
                        csv_sem = 0.0
                    chart_sem = chart_errbars[ti]
                    if abs(csv_sem - chart_sem) > 0.01:
                        errors.append(
                            f"SEM chart_{chart_idx} {gn} {hemi}/{seg}/{tgt}: "
                            f"CSV={csv_sem:.6f} Chart={chart_sem}"
                        )

            chart_idx += 1

print()
if errors:
    print(f"❌ Found {len(errors)} issues:")
    for e in errors:
        print(f"  {e}")
else:
    print("✅ ALL 18 charts verified — means, data points, SEMs all match CSV.")

# ── Wide-format CSV check ──
print()
wide_errors = []
for metric_key, suffix in [('streamline_count', 'count'), ('sift2_weighted', 'sift2'), ('density', 'density')]:
    with open(f'{REPORTS}/tract_metrics_{suffix}_v20260306_newmask.csv') as f:
        wreader = csv.DictReader(f)
        wdata = list(wreader)
        wcols = wreader.fieldnames
    for wrow in wdata:
        subj = wrow['subject']
        for col in wcols[2:]:
            parts = col.split('_')
            hemi, seg, tgt = parts[0], parts[1], parts[2]
            wval = float(wrow[col]) if wrow[col] else 0
            csv_match = [r for r in data if r['subject'] == subj
                         and r['hemisphere'] == hemi and r['segment'] == seg
                         and r['target'] == tgt]
            if csv_match:
                mval = float(csv_match[0][metric_key])
                if abs(wval - mval) > 0.001:
                    wide_errors.append(f"{suffix} {subj} {col}: wide={wval} main={mval}")

if wide_errors:
    print(f"❌ {len(wide_errors)} wide-format mismatches:")
    for e in wide_errors[:10]:
        print(f"  {e}")
else:
    print("✅ All wide-format CSVs match main CSV.")

# ── SIFT2/Count ratio check ──
print()
print("=== SIFT2 / Count ratio anomalies (>30× or <0.01× per streamline) ===")
anom = 0
for row in data:
    cnt = int(float(row['streamline_count']))
    s2 = float(row['sift2_weighted'])
    if cnt > 0 and s2 > 0:
        ratio = s2 / cnt
        if ratio > 30 or ratio < 0.01:
            anom += 1
            print(f"  ⚠  {row['subject']} {row['hemisphere']}/{row['segment']}/{row['target']}: "
                  f"count={cnt} sift2={s2:.1f} ratio={ratio:.1f}")
if anom == 0:
    print("  None found.")

# ── Check SIFT2 computed from raw vs count from cleaned ──
# The generation script uses n_streamlines from cleaned tracts but SIFT2 from raw tracts.
# Flag cases where the raw tract SIFT2 is very different from what we'd expect.
print()
print("=== Count(cleaned) vs SIFT2(raw) per-group summary ===")
for gn in groups:
    g_rows = [r for r in data if r['group'] == gn]
    cnts = [float(r['streamline_count']) for r in g_rows if float(r['streamline_count']) > 0]
    s2s = [float(r['sift2_weighted']) for r in g_rows if float(r['sift2_weighted']) > 0]
    if cnts and s2s:
        ratios = []
        for r in g_rows:
            c, s = float(r['streamline_count']), float(r['sift2_weighted'])
            if c > 0 and s > 0:
                ratios.append(s / c)
        if ratios:
            print(f"  {gn}: median SIFT2/count ratio = {np.median(ratios):.2f}, "
                  f"range [{min(ratios):.2f}, {max(ratios):.2f}], n={len(ratios)} tracts")
