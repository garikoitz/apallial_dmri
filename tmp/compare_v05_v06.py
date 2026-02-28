#!/usr/bin/env python3
"""Compare tract metrics v05 vs v06."""
import csv

def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))

v05 = load_csv('reports/tract_metrics_v05_cleaned.csv')
v06 = load_csv('reports/tract_metrics_v06_cleaned.csv')

def make_key(r):
    return (r['subject'], r['hemisphere'], r['segment'], r['target'])

v05_map = {make_key(r): r for r in v05}
v06_map = {make_key(r): r for r in v06}

# Compare affected subjects
affected = ['sub-769-C-O5-R1M', 'sub-Apallial17NmF', 'sub-Apallial20R1F']

print("=" * 130)
print("CHANGES IN AFFECTED SUBJECTS (v05 → v06)")
print("=" * 130)
print(f"{'Subject':>20} | {'Hemi':>5} | {'Segment':>10} | {'Target':>17} | {'v05_Cnt':>7} | {'v06_Cnt':>7} | {'v05_SIFT2':>10} | {'v06_SIFT2':>10} | {'v05_AFD':>8} | {'v06_AFD':>8}")
print("-" * 130)

n_changed = 0
for key in sorted(v06_map.keys()):
    subj = key[0]
    if subj not in affected:
        continue
    r5 = v05_map.get(key, {})
    r6 = v06_map[key]
    
    c5 = int(r5.get('streamline_count', 0))
    c6 = int(r6.get('streamline_count', 0))
    s5 = float(r5.get('sift2_sum', 0))
    s6 = float(r6.get('sift2_sum', 0))
    a5 = r5.get('mean_afd', '')
    a6 = r6.get('mean_afd', '')
    
    if c5 != c6 or abs(s5 - s6) > 0.01:
        n_changed += 1
        short = subj.replace('sub-', '')[:18]
        a5_str = f"{float(a5):.4f}" if a5 else "-"
        a6_str = f"{float(a6):.4f}" if a6 else "-"
        print(f"{short:>20} | {key[1]:>5} | {key[2]:>10} | {key[3]:>17} | {c5:>7} | {c6:>7} | {s5:>10.2f} | {s6:>10.2f} | {a5_str:>8} | {a6_str:>8}")

print(f"\nTotal changed rows: {n_changed}")

# Group summaries
print("\n" + "=" * 80)
print("GROUP SUMMARIES")
print("=" * 80)

for version, data, label in [('v05', v05, 'v05_cleaned'), ('v06', v06, 'v06_cleaned')]:
    print(f"\n--- {label} ---")
    for grp in ['CTRL', 'AP', 'APTR']:
        rows = [r for r in data if r['group'] == grp]
        total_count = sum(int(r['streamline_count']) for r in rows)
        total_sift2 = sum(float(r['sift2_sum']) for r in rows)
        afds = [float(r['mean_afd']) for r in rows if r['mean_afd']]
        mean_afd = sum(afds) / len(afds) if afds else 0
        print(f"  {grp:>5}: count={total_count:>8,} | SIFT2={total_sift2:>12,.2f} | mean_AFD={mean_afd:.4f} ({len(afds)} tracts with AFD)")

# Unaffected subjects: verify identical
print("\n" + "=" * 80)
print("UNAFFECTED SUBJECTS CHECK")
print("=" * 80)
n_diff = 0
for key in v06_map:
    if key[0] in affected:
        continue
    r5 = v05_map.get(key, {})
    r6 = v06_map[key]
    c5 = int(r5.get('streamline_count', 0))
    c6 = int(r6.get('streamline_count', 0))
    s5 = float(r5.get('sift2_sum', 0))
    s6 = float(r6.get('sift2_sum', 0))
    if c5 != c6 or abs(s5 - s6) > 0.01:
        n_diff += 1
        print(f"  DIFF: {key}")

if n_diff == 0:
    print("  ✓ All unaffected subjects identical between v05 and v06")
else:
    print(f"  ⚠ {n_diff} unexpected differences!")
