#!/usr/bin/env python3
"""
Open mrview for all 4 APTR mice with T2w_highres_flip base image
and colored fixel cluster overlays.

Colors:
  - Green shades  → DV clusters
  - Red shades    → AP clusters
  - Blue shades   → ML clusters

Usage:
    python tmp/open_fixel_clustering_APTR_mice.py
    python tmp/open_fixel_clustering_APTR_mice.py --opacity 0.7
    python tmp/open_fixel_clustering_APTR_mice.py --subjects sub-AP72-862-901-NmF sub-Apallial17NmF
"""

import argparse
import json
import os
from pathlib import Path

DERIV = Path('/home/glerma/glerma/PROYECTOS/mbrain/data/derivatives/mrtrix_newmask')

SUBJECTS = {
    'sub-AP72-862-901-NmF': 'ses-251210093342',
    'sub-AP72-862-901-R1M': 'ses-251210100043',
    'sub-Apallial17NmF':    'ses-240816091633',
    'sub-Apallial20R1F':    'ses-240909094232',
}

# Color palettes (R, G, B) in 0–255
GREEN_SHADES = [
    (0, 255, 0),      # bright green
    (0, 191, 0),      # medium green
    (51, 230, 51),    # light green
    (0, 140, 0),      # dark green
    (102, 255, 102),  # pale green
    (0, 217, 77),     # emerald
]

RED_SHADES = [
    (255, 0, 0),      # bright red
    (217, 0, 0),      # dark red
    (255, 77, 77),    # salmon
    (179, 0, 0),      # maroon
    (255, 38, 38),    # medium red
    (230, 51, 0),     # red-orange
]

BLUE_SHADES = [
    (0, 0, 255),      # bright blue
    (0, 0, 191),      # medium blue
    (77, 77, 255),    # light blue
    (0, 0, 140),      # dark blue
    (51, 128, 255),   # sky blue
    (0, 77, 230),     # royal blue
]


def get_rgb(direction: str, shade_idx: int) -> tuple:
    """Return (R, G, B) for a given direction and shade index."""
    if direction == 'DV':
        return GREEN_SHADES[shade_idx % len(GREEN_SHADES)]
    elif direction == 'AP':
        return RED_SHADES[shade_idx % len(RED_SHADES)]
    elif direction == 'ML':
        return BLUE_SHADES[shade_idx % len(BLUE_SHADES)]
    return (255, 255, 255)


def open_subject(sub: str, ses: str, opacity: float = 0.5) -> None:
    """Launch mrview for one subject with T2w base + colored cluster overlays."""
    dwi_dir = DERIV / sub / ses / 'dwi'
    clust_dir = dwi_dir / 'fixel_clusters'
    t2w = dwi_dir / 't2w_highres_flip.mif'
    info_path = clust_dir / 'cluster_info.json'

    if not t2w.exists():
        print(f'  WARNING: {t2w} not found, skipping {sub}')
        return
    if not info_path.exists():
        print(f'  WARNING: {info_path} not found, skipping {sub}')
        return

    with open(info_path) as f:
        meta = json.load(f)

    cmd_parts = [f'mrview "{t2w}"']
    shade_counter = {'DV': 0, 'AP': 0, 'ML': 0}

    for ci in meta['clusters']:
        cid = ci['cluster_id']
        d = ci['main_direction']
        r, g, b = get_rgb(d, shade_counter[d])
        shade_counter[d] += 1
        overlay = clust_dir / f'cluster_{cid}_{d}_colored.nii.gz'
        if not overlay.exists():
            # Fall back to .mif mask
            overlay = clust_dir / f'cluster_{cid}_{d}.mif'
        cmd_parts.append(
            f'-overlay.load "{overlay}" '
            f'-overlay.opacity {opacity} '
            f'-overlay.colour {r},{g},{b} '
            f'-overlay.threshold_min 0.5'
        )

    cmd = ' '.join(cmd_parts) + ' &'
    print(f'Opening {sub} …')
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser(
        description='Open mrview with T2w base + colored fixel cluster overlays for APTR mice')
    parser.add_argument('--subjects', nargs='+', default=None,
                        help='Subset of subject IDs to open (default: all 4)')
    parser.add_argument('--opacity', type=float, default=0.5,
                        help='Overlay opacity (default: 0.5)')
    args = parser.parse_args()

    if args.subjects:
        subs = {s: SUBJECTS[s] for s in args.subjects if s in SUBJECTS}
    else:
        subs = SUBJECTS

    for sub, ses in subs.items():
        open_subject(sub, ses, opacity=args.opacity)

    print(f'\nDone — {len(subs)} mrview window(s) launched.')


if __name__ == '__main__':
    main()
