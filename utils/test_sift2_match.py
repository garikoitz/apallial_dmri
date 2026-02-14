#!/usr/bin/env python3
"""
Test SIFT2 hash matching for different tract file types.

Verifies that passthrough tracts (raw, shrunk, cleaned, cleaned-shrunk) correctly
match against the wholebrain tractogram hash index for SIFT2 weight extraction.

This is a diagnostic tool to confirm that:
  - Raw tracts (from tckedit) hash-match the wholebrain (expected: yes)
  - Shrunk tracts (from tckedit with eroded ROIs) hash-match (expected: yes)
  - MATLAB-cleaned tracts hash-match (expected: NO — AFQ resamples vertices)
  - The count_tck_streamlines fallback correctly counts cleaned streamlines

Usage:
    python utils/test_sift2_match.py                          # Default subject
    python utils/test_sift2_match.py --subject sub-129-921-R1F --group CTRL

Output:
    Prints match counts and SIFT2 sums for each tract variant.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'functions'))

import numpy as np
from generate_tract_metrics_report import (
    parse_tck_header, build_wholebrain_hash_index,
    match_passthrough_weights, load_sift2_weights, count_tck_streamlines
)
from config_loader import load_config, get_subject_groups


def get_default_subject():
    """Return (subject_id, session, group) for the first AP subject."""
    cfg = load_config()
    groups = cfg['groups']
    for grp_name in ['AP', 'CTRL', 'APTR']:
        if grp_name in groups and groups[grp_name]:
            s = groups[grp_name][0]
            return s['id'], s['session'], grp_name
    return None, None, None


def run_test(subject_id, session, group, deriv_dir):
    """Test SIFT2 matching for all tract variants of one subject."""
    base = os.path.join(deriv_dir, subject_id, session, 'dwi', 'wholebrain')
    wb = os.path.join(base, 'wholebrain_1000k.tck')
    sift2_file = os.path.join(base, 'sift2_weights.txt')

    if not os.path.exists(wb):
        print(f"ERROR: Wholebrain not found: {wb}")
        return
    if not os.path.exists(sift2_file):
        print(f"ERROR: SIFT2 weights not found: {sift2_file}")
        return

    print(f"Subject: {subject_id} ({group})")
    print(f"Loading wholebrain hash index...")
    w = load_sift2_weights(sift2_file)
    hm = build_wholebrain_hash_index(wb)

    # Pick a tract name based on group
    if group == 'CTRL':
        tract = 'isocortex_left_anterior_to_thalamus_left'
    else:
        tract = 'seg9_graft_left_anterior_dorsal_to_Thalamus_left'

    tests = [
        ('Raw non-shrink',    os.path.join(base, 'passthrough_tracts', f'{tract}.tck')),
        ('Cleaned non-shrink', os.path.join(base, 'passthrough_tracts', f'{tract}_cleaned.tck')),
    ]
    if group in ('AP', 'APTR'):
        tests.extend([
            ('Raw shrink',        os.path.join(base, 'passthrough_tracts_shrink', f'{tract}.tck')),
            ('Cleaned shrink',    os.path.join(base, 'passthrough_tracts_shrink', f'{tract}_cleaned.tck')),
        ])

    print(f"\nTract: {tract}")
    print(f"{'Label':25s}  {'matched':>8s}  {'tck_count':>9s}  {'sift2_sum':>10s}")
    print("-" * 60)
    for label, path in tests:
        if os.path.exists(path):
            n_matched, s, m, _ = match_passthrough_weights(path, hm, w)
            n_count = count_tck_streamlines(path)
            print(f"{label:25s}  {n_matched:8d}  {n_count:9d}  {s:10.4f}")
        else:
            print(f"{label:25s}  {'FILE NOT FOUND':>30s}")


def main():
    parser = argparse.ArgumentParser(
        description='Test SIFT2 hash matching for raw/shrink/cleaned tracts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--subject', type=str, help='Subject ID')
    parser.add_argument('--session', type=str, help='Session ID')
    parser.add_argument('--group', type=str, choices=['CTRL', 'AP', 'APTR'],
                        help='Group name')
    args = parser.parse_args()

    cfg = load_config()
    deriv_dir = cfg['paths']['derivatives_dir']

    if args.subject:
        # Find session from config if not provided
        if not args.session:
            groups = cfg['groups']
            for grp_name, subjs in groups.items():
                for s in subjs:
                    if s['id'] == args.subject:
                        args.session = s['session']
                        if not args.group:
                            args.group = grp_name
                        break
        if not args.session:
            print(f"ERROR: Subject {args.subject} not found in config")
            sys.exit(1)
        run_test(args.subject, args.session, args.group, deriv_dir)
    else:
        # Test one subject per group
        groups = cfg['groups']
        for grp_name in ['CTRL', 'AP', 'APTR']:
            if grp_name in groups and groups[grp_name]:
                s = groups[grp_name][0]
                print(f"\n{'='*60}")
                run_test(s['id'], s['session'], grp_name, deriv_dir)


if __name__ == '__main__':
    main()
