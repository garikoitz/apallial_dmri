#!/usr/bin/env python3
"""Test fast .tck parsing and SIFT2 weight extraction approach."""
import numpy as np
import time
import os

SUBJ_DIR = "/Users/glerma/soft/mbrain/data/derivatives/mrtrix_preproc_flip1/sub-129-921-R1F/ses-260115121231"
wholebrain_tck = f"{SUBJ_DIR}/dwi/wholebrain/wholebrain_1000k.tck"
sift2_file = f"{SUBJ_DIR}/dwi/wholebrain/sift2_weights.txt"
passthrough_dir = f"{SUBJ_DIR}/dwi/wholebrain/passthrough_tracts"

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
    
    # Use file: . OFFSET if specified
    if 'file' in hdr:
        parts = hdr['file'].split()
        if len(parts) >= 2:
            data_offset = int(parts[-1])
    
    return data_offset, hdr

def load_tck_streamline_hashes(filepath, n_vertices_for_hash=3):
    """
    Read .tck binary data and build a hash index mapping
    streamline hash -> streamline index.
    
    Uses vectorized numpy operations for speed.
    Hash uses first n_vertices_for_hash vertices (36 bytes for 3 verts).
    """
    t0 = time.time()
    data_offset, hdr = parse_tck_header(filepath)
    
    # Read all binary data at once
    with open(filepath, 'rb') as f:
        f.seek(data_offset)
        raw_data = np.fromfile(f, dtype='<f4')
    
    data = raw_data.reshape(-1, 3)
    print(f"  Loaded {len(data)} triplets in {time.time()-t0:.1f}s")
    
    # Find NaN separators (streamline boundaries) - vectorized
    t1 = time.time()
    nan_mask = np.isnan(data[:, 0])
    nan_indices = np.where(nan_mask)[0]
    inf_mask = np.isinf(data[:, 0])
    
    # Streamline starts: position 0, then after each NaN
    starts = np.concatenate([[0], nan_indices + 1])
    
    # Filter: starts within data range 
    starts = starts[starts < len(data)]
    
    # Filter: skip Inf terminators
    starts = starts[~inf_mask[starts]]
    
    # Filter: at least n_vertices_for_hash vertices before next NaN
    end_idx = np.searchsorted(nan_indices, starts)
    end_positions = np.where(end_idx < len(nan_indices), nan_indices[end_idx], len(data))
    n_verts = end_positions - starts
    valid_starts = starts[n_verts >= n_vertices_for_hash]
    
    print(f"  Found {len(valid_starts)} valid streamlines in {time.time()-t1:.1f}s")
    
    # Build hash map using raw bytes of first 3 vertices
    t2 = time.time()
    flat = data.ravel()  # flat view, no copy
    n_floats = n_vertices_for_hash * 3  # 9 floats = 36 bytes
    
    hash_map = {}
    for i, s in enumerate(valid_starts):
        key = flat[s * 3 : s * 3 + n_floats].tobytes()
        hash_map[key] = i
    
    print(f"  Built hash map in {time.time()-t2:.1f}s")
    print(f"  Total: {len(hash_map)} streamlines indexed in {time.time()-t0:.1f}s")
    return hash_map, len(valid_starts)

def load_passthrough_hashes(filepath, n_vertices_for_hash=3):
    """Load a passthrough tract and return hashes (vectorized)."""
    data_offset, _ = parse_tck_header(filepath)
    
    with open(filepath, 'rb') as f:
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
    valid_starts = starts[n_verts >= n_vertices_for_hash]
    
    flat = data.ravel()
    n_floats = n_vertices_for_hash * 3
    hashes = []
    for s in valid_starts:
        key = flat[s * 3 : s * 3 + n_floats].tobytes()
        hashes.append(key)
    
    return hashes

def load_sift2_weights(filepath):
    """Load SIFT2 weights from text file."""
    t0 = time.time()
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Skip comment lines
    lines = [l for l in content.strip().split('\n') if not l.startswith('#')]
    weights = np.array([float(x) for x in ' '.join(lines).split()], dtype=np.float64)
    print(f"  Loaded {len(weights)} SIFT2 weights in {time.time()-t0:.1f}s")
    return weights

# ====================
# Main test
# ====================

print("=" * 60)
print("Test: Fast SIFT2 weight extraction via hash matching")
print("=" * 60)

# Step 1: Load wholebrain hash index
print("\n1. Loading wholebrain streamline hash index...")
hash_map, n_wb_streamlines = load_tck_streamline_hashes(wholebrain_tck)
print(f"   → {n_wb_streamlines} streamlines indexed")

# Step 2: Load SIFT2 weights
print("\n2. Loading SIFT2 weights...")
weights = load_sift2_weights(sift2_file)
print(f"   → {len(weights)} weights loaded")

# Step 3: Match a few passthrough tracts
print("\n3. Matching passthrough tracts...")
passthrough_files = sorted([
    f for f in os.listdir(passthrough_dir)
    if f.endswith('.tck')
])[:5]  # test first 5

total_time = time.time()
for tck_name in passthrough_files:
    tck_path = os.path.join(passthrough_dir, tck_name)
    t0 = time.time()
    
    hashes = load_passthrough_hashes(tck_path)
    
    matched_indices = []
    unmatched = 0
    for h in hashes:
        idx = hash_map.get(h)
        if idx is not None:
            matched_indices.append(idx)
        else:
            unmatched += 1
    
    if matched_indices:
        tract_weights = weights[matched_indices]
        w_sum = tract_weights.sum()
        w_mean = tract_weights.mean()
        n_matched = len(matched_indices)
    else:
        w_sum = 0
        w_mean = 0
        n_matched = 0
    
    elapsed = time.time() - t0
    print(f"  {tck_name}")
    print(f"    Streamlines: {len(hashes)} total, {n_matched} matched, {unmatched} unmatched")
    print(f"    SIFT2 sum: {w_sum:.4f}, mean: {w_mean:.6f}")
    print(f"    Time: {elapsed:.2f}s")

print(f"\nTotal matching time: {time.time()-total_time:.1f}s (for {len(passthrough_files)} tracts)")
print("DONE")
