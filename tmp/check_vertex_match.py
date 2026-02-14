#!/usr/bin/env python3
"""Check coordinate differences between raw and MATLAB-cleaned .tck files."""
import struct, numpy as np

def read_all_streamlines(path, max_streamlines=500):
    with open(path, 'rb') as f:
        content = f.read()
    idx = content.find(b'END\n')
    if idx < 0:
        idx = content.find(b'END\r\n')
    header = content[:idx+4].decode('latin-1')
    offset = None
    for line in header.split('\n'):
        if line.startswith('file:'):
            offset = int(line.split()[-1])
            break
    if offset is None:
        offset = idx + 4
    data = content[offset:]
    n_floats = len(data) // 4
    floats = struct.unpack('<%df' % n_floats, data)
    streamlines = []
    current = []
    for i in range(0, n_floats, 3):
        if i + 2 >= n_floats:
            break
        x, y, z = floats[i], floats[i+1], floats[i+2]
        if np.isinf(x):
            if current:
                streamlines.append(current)
            break
        if np.isnan(x):
            if current:
                streamlines.append(current)
                current = []
                if len(streamlines) >= max_streamlines:
                    break
            continue
        current.append((x, y, z))
    return streamlines

raw = 'data/derivatives/mrtrix_preproc_flip1/sub-596-C-OO-NmM/ses-240911103349/dwi/wholebrain/passthrough_tracts/seg9_graft_left_anterior_dorsal_to_Thalamus_left.tck'
clean = raw.replace('.tck', '_cleaned.tck')

r = read_all_streamlines(raw, 250)
c = read_all_streamlines(clean, 100)
rl = [len(s) for s in r]
cl = [len(s) for s in c]
print(f'Raw: {len(r)} streamlines, lengths: min={min(rl)}, max={max(rl)}, mean={np.mean(rl):.1f}')
print(f'Clean: {len(c)} streamlines, lengths: min={min(cl)}, max={max(cl)}, mean={np.mean(cl):.1f}')

# Now try matching WITHOUT length filter
for ci in range(min(5, len(c))):
    csl = c[ci]
    cv = np.array(csl[:3])
    best_dist = float('inf')
    best_idx = -1
    for ri in range(len(r)):
        rsl = r[ri]
        rv = np.array(rsl[:min(3, len(rsl))])
        if rv.shape != cv.shape:
            continue
        dist = np.sqrt(np.sum((cv - rv) ** 2))
        if dist < best_dist:
            best_dist = dist
            best_idx = ri
    if best_idx >= 0:
        print(f'  clean[{ci}] len={len(csl)}: closest raw[{best_idx}] len={len(r[best_idx])}, dist={best_dist:.6f}')
    else:
        print(f'  clean[{ci}] len={len(csl)}: no match found')
