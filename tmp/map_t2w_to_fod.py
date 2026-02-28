#!/usr/bin/env python3
"""Map T2W voxel indices (axis 0) to FOD space for fixel analysis."""
import numpy as np

subjects = {
    'CTRL': {
        't2w_dims': [46, 134, 160], 't2w_vox': [0.3, 0.1, 0.1],
        't2w_xfm': [[1, 0.005751, -0.003826, -7.72],
                     [-0.006366, 0.9823, -0.1875, -1.288],
                     [0.00268, 0.1875, 0.9823, -6.534]],
        'fod_dims': [28, 100, 107], 'fod_vox': [0.5, 0.15, 0.149533],
        'fod_xfm': [[0.9999, 0.01531, 0.0005308, -7.865],
                     [-0.01498, 0.9843, -0.1756, -2.067],
                     [-0.003211, 0.1756, 0.9845, -6.489]],
        't2w_slice': 24,
    },
    'AP': {
        't2w_dims': [48, 134, 180], 't2w_vox': [0.3, 0.1, 0.1],
        't2w_xfm': [[1, 0.00521, -0.006078, -8.773],
                     [-0.005692, 0.9966, -0.08225, -2.881],
                     [0.005629, 0.08228, 0.9966, -5.436]],
        'fod_dims': [30, 100, 120], 'fod_vox': [0.5, 0.15, 0.15],
        'fod_xfm': [[1, 0.006495, -0.003082, -8.936],
                     [-0.006709, 0.9972, -0.07514, -3.844],
                     [0.002585, 0.07516, 0.9972, -5.362]],
        't2w_slice': 25,
    },
    'APTR': {
        't2w_dims': [46, 134, 160], 't2w_vox': [0.3, 0.1, 0.1],
        't2w_xfm': [[0.9998, 0.007131, 0.01721, -7.112],
                     [-0.003861, 0.9831, -0.183, -1.232],
                     [-0.01823, 0.1829, 0.983, -6.347]],
        'fod_dims': [26, 100, 107], 'fod_vox': [0.5, 0.15, 0.149533],
        'fod_xfm': [[0.9996, 0.02445, 0.01566, -6.792],
                     [-0.02143, 0.9851, -0.1704, -2.031],
                     [-0.01959, 0.17, 0.9853, -6.373]],
        't2w_slice': 24,
    },
}

for name, s in subjects.items():
    R_t2w = np.array(s['t2w_xfm'])[:, :3]
    t_t2w = np.array(s['t2w_xfm'])[:, 3]
    D_t2w = np.diag(s['t2w_vox'])

    # T2W voxel at center of image for y,z
    t2w_vox = np.array([s['t2w_slice'],
                        s['t2w_dims'][1] // 2,
                        s['t2w_dims'][2] // 2], dtype=float)
    world = R_t2w @ D_t2w @ t2w_vox + t_t2w

    R_fod = np.array(s['fod_xfm'])[:, :3]
    t_fod = np.array(s['fod_xfm'])[:, 3]
    D_fod = np.diag(s['fod_vox'])

    fod_vox = np.linalg.inv(D_fod) @ np.linalg.inv(R_fod) @ (world - t_fod)

    print(f"{name}: T2W slice {s['t2w_slice']} -> FOD axis0 = {fod_vox[0]:.2f} "
          f"(round: {round(fod_vox[0])})")
