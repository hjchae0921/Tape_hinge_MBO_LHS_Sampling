# -*- coding: utf-8 -*-
"""
Bezier 6CP C2 cutout generator -- extracted from
`bezier_c2_cutout_sampling.ipynb` so it is importable by the MBO driver.

DV order:  [y1, x1, Dx, y3, yP2, fP3, ecc]
"""
import numpy as np
import pandas as pd
from scipy.special import comb

# Physical constants (mm)
D       = 38.0
R       = D / 2.0
L_CYL   = 220.0
L_HALF  = L_CYL / 2.0

DV_NAMES = ['y1', 'x1', 'Dx', 'y3', 'yP2', 'fP3', 'ecc']

# 7-DV bounds, identical to the notebook.
DV_BOUNDS = np.array([
    [ 2.0,  12.0],   # y1
    [25.0,  35.0],   # x1
    [ 5.0,  26.0],   # Dx
    [ 1.0,   6.5],   # y3
    [ 0.5,  13.0],   # yP2
    [ 0.0,   3.0],   # fP3
    [ 0.0,   0.7],   # ecc
])


def bernstein(i, n, t):
    return comb(n, i, exact=True) * (t ** i) * ((1 - t) ** (n - i))


def generate_bezier_Q1(dv, num_points=100):
    """Quarter cutout curve (x>=0, y>=0). dv: (7,) array-like."""
    y1, x1, Dx, y3, yP2, fP3, _ecc = dv
    x2  = 2 * x1
    x3  = x2 + Dx
    xP3 = x2 + fP3 * Dx
    yP3 = 2 * y3

    cps = np.array([
        [0.0, y1 ],
        [x1,  y1 ],
        [x2,  yP2],
        [xP3, yP3],
        [x3,  y3 ],
        [x3,  0.0],
    ])

    deg = 5
    t = np.linspace(0, 1, num_points)
    pts = np.zeros((num_points, 2))
    for i in range(6):
        B = bernstein(i, deg, t)
        pts[:, 0] += cps[i, 0] * B
        pts[:, 1] += cps[i, 1] * B

    derived = {'x2': x2, 'x3': x3, 'xP3': xP3, 'yP3': yP3}
    return pts, cps, derived


def mirror_to_full(q1):
    """Q1 curve -> closed full cutout (4-fold symmetry)."""
    tr = q1
    tl = q1[::-1].copy(); tl[:, 0] *= -1
    bl = q1.copy();       bl[:, 0] *= -1; bl[:, 1] *= -1
    br = q1[::-1].copy(); br[:, 1] *= -1
    return np.vstack([tl, tr, br[1:], bl[1:], tl[0:1]])


def write_cutout_csv(filepath, dv, num_points=100):
    """
    Write specimen/cutout_points_<N>.csv in the same format the LHS
    notebook produced (X,Y + DVs only on row 0, NaN elsewhere).
    """
    y1, x1, Dx, y3, yP2, fP3, ecc_val = [float(v) for v in dv]
    x2  = 2 * x1
    x3  = x2 + Dx
    xP3 = x2 + fP3 * Dx
    yP3 = 2 * y3

    q1, _, _ = generate_bezier_Q1(dv, num_points=num_points)
    full = mirror_to_full(q1)

    df = pd.DataFrame(full, columns=['X', 'Y'])
    n_rows = len(df)
    nan_tail = [np.nan] * (n_rows - 1)

    df['Eccentricity'] = [ecc_val] + nan_tail
    df['y1']           = [y1]      + nan_tail
    df['x1']           = [x1]      + nan_tail
    df['Dx']           = [Dx]      + nan_tail
    df['y3']           = [y3]      + nan_tail
    df['yP2']          = [yP2]     + nan_tail
    df['fP3']          = [fP3]     + nan_tail
    df['x2_derived']   = [x2]      + nan_tail
    df['x3_derived']   = [x3]      + nan_tail
    df['xP3_derived']  = [xP3]     + nan_tail
    df['yP3_derived']  = [yP3]     + nan_tail

    df.to_csv(filepath, index=False)
    return filepath
