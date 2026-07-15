# -*- coding: utf-8 -*-
"""
Pareto-front scatter plot for the tape-hinge MBO.

Called from run_mbo.py at the end of every batch iteration.
Writes:
  plot/pareto_latest.png        -- always overwritten
  plot/pareto_iter_<NNN>.png    -- history snapshot per batch

Axes (both objectives are MAXIMIZED):
  x = 1 - MAX_FI_ALL
  y = MAX_SE
Reference point: (0.0, 400.0)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless backend
import matplotlib.pyplot as plt

# Pull REF_POINT from gp_model if available; otherwise fall back to the
# spec value so plotter can be imported / tested without torch+botorch.
try:
    from gp_model import REF_POINT
except Exception:
    REF_POINT = [0.0, 400.0]


def _pareto_mask(Y):
    """Boolean mask of non-dominated rows (both objectives maximize)."""
    n = Y.shape[0]
    if n == 0:
        return np.zeros(0, dtype=bool)
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        if not mask[i]:
            continue
        dom = ((Y[:, 0] >= Y[i, 0]) & (Y[:, 1] >= Y[i, 1])
               & ((Y[:, 0] > Y[i, 0]) | (Y[:, 1] > Y[i, 1])))
        if dom.any():
            mask[i] = False
    return mask


def plot_pareto(Y_init, Y_bo, hv, iter_idx, plot_dir):
    """
    Y_init : (n_init, 2) numpy array, initial LHS objectives (maximize form)
    Y_bo   : (n_bo,   2) numpy array, BO objectives (maximize form)
    hv     : current dominated hypervolume vs REF_POINT
    iter_idx : batch counter, used in history filename
    """
    os.makedirs(plot_dir, exist_ok=True)
    Y_init = np.asarray(Y_init, dtype=float).reshape(-1, 2)
    Y_bo   = np.asarray(Y_bo,   dtype=float).reshape(-1, 2)
    Y_all  = np.vstack([Y_init, Y_bo]) if Y_bo.size else Y_init

    # Pareto front is computed only on the FEASIBLE region (above the
    # reference point in both objectives). This is the same domain that
    # Hypervolume.compute() integrates over, so the red Pareto line is
    # consistent with the reported HV.
    if Y_all.size:
        feasible = ((Y_all[:, 0] >= REF_POINT[0])
                    & (Y_all[:, 1] >= REF_POINT[1]))
        Y_feas = Y_all[feasible]
    else:
        Y_feas = Y_all

    if Y_feas.size:
        mask = _pareto_mask(Y_feas)
        Y_par = Y_feas[mask]
        order = np.argsort(Y_par[:, 0])
        Y_par_s = Y_par[order]
    else:
        Y_par   = np.empty((0, 2))
        Y_par_s = Y_par

    fig, ax = plt.subplots(figsize=(8, 6))

    # Shaded region above the reference point (dominated HV area).
    xlim = [
        min(Y_all[:, 0].min() if Y_all.size else REF_POINT[0], REF_POINT[0]) - 0.05,
        max(Y_all[:, 0].max() if Y_all.size else REF_POINT[0], REF_POINT[0]) + 0.05,
    ]
    ylim = [
        min(Y_all[:, 1].min() if Y_all.size else REF_POINT[1], REF_POINT[1]) - 50.0,
        max(Y_all[:, 1].max() if Y_all.size else REF_POINT[1], REF_POINT[1]) + 50.0,
    ]
    ax.axvline(REF_POINT[0], color='#aaaaaa', lw=0.8, ls=':')
    ax.axhline(REF_POINT[1], color='#aaaaaa', lw=0.8, ls=':')

    # All points
    if Y_init.size:
        ax.scatter(Y_init[:, 0], Y_init[:, 1],
                   c='#9ca3af', s=44, alpha=0.85,
                   edgecolors='black', linewidths=0.4,
                   label='Initial LHS (n=%d)' % len(Y_init))
    if Y_bo.size:
        ax.scatter(Y_bo[:, 0], Y_bo[:, 1],
                   c='#1d4ed8', s=46, alpha=0.9,
                   edgecolors='black', linewidths=0.4,
                   label='BO (n=%d)' % len(Y_bo))

    # Pareto front
    if Y_par.size:
        ax.scatter(Y_par[:, 0], Y_par[:, 1],
                   facecolors='none', edgecolors='#dc2626',
                   s=160, linewidths=2.0,
                   label='Pareto (n=%d)' % len(Y_par))
        ax.step(Y_par_s[:, 0], Y_par_s[:, 1],
                where='post', color='#dc2626', alpha=0.55, lw=1.3)

    # Reference point
    ax.scatter([REF_POINT[0]], [REF_POINT[1]],
               marker='*', s=220, c='black',
               label='Reference (%g, %g)' % tuple(REF_POINT))

    ax.set_xlabel('1 - MAX_FI_ALL  (maximize)')
    ax.set_ylabel('MAX_SE  (maximize)')
    ax.set_title('Tape-hinge MBO Pareto  |  batch %03d  '
                 'n=%d  HV=%.4f' % (iter_idx, len(Y_all), hv))
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower left', fontsize=9, framealpha=0.95)
    plt.tight_layout()

    latest  = os.path.join(plot_dir, 'pareto_latest.png')
    history = os.path.join(plot_dir, 'pareto_iter_%03d.png' % iter_idx)
    fig.savefig(latest,  dpi=130)
    fig.savefig(history, dpi=130)
    plt.close(fig)
    return latest, history
