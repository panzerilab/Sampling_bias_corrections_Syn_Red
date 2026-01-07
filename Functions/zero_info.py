#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gaussian PID histograms with GPID — parallel, no argparse.

- Spawns worker processes (multiprocessing) to run simulations in parallel
- Each worker:
  * draws N_TRIALS samples from N(0, cov_true) over [S, R1, R2]
  * computes sample covariance
  * calls exact_gauss_tilde_pid to get MI/PID atoms
- Aggregates results and plots normalized line-histograms:
  Joint, Single (pooled), Syn, Red, Unq (pooled)
- Saves arrays (.npz) and the figure (.svg)
"""

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count, get_start_method
from numpy.random import SeedSequence, default_rng
from gpid.tilde_pid import exact_gauss_tilde_pid

# ======================
# --- USER SETTINGS ----
# ======================
N_TRIALS     = 64       # samples per repetition
N_SIMS       = 10000     # repetitions (each is one PID from a fresh sample covariance)
D            = 10         # dimension per block (S, R1, R2)
CASE         = "independent"   # one of: "independent", "shared_info", "high_synergy"
FIGURES_DIR  = "Figures"
RESULTS_DIR  = "Results"
RNG_SEED     = 20        # master seed (reproducible across runs)
XMAX_BITS    = 0.5       # x-axis upper limit in bits
NBINS        = 40        # histogram bins
N_WORKERS    = None      # None -> use all cores, or set an int (e.g., 8)

# ======================
# Covariance builders
# ======================
def cov_independent(d: int) -> np.ndarray:
    """S, R1, R2 all independent standard Gaussian."""
    return np.eye(3 * d)

def cov_shared_s_info(d: int, rho: float = 0.6) -> np.ndarray:
    """Redundancy-ish: S correlated with R1 and R2, but R1 and R2 not correlated directly."""
    I = np.eye(d)
    Z = np.zeros((d, d))
    top = np.concatenate([I,       rho*I,   rho*I], axis=1)
    mid = np.concatenate([rho*I,   I,       Z    ], axis=1)
    bot = np.concatenate([rho*I,   Z,       I    ], axis=1)
    return np.concatenate([top, mid, bot], axis=0)

def cov_high_synergy(d: int, rho: float = 0.7) -> np.ndarray:
    """
    Toy 'synergy-ish' example: no direct S–R1/R2 correlation, but R1–R2 cross-correlation.
    (Heuristic; synergy here is illustrative, not guaranteed.)
    """
    C = np.eye(3 * d)
    C[d:2*d, 2*d:3*d] = rho * np.eye(d)
    C[2*d:3*d, d:2*d] = rho * np.eye(d)
    return C

CASES = {
    "independent":  cov_independent,
    "shared_info":  cov_shared_s_info,
    "high_synergy": cov_high_synergy,
}

# ======================
# Plot helper
# ======================
def line_hist(ax, data, nbins, color, label):
    x = np.asarray(data).ravel()
    x = x[np.isfinite(x)]
    if x.size == 0:
        return
    f, edges = np.histogram(x, bins=nbins)  # let NumPy pick range
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = f.sum()
    if total > 0:
        ax.plot(centers, f / total, lw=2, color=color, label=label)

# ======================
# Parallel worker setup
# ======================
# We'll stash read-only state in module-level globals so each worker
# doesn't get huge argument copies each call.
_G = {
    "cov_true": None,
    "dm": None, "dx": None, "dy": None,
    "n_trials": None,
}

def _worker_init(cov_true, dm, dx, dy, n_trials):
    _G["cov_true"] = cov_true
    _G["dm"] = dm
    _G["dx"] = dx
    _G["dy"] = dy
    _G["n_trials"] = n_trials

def _simulate_one(seed_int: int):
    """
    One repetition using its own RNG stream.
    Returns tuple:
      (imxy, imx, imy, ri, uix, uiy, si)
    All in bits (as returned by exact_gauss_tilde_pid).
    """
    rng = default_rng(seed_int)
    cov_true = _G["cov_true"]
    dm, dx, dy = _G["dm"], _G["dx"], _G["dy"]
    n_trials = _G["n_trials"]

    data = rng.multivariate_normal(mean=np.zeros(dm+dx+dy), cov=cov_true, size=n_trials)
    cov_sample = np.cov(data, rowvar=False, bias=False)

    # exact_gauss_tilde_pid expects covariance over [M, X, Y] = [S, R1, R2]
    ret = exact_gauss_tilde_pid(cov_sample, dm, dx, dy)

    imxy = ret[2]              # I(M;X,Y)
    uix, uiy, ri, si = ret[-4:]  # UIx, UIy, RI, SI
    imx = ret[0]              #uix + ri
    imy = ret[1]              #uiy + ri
    return (imxy, imx, imy, ri, uix, uiy, si)

# ======================
# Main
# ======================
def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    dm = dx = dy = D
    cov_true = CASES[CASE](D)

    # Prepare worker initializer and seeds (reproducible, independent streams)
    ss = SeedSequence(RNG_SEED)
    # Generate one uint64 seed per simulation
    seeds = ss.generate_state(N_SIMS, dtype=np.uint64)

    # Decide worker count
    n_workers = N_WORKERS or cpu_count()
    # macOS uses 'spawn' by default; that's fine. Just guarding for clarity:
    # print("Start method:", get_start_method())

    # Run in parallel
    try:
        with Pool(processes=n_workers,
                  initializer=_worker_init,
                  initargs=(cov_true, dm, dx, dy, N_TRIALS)) as pool:
            # Use imap for streaming results without holding everything in RAM at once
            # Chunk size heuristic:
            chunk = max(1, math.ceil(N_SIMS / (n_workers * 8)))
            results_iter = pool.imap(_simulate_one, seeds, chunksize=chunk)

            # Preallocate
            MI_joint   = np.empty(N_SIMS, dtype=float)
            MI_single  = np.empty((2, N_SIMS), dtype=float)
            PID_atoms  = np.empty((4, N_SIMS), dtype=float)  # [RI, UIx, UIy, SI]

            for i, tup in enumerate(results_iter):
                imxy, imx, imy, ri, uix, uiy, si = tup
                MI_joint[i]     = imxy
                MI_single[0, i] = imx
                MI_single[1, i] = imy
                PID_atoms[:, i] = (ri, uix, uiy, si)

    except Exception as e:
        # Fallback to serial if something odd happens (rare)
        print(f"[parallel run failed: {e}] — falling back to serial.", file=sys.stderr)
        rng = np.random.default_rng(RNG_SEED)
        MI_joint   = np.zeros(N_SIMS)
        MI_single  = np.zeros((2, N_SIMS))
        PID_atoms  = np.zeros((4, N_SIMS))
        for i in range(N_SIMS):
            data = rng.multivariate_normal(np.zeros(dm+dx+dy), cov_true, size=N_TRIALS)
            cov_sample = np.cov(data, rowvar=False, bias=False)
            ret = exact_gauss_tilde_pid(cov_sample, dm, dx, dy)
            imxy = ret[2]
            uix, uiy, ri, si = ret[-4:]
            MI_joint[i]     = imxy
            MI_single[0, i] = uix + ri
            MI_single[1, i] = uiy + ri
            PID_atoms[:, i] = (ri, uix, uiy, si)

    # ---------- Plot ----------
    # Colors (matching your MATLAB flavor)
    color_joint  = np.array([92.9, 70.2, 14.5]) / 100.0  # brown-ish
    color_syn    = np.array([47.8, 67.8, 20.4]) / 100.0  # green-ish
    color_unq    = np.array([50.2, 50.2, 50.2]) / 100.0  # dark gray
    color_red    = np.array([0.0, 0.4470, 0.7410])       # blue
    color_single = np.array([0.3, 0.3, 0.3])             # gray

    singles_all = MI_single.ravel()
    uniques_all = PID_atoms[1:3, :].ravel()  # UIx, UIy pooled
    joint_all   = MI_joint
    syn_all     = PID_atoms[3, :]
    red_all     = PID_atoms[0, :]

    fig, ax = plt.subplots(figsize=(8, 5))
    line_hist(ax, joint_all,   NBINS, color_joint,  "Joint")
    line_hist(ax, singles_all, NBINS, color_single, "Single")
    line_hist(ax, syn_all,     NBINS, color_syn,    "Syn")
    line_hist(ax, red_all,     NBINS, color_red,    "Red")
    line_hist(ax, uniques_all, NBINS, color_unq,    "Unq")

    ax.set_xlim(-0.01, XMAX_BITS)
    ax.set_xlabel("Information [bits]", fontsize=12)
    ax.set_ylabel("Normalized frequency", fontsize=12)
    ax.legend(frameon=False, fontsize=11)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    fig_path = os.path.join(FIGURES_DIR, f"Gaussian_PID_hist_{CASE}{N_TRIALS}trials_{D}dimensions_parallel.svg")
    fig.savefig(fig_path)
    print(f"Saved histogram: {fig_path}")

    # ---------- Save arrays ----------
    npz_path = os.path.join(RESULTS_DIR, f"Gaussian_PID_{CASE}{N_TRIALS}trials_{D}dimensions_parallel.npz")
    np.savez_compressed(
        npz_path,
        MI_joint=MI_joint,
        MI_single=MI_single,
        PID_atoms=PID_atoms,
        d=D,
        n_trials=N_TRIALS,
        n_sims=N_SIMS,
        case=CASE
    )
    print(f"Saved arrays: {npz_path}")

if __name__ == "__main__":
    main()
