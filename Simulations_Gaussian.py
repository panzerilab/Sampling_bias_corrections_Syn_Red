#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import print_function, division
import warnings
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy.linalg as la
from sklearn.model_selection import KFold
# from gpid.generate import swap_x_and_y
from gpid.tilde_pid import exact_gauss_tilde_pid
from gpid.utils import whiten, solve
import numpy.linalg as npla
import concurrent.futures
import seaborn as sns
import scipy.io as sio
import os
import time
from scipy.stats import wishart, ortho_group
from joblib import Parallel, delayed
import tools
import multiprocessing

# ===============================
# Parallelized job function
# ===============================

def run_simulation(ntrials, M, T, n_repetitions, alpha, bias_title, routine, case_type):
    """
    Run one experiment with the given parameters and write results.
    """
    print(f"Running: ntrials={ntrials}, alpha={alpha}, M={M}")
    N = M * 3
    dm = dx = dy = M
    # Generate covariance matrix and data
    p=q=0.1 #alpha
    r=0
    cov, dm, dx, dy = tools.gen_cov_matrix(N, M, p, q, r, case_type)
    np.fill_diagonal(cov, alpha * np.diagonal(cov))

    # Compute GT terms for entropy and PID (full and conditionally independent)
    HM_GT = 0.5 / np.log(2) * (la.slogdet(cov[:dm, :dm])[1] + dm * (1 + np.log(2 * np.pi)))
    HXY_GT = 0.5 / np.log(2) * (la.slogdet(cov[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY_GT = 0.5 / np.log(2) * (la.slogdet(cov)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    ret = exact_gauss_tilde_pid(cov, dm, dx, dy)
    imxy_GT, uix_GT, uiy_GT, ri_GT, si_GT = ret[2], *ret[-4:]
    imx_GT = uix_GT + ri_GT
    imy_GT = uiy_GT + ri_GT

    # Build the conditionally independent (IND) version
    cov_ind = np.copy(cov)
    sig_m = cov[:dm, :dm]
    sig_xy_m = cov[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
    cov_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]
    ret_ind = exact_gauss_tilde_pid(cov_ind, dm, dx, dy)
    imxy_ind_GT, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]
    HM_ind_GT = 0.5 / np.log(2) * (la.slogdet(cov_ind[:dm, :dm])[1] + dm * (1 + np.log(2 * np.pi)))
    HXY_ind_GT = 0.5 / np.log(2) * (la.slogdet(cov_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY_ind_GT = 0.5 / np.log(2) * (la.slogdet(cov_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

    # n_occurrences = 100
    n_occurrences = n_repetitions
    seed_list = np.arange(n_repetitions).tolist()   # Seeds from 10 to 100 with step 10

    imxy_iter = np.zeros(n_occurrences)
    uix_iter = np.zeros(n_occurrences)
    uiy_iter = np.zeros(n_occurrences)
    imx_iter = np.zeros(n_occurrences)
    imy_iter = np.zeros(n_occurrences)
    ri_iter = np.zeros(n_occurrences)
    si_iter = np.zeros(n_occurrences)
    HM_iter = np.zeros(n_occurrences)
    HXY_iter = np.zeros(n_occurrences)
    HMXY_iter = np.zeros(n_occurrences)
    imxy_ind_iter = np.zeros(n_occurrences)
    HM_ind_iter = np.zeros(n_occurrences)
    HXY_ind_iter = np.zeros(n_occurrences)
    HMXY_ind_iter = np.zeros(n_occurrences)

    # biased_results = []
    # GT_results = []
    for i, seed in enumerate(seed_list):
        # print(f"Routine: {bias_title}, Seed index: {i}, Seed: {seed}")
        rng = np.random.default_rng(seed)
        input_data = rng.multivariate_normal(np.zeros(cov.shape[0]), cov, size=ntrials)

        # Run the specified bias-correction routine
        (imxy, uix, uiy, imx, imy, ri, si,
         HM, HXY, HMXY, imxy_ind, HM_ind, HXY_ind, HMXY_ind) = routine(input_data, T, ntrials, seed)

        imxy_iter[i] = imxy
        uix_iter[i] = uix
        uiy_iter[i] = uiy
        imx_iter[i] = imx
        imy_iter[i] = imy
        ri_iter[i] = ri
        si_iter[i] = si
        HM_iter[i] = HM
        HXY_iter[i] = HXY
        HMXY_iter[i] = HMXY
        imxy_ind_iter[i] = imxy_ind
        HM_ind_iter[i] = HM_ind
        HXY_ind_iter[i] = HXY_ind
        HMXY_ind_iter[i] = HMXY_ind

    sampled_results = np.stack((imxy_iter, uix_iter, uiy_iter, imx_iter, imy_iter, ri_iter, si_iter,
         HM_iter, HXY_iter, HMXY_iter, imxy_ind_iter, HM_ind_iter, HXY_ind_iter, HMXY_ind_iter))
    gt_results = np.array([
        imxy_GT, uix_GT, uiy_GT, imx_GT, imy_GT,
        ri_GT, si_GT, HM_GT, HXY_GT, HMXY_GT,
        imxy_ind_GT, HM_ind_GT, HXY_ind_GT, HMXY_ind_GT
    ])

    return (sampled_results, gt_results)



def run_task(ntrials, M, T, n_repetitions, alpha, bias_title, routine, case_type):
    try:
        biased_results, GT_results = run_simulation(ntrials, M, T, n_repetitions, alpha, bias_title, routine, case_type)
        return (ntrials, alpha, M, bias_title, case_type, biased_results, GT_results, None)
    except Exception as e:
        return (ntrials, alpha, M, bias_title, case_type, None, None, e)


def main():
    # These simulations run the necessary simulations for all figures in 
    # the paper.
    # The first block does the simulations for all figures except 
    # the ones for dimension sweep until d=80, those are run after in the
    # same script and saved with the distinction d80 at the end of the
    # filename.
    # If you only want the simulations for the main figures, only out 
    # 'bit_of_all' in the list

    start_time = time.time()  # Record the start time

    # Parameters
    ntrials = [64, 128, 256, 512, 1024, 2058]
    M_vals_ALL = np.arange(4, 43, step=2)
    M_subset = [10, 20, 30]
    trial_subset = [64, 128, 256]
    alphas = np.array([1., 1.1, 1.2, 1.5, 2., 10000.])
    n_repetitions = 100
    T = 20  # Number of iterations for bias estimation

    # Define the bias-correction routines to run
    bias_routines = [
        ("No_bias_correction", tools.No_bias_correction_routine),
        ("informative_bias_correction", tools.informative_bias_correction_routine),
        ("shuffle_subtr", tools.shuffle_subtr_bias_correction_routine),
        ("Venka_bias_correction", tools.uniform_bias_correction_routine)
    ]
    bias_titles = [title for title, _ in bias_routines]

    # Collect all used M and ntrials values
    M_vals_set = sorted(set(M_subset + list(M_vals_ALL)))
    ntrials_set = sorted(set(ntrials + trial_subset))

    # Define all case types
    cases = ['bit_of_all', 'both_unique', 'high_synergy', 'zero_synergy']
    M_max_map = {
        64: 20,
        128: 42,
        256: 86
    }
    for ca in cases:
        print(f"\n===== Starting case: {ca} =====")
        tasks = []

        # Case 1: M in [10, 20, 30], iterate all ntrials
        for M_val in M_subset:
            for ntr in ntrials:
                for alpha in alphas:
                    for bias_title, routine in bias_routines:
                        tasks.append((ntr, M_val, T, n_repetitions, alpha, bias_title, routine, ca))

        # Case 2: ntrials in [64, 128, 256], iterate only up to specific M values
        for ntr in trial_subset:
            max_M = M_max_map[ntr]
            for M_val in M_vals_ALL:
                if M_val > max_M or M_val in M_subset:
                    continue  # Skip if above allowed dimension or already covered in subset
                for alpha in alphas:
                    for bias_title, routine in bias_routines:
                        tasks.append((ntr, M_val, T, n_repetitions, alpha, bias_title, routine, ca))

        # Run tasks in parallel
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(run_task)(*task) for task in tasks
        )

        # Initialize result arrays
        biased_results_ALL = np.zeros((
            14, len(ntrials_set), len(M_vals_set), len(bias_routines), len(alphas), n_repetitions
        ))
        GT_results_ALL = np.zeros((
            14, len(ntrials_set), len(M_vals_set), len(bias_routines), len(alphas)
        ))

        # Collect results
        for (ntr, alpha, M, bias_title, case_type, biased_results, GT_results, error) in results:
            if error is not None:
                print(f"Error in job for (case={case_type}, ntrials={ntr}, alpha={alpha}, M={M}, bias_title={bias_title}): {error}")
                continue

            idx_ntr = np.where(np.array(ntrials_set) == ntr)[0][0]
            idx_m = np.where(np.array(M_vals_set) == M)[0][0]
            idx_alpha = np.where(alphas == alpha)[0][0]
            idx_corr = bias_titles.index(bias_title)

            biased_results_ALL[:, idx_ntr, idx_m, idx_corr, idx_alpha, :] = biased_results
            GT_results_ALL[:, idx_ntr, idx_m, idx_corr, idx_alpha] = GT_results

        # Save results
        case_dir = "Results/Gauss/"
        os.makedirs(case_dir, exist_ok=True)
        if ca != 'bit_of_all':
            save_path = os.path.join(case_dir, f"Finalresults_across_M_and_ntrials_{ca}.mat")
        else:
            save_path = os.path.join(case_dir, f"Finalresults_across_M_and_ntrials.mat")

        sio.savemat(save_path, {
            'sampled_results': biased_results_ALL,
            'GT_results': GT_results_ALL,
            'M_vals': M_vals_set,
            'ntrials_vals': ntrials_set,
            'alphas': alphas,
            'bias_titles': bias_titles
        })

        print(f"Saved results for case '{ca}' at {save_path}")

    M_vals_ALL = np.arange(4, 83, step=4)
    M_subset = [10, 20, 30]
    trial_subset = [256]
    for ca in ['bit_of_all']:
        print(f"\n===== Starting case: {ca} =====")
        tasks = []

        # Case 1: M in [10, 20, 30], iterate all ntrials
        for M_val in M_subset:
            for ntr in ntrials:
                for alpha in alphas:
                    for bias_title, routine in bias_routines:
                        tasks.append((ntr, M_val, T, n_repetitions, alpha, bias_title, routine, ca))

        # Case 2: ntrials in [64, 128, 256], iterate only up to specific M values
        for ntr in trial_subset:
            max_M = M_max_map[ntr]
            for M_val in M_vals_ALL:
                if M_val > max_M or M_val in M_subset:
                    continue  # Skip if above allowed dimension or already covered in subset
                for alpha in alphas:
                    for bias_title, routine in bias_routines:
                        tasks.append((ntr, M_val, T, n_repetitions, alpha, bias_title, routine, ca))

        # Run tasks in parallel
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(run_task)(*task) for task in tasks
        )

        # Initialize result arrays
        biased_results_ALL = np.zeros((
            14, len(ntrials_set), len(M_vals_set), len(bias_routines), len(alphas), n_repetitions
        ))
        GT_results_ALL = np.zeros((
            14, len(ntrials_set), len(M_vals_set), len(bias_routines), len(alphas)
        ))

        # Collect results
        for (ntr, alpha, M, bias_title, case_type, biased_results, GT_results, error) in results:
            if error is not None:
                print(f"Error in job for (case={case_type}, ntrials={ntr}, alpha={alpha}, M={M}, bias_title={bias_title}): {error}")
                continue

            idx_ntr = np.where(np.array(ntrials_set) == ntr)[0][0]
            idx_m = np.where(np.array(M_vals_set) == M)[0][0]
            idx_alpha = np.where(alphas == alpha)[0][0]
            idx_corr = bias_titles.index(bias_title)

            biased_results_ALL[:, idx_ntr, idx_m, idx_corr, idx_alpha, :] = biased_results
            GT_results_ALL[:, idx_ntr, idx_m, idx_corr, idx_alpha] = GT_results

        # Save results
        case_dir = "Results/Gauss/"
        os.makedirs(case_dir, exist_ok=True)
        
        save_path = os.path.join(case_dir, f"Finalresults_across_M_and_ntrials_d80.mat")

        sio.savemat(save_path, {
            'sampled_results': biased_results_ALL,
            'GT_results': GT_results_ALL,
            'M_vals': M_vals_set,
            'ntrials_vals': ntrials_set,
            'alphas': alphas,
            'bias_titles': bias_titles
        })

        print(f"Saved results for case '{ca}' at {save_path}")


    # Print total runtime
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nTotal computed time: {total_time:.2f} seconds")

    return


if __name__ == "__main__":
    main()
