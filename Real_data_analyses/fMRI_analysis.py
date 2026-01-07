#!/usr/bin/env python
# coding: utf-8

import hdf5storage
import matplotlib.pyplot as plt
import nilearn
import numpy as np
from nilearn import datasets
import nibabel as nib
from tqdm import tqdm
from joblib import Parallel, delayed
import main_code_bias_corrs_on_data as bc
from utility_functions import *

# Load Z-score data
mat = hdf5storage.loadmat('Zscores_HCP_SH1000.mat')
ts = np.stack([mat["Zscores_HCP"][i][0] for i in range(len(mat["Zscores_HCP"]))]).transpose((0, 2, 1))
ts_sub = np.stack([mat["Zscores_HCP"][i][0] for i in range(len(mat["Zscores_HCP"]))]).transpose((0, 2, 1))

# Load atlas and extract volume
atlas = datasets.fetch_atlas_schaefer_2018(n_rois=1000, yeo_networks=7, resolution_mm=1)
atlas_vol = nib.load(atlas["maps"]).get_fdata()

# Compute centroids for each region
centroids = {}
for label in range(1, int(max(np.unique(atlas_vol)) + 1)):
    coords = np.argwhere(atlas_vol == label)
    if coords.size == 0:
        continue
    centroid = coords.mean(axis=0)
    centroids[label] = centroid

centroid_list = [centroids.get(label, np.array([np.nan, np.nan, np.nan])) for label in
                 range(1, int(max(np.unique(atlas_vol)) + 1))]
centroid_array = np.stack(centroid_list)

# Process atlas labels
labels = atlas["labels"]
labels_str = [label.decode('utf-8') for label in labels]

lookup = {}
for idx, label in enumerate(labels_str):
    if label.startswith("7Networks_"):
        parts = label.split('_')
        hemisphere = parts[1]
        network = parts[2]
        combination = f"{hemisphere}_{network}"
        lookup.setdefault(combination, []).append(idx)

unique_combinations = list(lookup.keys())
comb_to_code = {comb: i for i, comb in enumerate(unique_combinations)}

code_array = np.empty(len(labels_str), dtype=int)
for idx, label in enumerate(labels_str):
    if label.startswith("7Networks_"):
        parts = label.split('_')
        hemisphere = parts[1]
        network = parts[2]
        combination = f"{hemisphere}_{network}"
        code_array[idx] = comb_to_code[combination]
    else:
        code_array[idx] = -1

yeo_labels = code_array.copy()


# Utility functions
def find_constant_or_nan_columns(matrix):
    matrix = np.asarray(matrix)
    constant_or_nan_cols = []

    for i in range(matrix.shape[1]):
        col = matrix[:, i]
        if np.isnan(col).any() or np.all(col == col[0]):
            constant_or_nan_cols.append(i)

    return constant_or_nan_cols

def get_nonconstant_nonan_column_permutation(matrix, seed=None):
    matrix = np.asarray(matrix)
    np.random.seed(seed)
    valid_cols = []

    for i in range(matrix.shape[1]):
        col = matrix[:, i]
        if not np.isnan(col).any() and not np.all(col == col[0]):
            valid_cols.append(i)

    permuted_cols = np.random.permutation(valid_cols)
    return permuted_cols.tolist()

def process_seed_triplet(seeds, ts, centroid_array, yeo_labels, orders=None, time_length=256):
    if orders is None:
        orders = list(range(1, 21))
    n_orders = len(orders)

    total_timepoints = ts.shape[0]
    delay = 1

    if time_length == 'all':
        time_start = 0
        time_length = total_timepoints - delay
    else:
        time_start = 300

    syns_measures = np.zeros((n_orders, 4))
    reds_measures = np.zeros((n_orders, 4))
    uniques_measures = np.zeros((n_orders, 4))
    joints_measures = np.zeros((n_orders, 4))

    unique_neighbors = get_unique_neighbors(seeds, centroid_array, yeo_labels, max_dim=30)
    if len(np.intersect1d(unique_neighbors[seeds[0]][30], unique_neighbors[seeds[1]][30])) > 0:
        print('There are repeated ROIs in neighbors')

    for n, order in enumerate(orders):
        ts1 = ts[:-delay]
        ts2 = ts[:-delay]
        ts3 = ts[delay:]

        input_data = np.hstack([ts1, ts2, ts3])
        usable_rois = get_nonconstant_nonan_column_permutation(input_data[time_start:time_start + time_length, :])

        ts1 = input_data[time_start:time_start + time_length, usable_rois[:order]]
        ts2 = input_data[time_start:time_start + time_length, usable_rois[order:2*order]]
        ts3 = input_data[time_start:time_start + time_length, usable_rois[2*order:3*order]]
        input_data_small = np.hstack([ts1, ts2, ts3])

        warn_if_constant_or_nan(input_data_small)

        temp = bc.execute_bias_correction_routines(
            input_data_small,
            ["no_bias_correction", "informative_bias_correction", "shuffsub_bias_correction",
             "Venkatesh_bias_correction"],
            num_iterations=20,
            parallel=False,
            save=False,
            save_format="matlab",
            cov_method="shrink"
        )

        try:
            syns_measures[n, :] = np.array([temp[i][1][0]["si"] for i in range(4)])
            reds_measures[n, :] = np.array([temp[i][1][0]["ri"] for i in range(4)])
            uniques_measures[n, :] = np.array([temp[i][1][0]["uix"] + temp[i][1][0]["uiy"] for i in range(4)])
            joints_measures[n, :] = np.array([temp[i][1][0]["imxy"] for i in range(4)])
        except TypeError as e:
            print(f"TypeError at order {order}: {e}")
            syns_measures[n, :] = np.nan
            reds_measures[n, :] = np.nan
            uniques_measures[n, :] = np.nan
            joints_measures[n, :] = np.nan

    return {
        'seeds': seeds,
        'syns_measures': syns_measures,
        'reds_measures': reds_measures,
        'uniques_measures': uniques_measures,
        'joints_measures': joints_measures,
    }


# === FINAL RUN ===
subs = list(range(100))
orders_to_compute = np.arange(1, 31, step=1)
time_windows = [64, 128, 256, 512, 'all']
seeds_triplets = generate_triplets(yeo_labels, n_triplets=30, first_network=6)  # VIS Net = 6

N_time = len(time_windows)
N_triplets = len(seeds_triplets)
N_orders = len(orders_to_compute)
N_methods = 4

Nsyns_measures_tot = np.zeros((N_triplets, len(subs), N_orders, N_methods, N_time))
Nreds_measures_tot = np.zeros_like(Nsyns_measures_tot)
Nuniques_measures_tot = np.zeros_like(Nsyns_measures_tot)
Njoints_measures_tot = np.zeros_like(Nsyns_measures_tot)

for N, seeds in enumerate(seeds_triplets):
    print(f"\nProcessing seed triplet {N + 1}/{N_triplets}")

    for t_idx, t_win in enumerate(time_windows):
        print(f"  Time window: {t_win}")

        results = Parallel(n_jobs=-1)(
            delayed(process_seed_triplet)(
                seeds,
                ts_sub[s],
                centroid_array,
                yeo_labels,
                orders_to_compute,
                time_length=t_win
            ) for s in tqdm(subs, desc=f"Triplet {N + 1}, TimeWin {t_win}")
        )

        Nsyns_measures_tot[N, :, :, :, t_idx] = np.stack([res['syns_measures'] for res in results])
        Nreds_measures_tot[N, :, :, :, t_idx] = np.stack([res['reds_measures'] for res in results])
        Nuniques_measures_tot[N, :, :, :, t_idx] = np.stack([res['uniques_measures'] for res in results])
        Njoints_measures_tot[N, :, :, :, t_idx] = np.stack([res['joints_measures'] for res in results])

# Save all outputs
np.savez('pid_vals_allwindows_shrink_random_shortdelay.npz',
         Nsyns_measures_tot=Nsyns_measures_tot,
         Nreds_measures_tot=Nreds_measures_tot,
         Nuniques_measures_tot=Nuniques_measures_tot,
         Njoints_measures_tot=Njoints_measures_tot,
         time_windows=time_windows
         )
