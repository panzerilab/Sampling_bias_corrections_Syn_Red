import hdf5storage
import matplotlib.pyplot as plt
import nilearn
import numpy as np
from nilearn import datasets
import nibabel as nib
from tqdm import tqdm
from joblib import Parallel, delayed
import main_code_bias_corrs_on_data as bc
import warnings
from sklearn.covariance import LedoitWolf


def find_invertible_window(X, t):
    """
    Find a time window of size t where the covariance matrix is invertible.

    Parameters:
    - X: A 2D numpy array with shape (T, d), where T is the total number of time steps,
        and d is the number of features/dimensions.
    - t: The size of the time window.

    Returns:
    - start_index: The start index of the invertible window.
    - end_index: The end index of the invertible window.
    - cov: The covariance matrix of the selected window.
    """
    T, d = X.shape
    assert t >= d, "Window size t must be at least the number of dimensions d"

    for i in range(T - t + 1):
        # Select the window of data
        window = X[i:i + t,:]

        # Calculate the covariance matrix (without normalizing over rows, using `rowvar=False`)
        cov = np.cov(window)

        # Check if the covariance matrix is invertible (rank equals the number of dimensions)
        if np.linalg.matrix_rank(cov) == d:
            return i, i + t, cov  # Return start and end of window, and the invertible covariance matrix

    # If no invertible window is found, return None
    return None, None, None

def get_unique_neighbors(seeds, centroids, yeo_labels, max_dim=30):
    network_seeds = {}
    for seed in seeds:
        net = yeo_labels[seed]
        network_seeds.setdefault(net, []).append(seed)

    results = {seed: {} for seed in seeds}

    for net, seed_list in network_seeds.items():
        candidates = np.where(yeo_labels == net)[0]

        seed_sorted = {}
        for seed in seed_list:
            distances = np.linalg.norm(centroids[candidates] - centroids[seed], axis=1)
            order = np.argsort(distances)
            seed_sorted[seed] = list(candidates) #list(candidates[order])

        assigned = {seed: [] for seed in seed_list}
        for seed in seed_list:
            if seed in seed_sorted[seed]:
                assigned[seed].append(seed)
            else:
                assigned[seed].insert(0, seed)

        for d in range(2, max_dim + 1):
            global_assigned = set(sum(assigned.values(), []))
            for seed in seed_list:
                if len(assigned[seed]) >= max_dim:
                    continue
                for candidate in seed_sorted[seed]:
                    if candidate not in global_assigned:
                        assigned[seed].append(candidate)
                        global_assigned.add(candidate)
                        break
                results[seed][d] = np.array(assigned[seed].copy())

    for seed in seeds:
        results[seed][1] = [seed]

    return results

def generate_triplets(yeo_labels, n_triplets=100, first_network=1):
    indices_first = np.where(yeo_labels == first_network)[0]
    indices_other = np.where(yeo_labels != first_network)[0]

    if len(indices_first) < 2:
        raise ValueError("Not enough indices in the first network to sample two elements.")
    if len(indices_other) == 0:
        raise ValueError("No indices found outside the first network.")

    triplets = []
    for _ in range(n_triplets):
        chosen_first = np.random.choice(indices_first, size=2, replace=False)
        chosen_other = np.random.choice(indices_other, size=1)
        triplets.append([chosen_first[0], chosen_first[1], chosen_other[0]])

    return np.array(triplets)

def warn_if_constant_or_nan(matrix):
    matrix = np.asarray(matrix)  # Ensure it's a numpy array

    if np.isnan(matrix).any():
        warnings.warn("Matrix contains NaN values.", UserWarning)

    for idx, col in enumerate(matrix.T):  # Transpose and iterate over columns
        if np.all(col == col[0]):
            warnings.warn(f"Column {idx} is constant.", UserWarning)
