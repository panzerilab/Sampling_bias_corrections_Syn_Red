#!/usr/bin/env python3
from __future__ import print_function, division
import os
import numpy as np
import time
from joblib import Parallel, delayed
import scipy.io as sio
from scipy.io import loadmat
import pickle
import json

from Gaussian_bias_correction_methods_routines import (
    No_bias_correction_routine,
    informative_bias_correction_routine,
    no_info_bias_correction_routine,
    shuffle_subtr_bias_correction_routine,
    uniform_bias_correction_routine
    )
from utils_parallelization import (
    run_bias_corrections_on_data,
    run_task
    )


def main(input_data):
    """
    Main function to execute multiple bias-correction routines in parallel.
    
    Parameters:
        input_data (np.ndarray): The input dataset.
    
    Saves:
        MATLAB .mat files with computed bias-correction results.
    """
    start_time = time.time()  # Record the start time    
    T = 10  # Number of iterations for bias estimation (user-defined parameter)

    # Define the bias-correction routines and their corresponding titles
    bias_routines = [
        ("No_bias_correction", No_bias_correction_routine),
        ("informative_bias_correction", informative_bias_correction_routine),
        ("zero_info_bias_correction", no_info_bias_correction_routine),
        ("shuffle_subtr", shuffle_subtr_bias_correction_routine),
        ("Venka_bias_correction", uniform_bias_correction_routine)
    ]    

    # Build tasks for parallel processing
    tasks = []
    for bias_title, routine in bias_routines:
        tasks.append((input_data, bias_title, routine, T))

    # Run tasks in parallel using joblib
    results = Parallel(n_jobs=-1, verbose=10)(
        delayed(run_task)(*task) for task in tasks
    )

    # Save each result to a MATLAB-compatible .mat file
    for (bias_title, unbiased_results, error) in results:
        if error is not None:
            print(f"Error in job for bias_title={bias_title}): {error}")
            continue

        # Store results in a dictionary format suitable for MATLAB
        unbiased_results_ALL_mat = {bias_title: unbiased_results}
        file_name = f"pid_df_{bias_title}_DATA.mat"    

        # Save results as a .mat file
        sio.savemat(file_name, unbiased_results_ALL_mat)        

    # Print total runtime
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nTotal computed time: {total_time:.2f} seconds")




def execute_bias_correction_routines(input_data, bias_routines, num_iterations=10, parallel=False, save=False, save_format="python", cov_method='numpy'):
    """
    Executes multiple bias-correction routines in parallel based on the provided list of bias-correction names.

    Parameters:
        input_data (np.ndarray): The input dataset.
        bias_routines (list of str): A list of strings where each string is the name of a bias correction method
                                     to be applied. The available bias correction methods are predefined.
        num_iterations (int, optional): The number of iterations for bias estimation. Default is 10.
        save_format (str, optional): The format to save results. Can be 'MATLAB' for .mat files or 'Python' for a
                                     common Python format (e.g., .pkl or .json). Default is 'MATLAB'.

    Saves:
        The computed bias-correction results in the specified format.
    """
    # Define the bias-correction routines and their corresponding functions
    bias_routine_functions = {
        "no_bias_correction": No_bias_correction_routine,
        "informative_bias_correction": informative_bias_correction_routine,
        "zero_info_bias_correction": no_info_bias_correction_routine,
        "shuffsub_bias_correction": shuffle_subtr_bias_correction_routine,
        "Venkatesh_bias_correction": uniform_bias_correction_routine
    }

    # Check if all requested bias routines are valid
    for routine in bias_routines:
        if routine not in bias_routine_functions:
            raise ValueError(f"Invalid bias correction routine: {routine}")

    start_time = time.time()  # Record the start time

    # Build tasks for parallel processing
    tasks = []
    for bias_title in bias_routines:
        routine_function = bias_routine_functions[bias_title]
        tasks.append((input_data, bias_title, routine_function, num_iterations,cov_method))

    if parallel:
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(run_task)(*task) for task in tasks
        )
    else:
        results = [run_task(*task) for task in tasks]

    if save:
        # Create Results folder if it doesn't exist
        results_dir = "Results"
        os.makedirs(results_dir, exist_ok=True)

        # Save each result in the specified format
        for (bias_title, unbiased_results, error) in results:
            if error is not None:
                print(f"Error in job for bias_title={bias_title}): {error}")
                continue

            unbiased_results_ALL = {bias_title: unbiased_results}
            file_path = os.path.join(results_dir, f"pid_df_{bias_title}_DATA")

            if save_format == "matlab":
                sio.savemat(file_path + ".mat", unbiased_results_ALL)
            elif save_format == "python":
                with open(file_path + ".pkl", 'wb') as f:
                    pickle.dump(unbiased_results_ALL, f)
            elif save_format == "json":
                with open(file_path + ".json", 'w') as f:
                    json.dump(unbiased_results_ALL, f)
            else:
                raise ValueError("Unsupported save format. Choose 'matlab', 'python', or 'json'.")

    # Print total runtime
    end_time = time.time()
    total_time = end_time - start_time
    print(f"\nTotal computed time: {total_time:.2f} seconds")
    return results

if __name__ == '__main__':
    ###################
    # Load continuous data that should approximate a Gaussian distribution.
    # The input size should be (number of timepoints) × (number of dimensions of the total system).
    # Two source variables and two target variables should have the same dimensionality (dim) each.
    # For example, if dim = 4, the total number of dimensions should be 12.
    # input_data = loadmat('continuous_data.mat')

    # Example of generating synthetic multivariate Gaussian data.
    # Uncomment if you want to understand the structure of the input data.
    M=4 # Dimensionality of each variable
    seed = 22 # Random seed for reproducibility
    cov_each_var = np.eye(M)    # Identity covariance matrix for each variable
    cov = np.kron(np.eye(3), cov_each_var)    # Block-diagonal covariance matrix
    rng = np.random.default_rng(seed)      # Random number generator           
    z = rng.multivariate_normal(np.zeros(cov.shape[0]), cov, size=1800) # Generate Gaussian data
    input_data = z
    ###################
    
    # # # Example of generating synthetic multivariate Gaussian data from the 5 Venkatesh cases.
    # seed = 22
    # M=20
    # N = M * 3
    # dm = dx = dy = M
    # p=q=0.1
    # r=0
    # # Generate covariance matrix and data
    # cov, dm, dx, dy = gen_cov_matrix(N, M, p, q, r, 'both_unique')
    # rng = np.random.default_rng(seed)
    # input_data = rng.multivariate_normal(np.zeros(cov.shape[0]), cov, size=64)
    
    # Run the main function with the loaded data
    #results = main(input_data)
    results = execute_bias_correction_routines(input_data, ['no_bias_correction', 'shuffsub_bias_correction'], num_iterations=20, save_format="MATLAB",cov_method='lasso')