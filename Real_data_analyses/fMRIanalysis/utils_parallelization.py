#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def run_bias_corrections_on_data(input_data, bias_title, routine, T, cov_method='numpy'):
    """
    Run a bias-correction routine on the input data.

    Parameters:
        input_data (np.ndarray): The input dataset.
        bias_title (str): Name of the bias-correction routine.
        routine (callable): The bias-correction function to use.
        T (int): Number of iterations for bias estimation.

    Returns:
        list: A list containing a dictionary with the bias-correction results.
    """
    seed = 22 # Set random seed for reproducibility
    n_samples = input_data.shape[0]
    
    # Run the specified bias-correction routine
    (imxy, uix, uiy, imx, imy, ri, si,
     HM, HXY, HMXY, imxy_ind, HM_ind, HXY_ind, HMXY_ind) = routine(input_data, T, n_samples, seed,cov_method)
    
    # Store the results in a structured format
    unbiased_results = []        
    unbiased_results.append({
    "imxy": imxy, "uix": uix, "uiy": uiy, "imx": imx, "imy": imy, "ri": ri, "si": si,
    "hm": HM, "hxy": HXY, "hmxy": HMXY, "imxy_ind": imxy_ind,
    "hm_ind": HM_ind, "hxy_ind": HXY_ind, "hmxy_ind": HMXY_ind
    })
    
    return (unbiased_results)          
              

def run_task(input_data, bias_title, routine, T,cov_method='numpy'):
    """
    Run a single bias-correction task while handling potential errors.

    Parameters:
        input_data (np.ndarray): The input dataset.
        bias_title (str): Name of the bias-correction routine.
        routine (callable): The bias-correction function.
        T (int): Number of iterations for bias estimation.

    Returns:
        tuple: (bias_title, result_dict, error) where error is None if successful.
    """
    try:
        # Run the bias-correction routine on the input data
        unbiased_results = run_bias_corrections_on_data(input_data, bias_title, routine, T,cov_method)
        return (bias_title, unbiased_results, None)
    except Exception as e:
        # Catch any errors and return them along with the bias title
        return (bias_title, None, e)