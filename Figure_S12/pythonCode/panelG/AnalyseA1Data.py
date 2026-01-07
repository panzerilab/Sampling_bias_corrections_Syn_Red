from __future__ import print_function, division
import numpy as np
from scipy.io import loadmat
import pandas as pd
import matplotlib.pyplot as plt
from gpid.tilde_pid import exact_gauss_tilde_pid
import scipy.io
import numpy as np
import pickle
from itertools import combinations
import os
print(os.getcwd())
rng = np.random.default_rng(seed = 22)
def process_pair(k, pair_indices, R, S):
    dm, dx, dy = 1, 1, 1
    neuron1 = pair_indices[k, 0]
    neuron2 = pair_indices[k, 1]

    R1 = R[neuron1]
    num_trials = R1.size
    R2 = R[neuron2]
    mxy = np.vstack((S, R1, R2))
    cov = np.corrcoef(mxy)

    retunbiased = exact_gauss_tilde_pid(cov, dm, dx, dy, verbose=False, ret_t_sigt=False,
                                        plot=False, unbiased=True, sample_size=num_trials)

    ret = exact_gauss_tilde_pid(cov, dm, dx, dy)
    red_venk = ret[7]
    u1_venk = ret[5]
    u2_venk = ret[6]
    syn_venk = ret[8]
    Ijoint = u1_venk + u2_venk + red_venk + syn_venk

    red_venksh = np.zeros(10)
    u1_venksh = np.zeros(10)
    u2_venksh = np.zeros(10)
    syn_venksh = np.zeros(10)
    Ijointsh = np.zeros(10)

    for i in range(10):
        z_shuffle = np.copy(S)
        rng.shuffle(z_shuffle, axis=1)
        mxysh = np.vstack((z_shuffle, R1, R2))
        covsh = np.corrcoef(mxysh)
        retsh = exact_gauss_tilde_pid(covsh, dm, dx, dy)
        red_venksh[i] = retsh[7]
        u1_venksh[i] = retsh[5]
        u2_venksh[i] = retsh[6]
        syn_venksh[i] = retsh[8]
        Ijointsh[i] = u1_venksh[i] + u2_venksh[i] + red_venksh[i] + syn_venksh[i]
    IjointCorr = Ijoint - np.mean(Ijointsh)
    red_venkCorr = red_venk - np.mean(red_venksh)
    u1_venkCorr = u1_venk - np.mean(u1_venksh)
    u2_venkCorr = u2_venk - np.mean(u2_venksh)
    syn_venkCorr = syn_venk - np.mean(syn_venksh)

    red_GpidCorr = retunbiased[7]
    u1_GpidCorr = retunbiased[5]
    u2_GpidCorr = retunbiased[6]
    syn_GpidCorr = retunbiased[8]
    Ijoint_GpidCorr = u1_GpidCorr + u2_GpidCorr + red_GpidCorr + syn_GpidCorr

    return IjointCorr, u1_venkCorr, u2_venkCorr, red_venkCorr, syn_venkCorr, Ijoint_GpidCorr, u1_GpidCorr, u2_GpidCorr, red_GpidCorr, syn_GpidCorr


def main():
     mat = scipy.io.loadmat('A1_python.mat')
     R_all = mat['R_all']
     S_all = mat['S_all']
     subfield_names = R_all.dtype.names

     all_Ijoint = []
     all_u1_venk = []
     all_u2_venk = []
     all_red_venk = []
     all_syn_venk = []

     all_Ijoint_gpid = []
     all_u1_gpid = []
     all_u2_gpid = []
     all_red_gpid = []
     all_syn_gpid = []
     for subset_name in subfield_names:
         R = R_all[subset_name][0, 0]
         S = S_all[subset_name][0, 0]

         N, T = R.shape
         pair_indices = np.array(list(combinations(range(N), 2)))
         for k in range(len(pair_indices)):
             result = process_pair(k, pair_indices, R, S)
             all_Ijoint.append(result[0])
             all_u1_venk.append(result[1])
             all_u2_venk.append(result[2])
             all_red_venk.append(result[3])
             all_syn_venk.append(result[4])

             all_Ijoint_gpid.append(result[5])
             all_u1_gpid.append(result[6])
             all_u2_gpid.append(result[7])
             all_red_gpid.append(result[8])
             all_syn_gpid.append(result[9])

     mean_Ijoint = np.mean(all_Ijoint)
     sem_Ijoint = np.std(all_Ijoint, ddof=1) / np.sqrt(len(all_Ijoint))

     mean_u1_venk = np.mean(all_u1_venk)
     sem_u1_venk = np.std(all_u1_venk, ddof=1) / np.sqrt(len(all_u1_venk))
     mean_u2_venk = np.mean(all_u2_venk)
     sem_u2_venk = np.std(all_u2_venk, ddof=1) / np.sqrt(len(all_u2_venk))

     mean_red_venk = np.mean(all_red_venk)
     sem_red_venk = np.std(all_red_venk, ddof=1) / np.sqrt(len(all_red_venk))

     mean_syn_venk = np.mean(all_syn_venk)
     sem_syn_venk = np.std(all_syn_venk, ddof=1) / np.sqrt(len(all_syn_venk))

     mean_Ijoint_gpid = np.mean(all_Ijoint_gpid)
     sem_Ijoint_gpid = np.std(all_Ijoint_gpid, ddof=1) / np.sqrt(len(all_Ijoint_gpid))

     mean_u1_gpid = np.mean(all_u1_gpid)
     sem_u1_gpid = np.std(all_u1_gpid, ddof=1) / np.sqrt(len(all_u1_gpid))
     mean_u2_gpid = np.mean(all_u2_gpid)
     sem_u2_gpid = np.std(all_u2_gpid, ddof=1) / np.sqrt(len(all_u2_gpid))

     mean_red_gpid = np.mean(all_red_gpid)
     sem_red_gpid = np.std(all_red_gpid, ddof=1) / np.sqrt(len(all_red_gpid))

     mean_syn_gpid = np.mean(all_syn_gpid)
     sem_syn_gpid = np.std(all_syn_gpid, ddof=1) / np.sqrt(len(all_syn_gpid))

     results = {
         'mean_Ijoint': mean_Ijoint,
         'sem_Ijoint': sem_Ijoint,
         'mean_u1_Corr': mean_u1_venk,
         'sem_u1_Corr': sem_u1_venk,
         'mean_u2_Corr': mean_u2_venk,
         'sem_u2_Corr': sem_u2_venk,
         'mean_red_Corr': mean_red_venk,
         'sem_red_Corr': sem_red_venk,
         'mean_syn_Corr': mean_syn_venk,
         'sem_syn_Corr': sem_syn_venk,
         'mean_Ijoint_gpid': mean_Ijoint_gpid,
         'sem_Ijoint_gpid': sem_Ijoint_gpid,
         'mean_u1_gpid': mean_u1_gpid,
         'sem_u1_gpid': sem_u1_gpid,
         'mean_u2_gpid': mean_u2_gpid,
         'sem_u2_gpid': sem_u2_gpid,
         'mean_red_gpid': mean_red_gpid,
         'sem_red_gpid': sem_red_gpid,
         'mean_syn_gpid': mean_syn_gpid,
         'sem_syn_gpid': sem_syn_gpid
     }
     df = pd.DataFrame([results])
     df.to_csv('resultsA1_shuffCorr.csv', index=False)
     #with open('resultsFrancis_biased.pkl', 'wb') as file:
     #    pickle.dump(results, file)


'''
    with open('resultsFrancis_unbiased.pkl', 'rb') as file:
        loaded_results = pickle.load(file)

    mean_Ijoint = loaded_results['mean_Ijoint']
    sem_Ijoint = loaded_results['sem_Ijoint']
    mean_u1_venk = loaded_results['mean_u1_venk']
    sem_u1_venk = loaded_results['sem_u1_venk']
    mean_u2_venk = loaded_results['mean_u2_venk']
    sem_u2_venk = loaded_results['sem_u2_venk']
    mean_red_venk = loaded_results['mean_red_venk']
    sem_red_venk = loaded_results['sem_red_venk']
    mean_syn_venk = loaded_results['mean_syn_venk']
    sem_syn_venk = loaded_results['sem_syn_venk']

    customColors_info_PID = [
        [0.9290, 0.6940, 0.1250],  # Color for 'Joint'
        [0, 0.4470, 0.7410],  # Color for 'Red'
        [0.4660, 0.6740, 0.1880],  # Color for 'Syn'
        [0.5, 0.5, 0.5],  # Color for 'U1'
        [0.4, 0.26, 0.13]  # Color for 'U2'
        ]

    means = [mean_Ijoint, mean_u1_venk, mean_u2_venk, mean_red_venk, mean_syn_venk]
    sems = [sem_Ijoint, sem_u1_venk, sem_u2_venk, sem_red_venk, sem_syn_venk]
    labels = ['Joint', 'Red', 'Syn', 'Unq1', 'Unq2']
    fig, ax = plt.subplots(figsize=(4, 4))
    bars = ax.bar(labels, means, yerr=sems, capsize=5, color=customColors_info_PID, edgecolor='black')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.savefig('RealData_results/Francis_Results_unbiased.svg', format='svg', bbox_inches='tight')
'''
if __name__ == '__main__':
    main()
