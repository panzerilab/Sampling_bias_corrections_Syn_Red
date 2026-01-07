from __future__ import print_function, division
import numpy as np
from scipy.io import loadmat
from scipy.io import savemat
from scipy.stats import poisson
import pandas as pd
import matplotlib.pyplot as plt
import math
from gpid.tilde_pid import exact_gauss_tilde_pid
from gpid.estimate import approx_pid_from_cov
from gpid.mmi_pid import mmi_pid
from scipy.linalg import det

data = loadmat('../BiasPaper2025/FigureR1A/data_4Stim.mat')

X1 = data['X1']
X2 = data['X2']
Y = data['Y'].flatten()

dm = 1
dx = 1
dy = 1

Bsweep = data['Bsweep'].flatten()
alphasweep = data['alphasweep'].flatten()
info_gauss = np.zeros((len(alphasweep), len(Bsweep)))

num_stimuli = 4
num_trials = X1.shape[2] // num_stimuli

for a_idx, alpha in enumerate(alphasweep):
    for b_idx, B in enumerate(Bsweep):

        X1_current = X1[a_idx, b_idx, :]
        X2_current = X2[a_idx, b_idx, :]

        mxy = np.vstack((Y, X1_current, X2_current))
        cov = np.corrcoef(mxy)
        cov_2 = np.cov(mxy)

        ret = exact_gauss_tilde_pid(cov, dm, dx, dy)
        ret2 = exact_gauss_tilde_pid(cov_2, dm, dx, dy)
        red_venk = ret[7]
        u1_venk = ret[5]
        u2_venk = ret[6]
        syn_venk = ret[8]
        info_tmp = u1_venk + u2_venk + red_venk + syn_venk
        info_gauss[a_idx, b_idx] = info_tmp
savemat('Sim_results/info_Gpid_4Stim.mat', {'info_FullGauss': info_gauss})
print('A')