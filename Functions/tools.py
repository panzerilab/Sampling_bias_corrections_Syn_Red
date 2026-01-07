#!/usr/bin/env python3
from __future__ import print_function, division
import warnings
import pickle
import numpy as np
import pandas as pd
from scipy.special import psi
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


def random_rotation_mxy(cov, dm, dx, dy, seed=None):
    """
    Perform random rotations on the M, X and Y components of a covariance matrix.
    Returns a copy.
    """
    cov = cov.copy()

    # rng = np.random.default_rng(seed)

    # R = ortho_group.rvs(dm, seed=rng)
    # cov[:dm, :] = R @ cov[:dm, :]
    # cov[:, :dm] = cov[:, :dm] @ R.T
    # R = ortho_group.rvs(dx, seed=rng)
    # cov[dm:dm+dx, :] = R @ cov[dm:dm+dx, :]
    # cov[:, dm:dm+dx] = cov[:, dm:dm+dx] @ R.T
    # R = ortho_group.rvs(dy, seed=rng)
    # cov[dm+dx:, :] = R @ cov[dm+dx:, :]
    # cov[:, dm+dx:] = cov[:, dm+dx:] @ R.T

    rng = np.random.default_rng(seed=22)

    R = ortho_group.rvs(dm, random_state=rng)
    cov[:dm, :] = R @ cov[:dm, :]
    cov[:, :dm] = cov[:, :dm] @ R.T
    R = ortho_group.rvs(dx, random_state=rng)
    cov[dm:dm + dx, :] = R @ cov[dm:dm + dx, :]
    cov[:, dm:dm + dx] = cov[:, dm:dm + dx] @ R.T
    R = ortho_group.rvs(dy, random_state=rng)
    cov[dm + dx:, :] = R @ cov[dm + dx:, :]
    cov[:, dm + dx:] = cov[:, dm + dx:] @ R.T

    return cov


def move_cov_block_to_end(cov, i, j):
    """
    Moves the variables indexed by i:j to the end of the covariance matrix.
    Assumes the covariance matrix is symmetric.
    Negative indices will not work.
    """
    if i >= j:
        raise ValueError('Index i must be less than index j')
    cov_new = np.delete(cov, np.s_[i:j], axis=0)
    cov_new = np.delete(cov_new, np.s_[i:j], axis=1)
    diag_block = cov[i:j, i:j]
    nondiag_block = np.delete(cov[i:j, :], np.s_[i:j], axis=1)
    return np.block([[cov_new, nondiag_block.T], [nondiag_block, diag_block]])


def swap_x_and_y(cov, dm, dx, dy):
    return move_cov_block_to_end(cov, dm, dm + dx), dm, dy, dx


def merge_covs(cov1, cov2, dm1, dx1, dy1, dm2=None, dx2=None, dy2=None,
               random_rotn=False, seed=None):
    """
    Merges two covariance matrices of given M, X and Y dimensions.
    Performs a random rotation on M, X and Y after merging, if random_rotn is True.
    """
    if dm2 is None:
        dm2 = dm1
    if dx2 is None:
        dx2 = dx1
    if dy2 is None:
        dy2 = dy1

    dm = dm1 + dm2
    dx = dx1 + dx2
    dy = dy1 + dy2

    zero_block = np.zeros((cov1.shape[0], cov2.shape[0]))
    cov = np.block([[cov1, zero_block], [zero_block.T, cov2]])

    # Move required blocks to end: order of ops is important!
    cov = move_cov_block_to_end(cov, dm1, dm1 + dx1 + dy1)
    # Diag now reads: M1, M2, X2, Y2, X1, Y1
    cov = move_cov_block_to_end(cov, dm, dm + dx2)
    # Diag now reads: M1, M2, Y2, X1, Y1, X2
    cov = move_cov_block_to_end(cov, dm + dy2 + dx1, dm + dy2 + dx1 + dy1)
    # Diag now reads: M1, M2, Y2, X1, X2, Y1
    cov = move_cov_block_to_end(cov, dm, dm + dy2)
    # Diag now reads: M1, M2, X1, X2, Y1, Y2

    if random_rotn:
        cov = random_rotation_mxy(cov, dm, dx, dy, seed=seed)

    return cov, dm, dx, dy


def gen_adjacency_matrix(N, M, p, q, r, mode):
    rng = np.random.default_rng(seed=22)
    A = np.zeros((N, N))
    if mode == 'both_unique':
        B = np.array([[1, 0, 2, 0],
                      [0, 1, 0, 2],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
        subnet_sizes = np.array([M // 2, M // 2, M, M])
    elif mode == 'fully_redundant':
        B = np.array([[1, 2, 2],
                      [0, 1, 0],
                      [0, 0, 1]])
        subnet_sizes = np.array([M, M, M])
    elif mode == 'unique_plus_redundant':
        B = np.array([[1, 0, 0, 2, 0],
                      [0, 1, 0, 0, 2],
                      [0, 0, 1, 2, 2],
                      [0, 0, 0, 1, 0],
                      [0, 0, 0, 0, 1]])
        subnet_sizes = np.array([3, 3, M - 6, M, M])
    elif mode == 'zero_synergy':
        B = np.array([[1, 2, 0],
                      [0, 1, 2],
                      [0, 0, 1]])
        subnet_sizes = np.array([M, M, M])
    elif mode == 'high_synergy':
        B = np.array([[1, 2, 0],
                      [0, 1, 0],
                      [0, 0, 1]])
        subnet_sizes = np.array([M, M, M])
    else:
        raise ValueError('Unrecognized mode %s' % mode)
    assert subnet_sizes.sum() == N
    L = subnet_sizes.size
    for i in range(L):
        for j in range(L):
            subnet_margins = np.r_[0, np.cumsum(subnet_sizes)]
            Mi = subnet_margins[i]
            Ki = subnet_margins[i + 1]
            Mj = subnet_margins[j]
            Kj = subnet_margins[j + 1]
            m = subnet_sizes[i]
            k = subnet_sizes[j]
            if B[i, j] == 1:
                A[Mi: Ki, Mj: Kj] = (rng.uniform(size=(m, k)) <= p).astype(int)
            elif B[i, j] == 2:
                A[Mi: Ki, Mj: Kj] = (rng.uniform(size=(m, k)) <= q).astype(int)
            else:
                A[Mi: Ki, Mj: Kj] = (rng.uniform(size=(m, k)) <= r).astype(int)
    if mode == 'fully_redundant':
        A[:M, 2 * M:] = A[:M, M:2 * M]
    elif mode == 'unique_plus_redundant':
        A[M // 2:M, 2 * M:] = A[M // 2:M, M:2 * M]
    np.fill_diagonal(A, 0)
    return A, B


def gauss_weighted_adj(A):
    """Adds random gaussian weights to an adjacency matrix"""
    rng = np.random.default_rng(seed=22)
    return A * rng.normal(0, 1, size=A.shape)


def gen_cov_matrix(N, M, p, q, r, mode):
    if mode == 'bit_of_all':
        cov1, dm1, dx1, dy1 = gen_cov_matrix(N // 2, M // 2, p, q, r, mode='zero_synergy')
        cov1, dm1, dx1, dy1 = swap_x_and_y(cov1, dm1, dx1, dy1)
        cov2, dm2, dx2, dy2 = gen_cov_matrix(N // 2, M // 2, p, q, r, mode='high_synergy')
        return merge_covs(cov1, cov2, dm1, dx1, dy1, dm2, dx2, dy2, random_rotn=True)
    A, _ = gen_adjacency_matrix(N, M, p, q, r, mode)
    A = gauss_weighted_adj(A)
    dm, dx, dy = M, M, M
    if mode in ['both_unique', 'fully_redundant', 'high_synergy']:
        sigm = np.eye(dm)
        sigx_m = np.eye(dx)
        sigy_m = np.eye(dy)
        hx = A[:dm, dm:dm + dx].T
        if mode == 'fully_redundant':
            hy = hx
            sigw = 0.9 * np.eye(dx)  # Assumes dx == dy
        elif mode == 'high_synergy':
            hy = np.zeros((dy, dm))
            sigw = 0.8 * np.eye(dx)  # Assumes dx == dy
        else:
            hy = A[:dm, dm + dx:].T
            sigw = np.zeros((dx, dy))
        cov = np.block([[sigm, sigm @ hx.T, sigm @ hy.T],
                        [hx @ sigm, hx @ sigm @ hx.T + sigx_m, hx @ sigm @ hy.T + sigw],
                        [hy @ sigm, hy @ sigm @ hx.T + sigw.T, hy @ sigm @ hy.T + sigy_m]])
    elif mode == 'zero_synergy':
        hx = A[:dm, dm:dm + dx].T
        hyx = A[dm:dm + dx, dm + dx:].T
        hy = hyx @ hx
        sigm = np.eye(dm)
        sigx_m = np.eye(dx)
        covx = hx @ sigm @ hx.T + sigx_m
        sigy_x = np.eye(dx)
        cov = np.block([[sigm, sigm @ hx.T, sigm @ hy.T],
                        [hx @ sigm, covx, covx @ hyx.T],
                        [hy @ sigm, hyx @ covx, hyx @ covx @ hyx.T + sigy_x]])
    else:
        raise ValueError('Unrecognized mode %s' % mode)
    return cov, dm, dx, dy

def lower_info(cov_matrix,a):
    np.fill_diagonal(cov_matrix, a*np.diagonal(cov_matrix))
    return cov_matrix

def inf_from_cov(cov_hat):
    M = cov_hat.shape[1] // 3
    dm = dx = dy = M
    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HX = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:dm + dx, dm:dm + dx])[1] + (dx) * (1 + np.log(2 * np.pi)))
    HY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm + dx:, dm + dx:])[1] + (dy) * (1 + np.log(2 * np.pi)))
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    HMX = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm + dx, :dm + dx])[1] + (dm + dx) * (1 + np.log(2 * np.pi)))

    # Combine them if you want the joint covariance matrix of M and Y
    Sigma_MY_combined = np.block([
        [cov_hat[:dm, :dm], cov_hat[:dm, dm + dx:]],
        [cov_hat[dm + dx:, :dm], cov_hat[dm + dx:, dm + dx:]]
    ])
    HMY = 0.5 / np.log(2) * (la.slogdet(Sigma_MY_combined)[1] + (dm + dy) * (1 + np.log(2 * np.pi)))
    # Compute biased PID estimates
    try:
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        # imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
        imx, imy, imxy, union_info = ret[:4]
    except:  # Mainly to catch LinAlgWarning's
        # imxy, uix, uiy, ri, si = [np.nan, ] * 5
        imx, imy, imxy, union_info = [np.nan, ] * 4


    imx = HM + HX - HMX
    imy = HM + HY - HMY
    imxy = HM + HXY - HMXY
    uix = union_info - imy
    uiy = union_info - imx
    ri = imx + imy - union_info
    si = imxy - union_info
    return (imxy, imx, imy, ri, si, uix, uiy, union_info, HM, HXY, HMXY)

def No_bias_correction_routine(input_data, T, ntrials, seed):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M
    cov_hat = np.cov(input_data.T)

    imxy, imx, imy, ri, si, uix, uiy, union_info, HM, HXY, HMXY = inf_from_cov(cov_hat)

    cov_hat_ind = np.copy(cov_hat)
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

    # Venkatesh bias correction approach
    imxy_ind, imx_ind, imy_ind, ri_ind, si_ind, uix_ind, uiy_ind, union_info_ind, HM_ind, HXY_ind, HMXY_ind = inf_from_cov(cov_hat_ind)

    return imxy, si, ri, uix, uiy, union_info, cov_hat_ind, imxy_ind, imx, imy, HM, HXY, HMXY, HM_ind, HXY_ind, HMXY_ind


def informative_bias_correction_routine(input_data, T, ntrials, seed):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M

    cov_hat = np.cov(input_data.T)
    imxy, imx, imy, ri, si, uix, uiy, union_info, HM, HXY, HMXY = inf_from_cov(cov_hat)

    ks = np.arange(1, M + 1)
    ks2 = np.arange(1, 2 * M + 1)
    ks3 = np.arange(1, 3 * M + 1)

    # Compute bias terms using vectorized psi
    biasH_M_goodman = M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks) / 2))
    biasH_XY_goodman = 2 * M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks2) / 2))
    biasH_MXY_goodman = 3 * M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks3) / 2))

    # Compute analytical bias terms
    log2 = np.log(2)
    Jointinfo_analytical_bias = (1 / (2 * log2)) * (biasH_M_goodman + biasH_XY_goodman - biasH_MXY_goodman)
    Info_analytical_bias = (1 / (2 * log2)) * (2 * biasH_M_goodman - biasH_XY_goodman)

    cov_hat_ind = np.copy(cov_hat)
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))


    # Now estimate the bias using T re-samplings from the estimated covariance
    rng = np.random.default_rng(seed)
    imxy_rnd_all = np.zeros(T)
    imx_rnd_all = np.zeros(T)
    imy_rnd_all = np.zeros(T)
    uix_rnd_all= np.zeros(T)
    uiy_rnd_all = np.zeros(T)
    ri_rnd_all = np.zeros(T)
    si_rnd_all = np.zeros(T)
    union_rnd_all = np.zeros(T)
    HM_rnd_all = np.zeros(T)
    HXY_rnd_all = np.zeros(T)
    HMXY_rnd_all = np.zeros(T)
    imxy_ind_rnd_all = np.zeros(T)
    HM_ind_rnd_all = np.zeros(T)
    HXY_ind_rnd_all = np.zeros(T)
    HMXY_ind_rnd_all = np.zeros(T)

    for i in range(T):
        z = rng.multivariate_normal(np.zeros(cov_hat.shape[0]), cov_hat, size=ntrials)
        cov_hat_rnd = np.cov(z.T)

        imxy_rnd, imx_rnd, imy_rnd, ri_rnd, si_rnd, uix_rnd, uiy_rnd, union_info_rnd, HM_rnd, HXY_rnd, HMXY_rnd = inf_from_cov(
            cov_hat_rnd)

        cov_rnd_ind = np.copy(cov_hat_rnd)
        sig_m = cov_hat_rnd[:dm, :dm]
        sig_xy = cov_hat_rnd[dm:, dm:]
        sig_xy_m = cov_hat_rnd[dm:, :dm]
        sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
        cov_rnd_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
        cov_rnd_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

        ret_ind_rnd = exact_gauss_tilde_pid(cov_rnd_ind, dm, dx, dy)
        imxy_ind_rnd, uix_ind_rnd, uiy_ind_rnd, ri_ind_rnd, si_ind_rnd = ret_ind_rnd[2], *ret_ind_rnd[-4:]
        HM_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_rnd_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
        HXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_rnd_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
        HMXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_rnd_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

        imxy_rnd_all[i] = imxy_rnd
        uix_rnd_all[i] = uix_rnd
        uiy_rnd_all[i] = uiy_rnd
        imx_rnd_all[i] = imx_rnd
        imy_rnd_all[i] = imy_rnd
        union_rnd_all[i] = union_info_rnd
        ri_rnd_all[i] = ri_rnd
        si_rnd_all[i] = si_rnd
        HM_rnd_all[i] = HM_rnd
        HXY_rnd_all[i] = HXY_rnd
        HMXY_rnd_all[i] = HMXY_rnd
        imxy_ind_rnd_all[i] = imxy_ind_rnd
        HM_ind_rnd_all[i] = HM_ind_rnd
        HXY_ind_rnd_all[i] = HXY_ind_rnd
        HMXY_ind_rnd_all[i] = HMXY_ind_rnd

        # Remove the bias
    imxy =imxy - Jointinfo_analytical_bias # 2*imxy - np.mean(imxy_rnd_all) #
    imx = imx - Info_analytical_bias #2*imx - np.mean(imx_rnd_all) #
    imy = imy - Info_analytical_bias #2*imy - np.mean(imy_rnd_all) #
    union = 2*union_info - np.mean(union_rnd_all)
    uix =  union  - imy
    uiy =  union  - imy
    ri = imx + imy - union  # ri - np.mean(ri_rnd_bias)
    si = imxy - union  # si - np.mean(si_rnd_bias)
    HM = 2*HM - np.mean(HM_rnd_all)
    HXY = 2*HXY - np.mean(HXY_rnd_all)
    HMXY = 2*HMXY - np.mean(HMXY_rnd_all)
    imxy_ind = 2*imxy_ind - np.mean(imxy_ind_rnd_all)
    HM_ind = 2*HM_ind - np.mean(HM_ind_rnd_all)
    HXY_ind = 2*HXY_ind - np.mean(HXY_ind_rnd_all)
    HMXY_ind = 2*HMXY_ind - np.mean(HMXY_ind_rnd_all)

    return imxy, si, ri, uix, uiy, union, cov_hat_ind, imxy_ind, imx, imy, HM, HXY, HMXY, HM_ind, HXY_ind, HMXY_ind

def shuffle_subtr_bias_correction_routine(input_data, T, ntrials, seed):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M

    cov_hat = np.cov(input_data.T)
    imxy, imx, imy, ri, si, uix, uiy, union_info, HM, HXY, HMXY = inf_from_cov(cov_hat)

    ks = np.arange(1, M + 1)
    ks2 = np.arange(1, 2 * M + 1)
    ks3 = np.arange(1, 3 * M + 1)

    # Compute bias terms using vectorized psi
    biasH_M_goodman = M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks) / 2))
    biasH_XY_goodman = 2 * M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks2) / 2))
    biasH_MXY_goodman = 3 * M * np.log(2 / (ntrials - 1)) + np.sum(psi((ntrials - ks3) / 2))

    # Compute analytical bias terms
    log2 = np.log(2)
    Jointinfo_analytical_bias = (1 / (2 * log2)) * (biasH_M_goodman + biasH_XY_goodman - biasH_MXY_goodman)
    Info_analytical_bias = (1 / (2 * log2)) * (2 * biasH_M_goodman - biasH_XY_goodman)

    cov_hat_ind = np.copy(cov_hat)
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

    # Now estimate the bias using T re-samplings from the estimated covariance
    rng = np.random.default_rng(seed)
    imxy_rnd_bias = np.zeros(T)
    imx_rnd_bias  = np.zeros(T)
    imy_rnd_bias  = np.zeros(T)
    uix_rnd_bias = np.zeros(T)
    uiy_rnd_bias = np.zeros(T)
    ri_rnd_bias = np.zeros(T)
    si_rnd_bias = np.zeros(T)
    union_rnd_bias = np.zeros(T)
    HM_rnd_bias = np.zeros(T)
    HXY_rnd_bias = np.zeros(T)
    HMXY_rnd_bias = np.zeros(T)
    imxy_ind_rnd_bias = np.zeros(T)
    HM_ind_rnd_bias = np.zeros(T)
    HXY_ind_rnd_bias = np.zeros(T)
    HMXY_ind_rnd_bias = np.zeros(T)

    for i in range(T):
        ######## Shuffle subtraction approach
        # Ground-truth PID terms from the current covariance estimate
        input_data_shuffled = input_data.copy()
        perm = rng.permutation(input_data.shape[0])
        input_data_shuffled[:, :dm] = input_data[perm, :dm]
        cov_shuff = np.cov(input_data_shuffled.T)
        ##########

        imxy_shuff, imx_shuff, imy_shuff, ri_shuff, si_shuff, uix_shuff, uiy_shuff, union_info_shuff, HM_shuff, HXY_shuff, HMXY_shuff = inf_from_cov(cov_shuff)
        # Use the conditionally-independent condition \Sigma_{X,Y|M} in the covariance matrix "cov"
        cov_ind_rnd = np.copy(cov_shuff)
        sig_m = cov_shuff[:dm, :dm]
        sig_xy = cov_shuff[dm:, dm:]
        sig_xy_m = cov_shuff[dm:, :dm]
        sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
        cov_ind_rnd[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
        cov_ind_rnd[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

        ret_ind_rnd = exact_gauss_tilde_pid(cov_ind_rnd, dm, dx, dy)
        imxy_ind_shuff, uix_ind_shuff, uiy_ind_shuff, ri_ind_shuff, si_ind_shuff = ret_ind_rnd[2], *ret_ind_rnd[-4:]
        # Compute the GT Entropy terms in the IND case using the covariance matrix "cov_ind"
        HM_ind_shuff = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
        HXY_ind_shuff = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
        HMXY_ind_shuff = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

        imxy_rnd_bias[i] = imxy_shuff
        uix_rnd_bias[i] = uix_shuff
        uiy_rnd_bias[i] = uiy_shuff
        imx_rnd_bias[i] = imx_shuff
        imy_rnd_bias[i] = imy_shuff
        ri_rnd_bias[i] = ri_shuff
        si_rnd_bias[i] = si_shuff
        union_rnd_bias[i] = union_info_shuff
        HM_rnd_bias[i] = HM_shuff
        HXY_rnd_bias[i] = HXY_shuff
        HMXY_rnd_bias[i] = HMXY_shuff
        imxy_ind_rnd_bias[i] = imxy_ind_shuff
        HM_ind_rnd_bias[i] = HM_ind_shuff
        HXY_ind_rnd_bias[i] = HXY_ind_shuff
        HMXY_ind_rnd_bias[i] = HMXY_ind_shuff

        # Remove the bias
    imxy = imxy - Jointinfo_analytical_bias # imxy - np.mean(imxy_rnd_bias) #
    imx = imx - Info_analytical_bias #imx - np.mean(imx_rnd_bias) #
    imy = imy - Info_analytical_bias #imy - np.mean(imy_rnd_bias) #
    union = union_info - np.mean(union_rnd_bias)
    uix = union - imy
    uiy = union - imx  #
    ri = imx+imy - union #ri - np.mean(ri_rnd_bias)
    si = imxy - union #si - np.mean(si_rnd_bias)
    HM = HM - np.mean(HM_rnd_bias)
    HXY = HXY - np.mean(HXY_rnd_bias)
    HMXY = HMXY - np.mean(HMXY_rnd_bias)
    imxy_ind = imxy_ind - np.mean(imxy_ind_rnd_bias)
    HM_ind = HM_ind - np.mean(HM_ind_rnd_bias)
    HXY_ind = HXY_ind - np.mean(HXY_ind_rnd_bias)
    HMXY_ind = HMXY_ind - np.mean(HMXY_ind_rnd_bias)

    return imxy, si, ri, uix, uiy, union, cov_hat_ind, imxy_ind, imx, imy, HM, HXY, HMXY, HM_ind, HXY_ind, HMXY_ind


def uniform_bias_correction_routine(input_data, T, ntrials, seed):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M
    cov_hat = np.cov(input_data.T)

    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    # Compute biased PID estimates
    try:
        # Venkatesh bias correction approach
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy, sample_size=ntrials, unbiased=True)

        # ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
        union = ret[3]
    except:  # Mainly to catch LinAlgWarning's
        imxy, uix, uiy, ri, si = [np.nan, ] * 5
        union = np.nan

    imx = uix + ri
    imy = uiy + ri

    cov_hat_ind = np.copy(cov_hat)
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

    # Venkatesh bias correction approach
    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy, sample_size=ntrials, unbiased=True)

    # ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

    return imxy, si, ri, uix, uiy, union, cov_hat_ind, imxy_ind, imx, imy, HM, HXY, HMXY, HM_ind, HXY_ind, HMXY_ind



