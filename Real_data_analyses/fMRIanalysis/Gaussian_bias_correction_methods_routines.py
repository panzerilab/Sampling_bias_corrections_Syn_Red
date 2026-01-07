#!/usr/bin/env python3
from __future__ import print_function, division
import numpy as np
import numpy.linalg as la
from gpid.tilde_pid import exact_gauss_tilde_pid
from gpid.utils import solve
from scipy.special import psi
from sklearn.covariance import LedoitWolf
from sklearn.covariance import GraphicalLassoCV
from sklearn.covariance import oas
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer


def inf_from_cov(cov_hat):
    M = cov_hat.shape[1] // 3
    dm = dx = dy = M

    Sigma_MY_combined = np.block([
        [cov_hat[:dm, :dm], cov_hat[:dm, dm + dx:]],
        [cov_hat[dm + dx:, :dm], cov_hat[dm + dx:, dm + dx:]]
    ])

    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HX = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:dm + dx, dm:dm + dx])[1] + (dx) * (1 + np.log(2 * np.pi)))
    HY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm + dx:, dm + dx:])[1] + (dy) * (1 + np.log(2 * np.pi)))
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    HMX = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm + dx, :dm + dx])[1] + (dm + dx) * (1 + np.log(2 * np.pi)))
    HMY = 0.5 / np.log(2) * (la.slogdet(Sigma_MY_combined)[1] + (dm + dy) * (1 + np.log(2 * np.pi)))
    # Compute biased PID estimates
    try:
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        # imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
        imx, imy, imxy, union_info = ret[:4]
    except:  # Mainly to catch LinAlgWarning's
        # imxy, uix, uiy, ri, si = [np.nan, ] * 5
        imx, imy, imxy, union_info = [np.nan, ] * 4


    # imx = HM + HX - HMX
    # imy = HM + HY - HMY
    # imxy = HM + HXY - HMXY
    uix = union_info - imy
    uiy = union_info - imx
    uix_uiy = uix + uiy
    ri = imx + imy - union_info
    si = imxy - union_info
    return (imxy, imx, imy, ri, si, uix_uiy, union_info, HM, HXY, HMXY)

def informative_bias_correction_routine(input_data, T, ntrials, seed, cov_method='numpy', precomp_cov=None):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M

    if precomp_cov is None:
        if cov_method == 'numpy':
            cov_hat = np.cov(input_data.T)
        elif cov_method == 'shrink':
            lw = LedoitWolf()
            lw.fit(input_data)  # A: rows are samples, columns are features
            cov_hat = lw.covariance_
        elif cov_method =='lasso':
            imputer = SimpleImputer(strategy='mean')
            X_clean = imputer.fit_transform(input_data)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_clean)
            model = GraphicalLassoCV()
            model.fit(X_scaled)
            cov_hat = model.covariance_
        elif cov_method == 'osa':
            cov_hat, _ = oas(input_data)
        else:
            cov_hat=precomp_cov

    imxy, imx, imy, ri, si, uix_uiy, union_info, HM, HXY, HMXY = inf_from_cov(cov_hat)

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
    uix_uiy_rnd_all = np.zeros(T)
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
        # if cov_method == 'numpy':
        #     cov_hat_rnd = np.cov(z.T)
        # elif cov_method == 'shrink':
        #     lw = LedoitWolf()
        #     lw.fit(z)  # A: rows are samples, columns are features
        #     cov_hat_rnd = lw.covariance_
        # elif cov_method == 'lasso':
        #     imputer = SimpleImputer(strategy='mean')
        #     X_clean = imputer.fit_transform(z)
        #     scaler = StandardScaler()
        #     X_scaled = scaler.fit_transform(X_clean)
        #     model = GraphicalLassoCV()
        #     model.fit(X_scaled)
        #     cov_hat_rnd = model.covariance_
        # elif cov_method == 'osa':
        #     cov_hat_rnd, _ = oas(z)


        imxy_rnd, imx_rnd, imy_rnd, ri_rnd, si_rnd, uix_uiy_rnd, union_info_rnd, HM_rnd, HXY_rnd, HMXY_rnd = inf_from_cov(
            cov_hat_rnd)

        cov_hat_ind = np.copy(cov_hat_rnd)
        sig_m = cov_hat_rnd[:dm, :dm]
        sig_xy = cov_hat_rnd[dm:, dm:]
        sig_xy_m = cov_hat_rnd[dm:, :dm]
        sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
        cov_hat_ind[dm:dm + dx, dm + dx:] = sig_xy[:dx, dx:]
        cov_hat_ind[dm + dx:, dm:dm + dx] = sig_xy[dx:, :dx]

        ret_ind_rnd = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
        imxy_ind_rnd, uix_ind_rnd, uiy_ind_rnd, ri_ind_rnd, si_ind_rnd = ret_ind_rnd[2], *ret_ind_rnd[-4:]
        HM_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
        HXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
        HMXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))

        imxy_rnd_all[i] = imxy_rnd
        uix_uiy_rnd_all[i] = uix_uiy_rnd
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
    imxy_corr = imxy - Jointinfo_analytical_bias # 2*imxy - np.mean(imxy_rnd_all) #
    imx_corr  = imx - Info_analytical_bias # 2*imx - np.mean(imx_rnd_all) #
    imy_corr  = imy - Info_analytical_bias # 2*imy - np.mean(imy_rnd_all) #
    union_corr = 2*union_info - np.mean(union_rnd_all)
    uix_corr = union_corr - imy_corr  # uix_uiy - np.mean(uix_uiy_rnd_bias)
    uiy_corr = union_corr - imx_corr  # uix_uiy - np.mean(uix_uiy_rnd_bias)
    ri_corr = imx_corr + imy_corr - union_corr  # ri - np.mean(ri_rnd_bias)
    si_corr = imxy_corr - union_corr  # si - np.mean(si_rnd_bias)
    HM = 2*HM - np.mean(HM_rnd_all)
    HXY = 2*HXY - np.mean(HXY_rnd_all)
    HMXY = 2*HMXY - np.mean(HMXY_rnd_all)
    imxy_ind = 2*imxy_ind - np.mean(imxy_ind_rnd_all)
    HM_ind = 2*HM_ind - np.mean(HM_ind_rnd_all)
    HXY_ind = 2*HXY_ind - np.mean(HXY_ind_rnd_all)
    HMXY_ind = 2*HMXY_ind - np.mean(HMXY_ind_rnd_all)

    return imxy_corr, uix_corr, uiy_corr, imx_corr, imy_corr, ri_corr, si_corr, HM, HXY, HMXY, imxy_ind, HM_ind, HXY_ind, HMXY_ind

def no_info_bias_correction_routine(input_data, T, ntrials, seed):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M

    cov_hat = np.cov(input_data.T)
    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )
    # Compute biased PID estimates
    try:
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
    except:  # Mainly to catch LinAlgWarning's
        imxy, uix, uiy, ri, si = [np.nan,] * 5

    # uix_uiy = uix + uiy
    imx = uix + ri
    imy = uiy + ri

    cov_hat_ind = np.copy(cov_hat)
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm+dx, dm+dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm+dx:, dm:dm+dx] = sig_xy[dx:, :dx]

    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )

    # Now estimate the bias using T re-samplings from the estimated covariance
    rng = np.random.default_rng(seed)
    imxy_rnd_bias = np.zeros(T)
    uix_rnd_bias = np.zeros(T)
    uiy_rnd_bias = np.zeros(T)
    imx_rnd_bias = np.zeros(T)
    imy_rnd_bias = np.zeros(T)
    ri_rnd_bias = np.zeros(T)
    si_rnd_bias = np.zeros(T)
    HM_rnd_bias = np.zeros(T)
    HXY_rnd_bias = np.zeros(T)
    HMXY_rnd_bias = np.zeros(T)
    imxy_ind_rnd_bias = np.zeros(T)
    HM_ind_rnd_bias = np.zeros(T)
    HXY_ind_rnd_bias = np.zeros(T)
    HMXY_ind_rnd_bias = np.zeros(T)

    for i in range(T):

        ######## Zero - info bias correction
        # Ground-truth PID terms from the current covariance estimate
        mean = np.zeros(M)
        cov_each_var = np.eye(M)
        cov = np.kron(np.eye(3), cov_each_var)
        ########

        HM_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
        HXY_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
        HMXY_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )
        ret = exact_gauss_tilde_pid(cov, dm, dx, dy)
        imxy_GT_rnd, uix_GT_rnd, uiy_GT_rnd, ri_GT_rnd, si_GT_rnd = ret[2], *ret[-4:]
        # uix_uiy_GT_rnd = uix_GT_rnd + uiy_GT_rnd
        imx_GT_rnd = uix_GT_rnd + ri_GT_rnd
        imy_GT_rnd = uiy_GT_rnd + ri_GT_rnd

        # Use the conditionally-independent condition \Sigma_{X,Y|M} in the covariance matrix "cov"
        cov_ind_rnd = np.copy(cov)
        sig_m = cov[:dm, :dm]
        sig_xy = cov[dm:, dm:]
        sig_xy_m = cov[dm:, :dm]
        sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
        cov_ind_rnd[dm:dm+dx, dm+dx:] = sig_xy[:dx, dx:]
        cov_ind_rnd[dm+dx:, dm:dm+dx] = sig_xy[dx:, :dx]

        # Compute the GT JointI using after whitening the covariance matrix "cov_ind"
        ret_ind_rnd = exact_gauss_tilde_pid(cov_ind_rnd, dm, dx, dy)
        imxy_ind_GT_rnd, uix_ind_rnd, uiy_ind_rnd, ri_ind_rnd, si_ind_rnd = ret_ind_rnd[2], *ret_ind_rnd[-4:]
        # Compute the GT Entropy terms in the IND case using the covariance matrix "cov_ind"
        HM_ind_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
        HXY_ind_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
        HMXY_ind_GT_rnd = 0.5 / np.log(2) * (la.slogdet(cov_ind_rnd)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )

        # informative bias correction approach
        z = rng.multivariate_normal(np.zeros(cov.shape[0]), cov, size=ntrials)
        cov_hat_rnd = np.cov(z.T)

        HM_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_rnd[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
        HXY_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_rnd[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
        HMXY_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_rnd)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )
        # Compute biased PID estimates
        try:
            ret = exact_gauss_tilde_pid(cov_hat_rnd, dm, dx, dy)
            imxy_rnd, uix_rnd, uiy_rnd, ri_rnd, si_rnd = ret[2], *ret[-4:]
        except:  # Mainly to catch LinAlgWarning's
            imxy_rnd, uix_rnd, uiy_rnd, ri_rnd, si_rnd = [np.nan,] * 5

        # uix_uiy_rnd = uix_rnd + uiy_rnd
        imx_rnd = uix_rnd + ri_rnd
        imy_rnd = uiy_rnd + ri_rnd

        cov_hat_ind = np.copy(cov_hat_rnd)
        sig_m = cov_hat_rnd[:dm, :dm]
        sig_xy = cov_hat_rnd[dm:, dm:]
        sig_xy_m = cov_hat_rnd[dm:, :dm]
        sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
        cov_hat_ind[dm:dm+dx, dm+dx:] = sig_xy[:dx, dx:]
        cov_hat_ind[dm+dx:, dm:dm+dx] = sig_xy[dx:, :dx]

        ret_ind_rnd = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
        imxy_ind_rnd, uix_ind_rnd, uiy_ind_rnd, ri_ind_rnd, si_ind_rnd = ret_ind_rnd[2], *ret_ind_rnd[-4:]
        HM_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) )
        HXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
        HMXY_ind_rnd = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )

        imxy_rnd_bias[i] = imxy_rnd - imxy_GT_rnd
        uix_rnd_bias[i] = uix_rnd - uix_GT_rnd
        uiy_rnd_bias[i] = uiy_rnd - uiy_GT_rnd
        imx_rnd_bias[i] = imx_rnd - imx_GT_rnd
        imy_rnd_bias[i] = imy_rnd - imy_GT_rnd
        ri_rnd_bias[i] = ri_rnd - ri_GT_rnd
        si_rnd_bias[i] = si_rnd - si_GT_rnd
        HM_rnd_bias[i] = HM_rnd - HM_GT_rnd
        HXY_rnd_bias[i] = HXY_rnd - HXY_GT_rnd
        HMXY_rnd_bias[i] = HMXY_rnd - HMXY_GT_rnd
        imxy_ind_rnd_bias[i] = imxy_ind_rnd - imxy_ind_GT_rnd
        HM_ind_rnd_bias[i] = HM_ind_rnd - HM_ind_GT_rnd
        HXY_ind_rnd_bias[i] = HXY_ind_rnd - HXY_ind_GT_rnd
        HMXY_ind_rnd_bias[i] = HMXY_ind_rnd - HMXY_ind_GT_rnd

    # Remove the bias
    imxy = imxy - np.mean(imxy_rnd_bias)
    uix = uix - np.mean(uix_rnd_bias)
    uiy = uiy - np.mean(uiy_rnd_bias)
    imx = imx - np.mean(imx_rnd_bias)
    imy = imy - np.mean(imy_rnd_bias)
    ri = ri - np.mean(ri_rnd_bias)
    si = si - np.mean(si_rnd_bias)
    HM = HM - np.mean(HM_rnd_bias)
    HXY = HXY - np.mean(HXY_rnd_bias)
    HMXY = HMXY - np.mean(HMXY_rnd_bias)
    imxy_ind = imxy_ind - np.mean(imxy_ind_rnd_bias)
    HM_ind = HM_ind - np.mean(HM_ind_rnd_bias)
    HXY_ind = HXY_ind - np.mean(HXY_ind_rnd_bias)
    HMXY_ind = HMXY_ind - np.mean(HMXY_ind_rnd_bias)

    return imxy, uix, uiy, imx, imy, ri, si, HM, HXY, HMXY,imxy_ind, HM_ind, HXY_ind, HMXY_ind

def shuffle_subtr_bias_correction_routine(input_data, T, ntrials, seed, cov_method='numpy', precomp_cov=None):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M
    if precomp_cov is None:
        if cov_method == 'numpy':
            cov_hat = np.cov(input_data.T)
        elif cov_method == 'shrink':
            lw = LedoitWolf()
            lw.fit(input_data)  # A: rows are samples, columns are features
            cov_hat = lw.covariance_
        elif cov_method =='lasso':
            imputer = SimpleImputer(strategy='mean')
            X_clean = imputer.fit_transform(input_data)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_clean)
            model = GraphicalLassoCV()
            model.fit(X_scaled)
            cov_hat = model.covariance_
        elif cov_method == 'osa':
            cov_hat, _ = oas(input_data)
        else:
            cov_hat=precomp_cov

    imxy, imx, imy, ri, si, uix_uiy, union_info, HM, HXY, HMXY = inf_from_cov(cov_hat)

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
    uix_uiy_rnd_bias = np.zeros(T)
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
        # if cov_method == 'numpy':
        #     cov_shuff = np.cov(input_data_shuffled.T)
        # elif cov_method == 'shrink':
        #     lw = LedoitWolf()
        #     lw.fit(input_data_shuffled)  # A: rows are samples, columns are features
        #     cov_shuff = lw.covariance_
        # elif cov_method == 'lasso':
        #     imputer = SimpleImputer(strategy='mean')
        #     X_clean = imputer.fit_transform(input_data_shuffled)
        #     scaler = StandardScaler()
        #     X_scaled = scaler.fit_transform(X_clean)
        #     model = GraphicalLassoCV()
        #     model.fit(X_scaled)
        #     cov_shuff = model.covariance_
        # elif cov_method == 'osa':
        #     cov_shuff, _ = oas(input_data_shuffled)

        imxy_shuff, imx_shuff, imy_shuff, ri_shuff, si_shuff, uix_uiy_shuff, union_info_shuff, HM_shuff, HXY_shuff, HMXY_shuff = inf_from_cov(cov_shuff)
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
        uix_uiy_rnd_bias[i] = uix_uiy_shuff
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
    imxy_corr =  imxy - Jointinfo_analytical_bias # imxy - np.mean(imxy_rnd_bias) #
    imx_corr  = imx - Info_analytical_bias # imx - np.mean(imx_rnd_bias) #
    imy_corr  = imy - Info_analytical_bias # imy - np.mean(imy_rnd_bias) #
    union_corr = union_info - np.mean(union_rnd_bias)
    uix_corr = union_corr - imy_corr  # uix_uiy - np.mean(uix_uiy_rnd_bias)
    uiy_corr = union_corr - imx_corr  # uix_uiy - np.mean(uix_uiy_rnd_bias)
    ri_corr = imx_corr+imy_corr - union_corr #ri - np.mean(ri_rnd_bias)
    si_corr = imxy_corr - union_corr #si - np.mean(si_rnd_bias)
    HM = HM - np.mean(HM_rnd_bias)
    HXY = HXY - np.mean(HXY_rnd_bias)
    HMXY = HMXY - np.mean(HMXY_rnd_bias)
    imxy_ind = imxy_ind - np.mean(imxy_ind_rnd_bias)
    HM_ind = HM_ind - np.mean(HM_ind_rnd_bias)
    HXY_ind = HXY_ind - np.mean(HXY_ind_rnd_bias)
    HMXY_ind = HMXY_ind - np.mean(HMXY_ind_rnd_bias)

    return imxy_corr, uix_corr, uiy_corr, imx_corr, imy_corr, ri_corr, si_corr, HM, HXY, HMXY, imxy_ind, HM_ind, HXY_ind, HMXY_ind

def uniform_bias_correction_routine(input_data, T, ntrials, seed, cov_method='numpy'):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M
    if cov_method == 'numpy':
        cov_hat = np.cov(input_data.T)
    elif cov_method == 'shrink':
        lw = LedoitWolf()
        lw.fit(input_data)  # A: rows are samples, columns are features
        cov_hat = lw.covariance_
    elif cov_method =='lasso':
        imputer = SimpleImputer(strategy='mean')
        X_clean = imputer.fit_transform(input_data)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_clean)
        model = GraphicalLassoCV()
        model.fit(X_scaled)
        cov_hat = model.covariance_
    elif cov_method == 'osa':
        cov_hat, _ = oas(input_data)
    # cov_hat = np.cov(input_data.T)



    # model = GraphicalLassoCV()
    # model.fit(input_data)
    # cov_hat = model.covariance_

    # cov_hat, _  = oas(input_data)

    Sigma_MY_combined = np.block([
        [cov_hat[:dm, :dm], cov_hat[:dm, dm + dx:]],
        [cov_hat[dm + dx:, :dm], cov_hat[dm + dx:, dm + dx:]]
    ])

    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HX = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:dm + dx, dm:dm + dx])[1] + (dx) * (1 + np.log(2 * np.pi)))
    HY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm + dx:, dm + dx:])[1] + (dy) * (1 + np.log(2 * np.pi)))
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    HMX = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm + dx, :dm + dx])[1] + (dm + dx) * (1 + np.log(2 * np.pi)))
    HMY = 0.5 / np.log(2) * (la.slogdet(Sigma_MY_combined)[1] + (dm + dy) * (1 + np.log(2 * np.pi)))
    # Compute biased PID estimates
    try:
        # Venkatesh bias correction approach
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy, sample_size=ntrials ,unbiased = True)

        # ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
    except:  # Mainly to catch LinAlgWarning's
        imxy, uix, uiy, ri, si = [np.nan,] * 5
    
    # uix_uiy = uix + uiy
    imx = uix + ri
    imy = uiy + ri
                        
    cov_hat_ind = np.copy(cov_hat)                                                                
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm+dx, dm+dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm+dx:, dm:dm+dx] = sig_xy[dx:, :dx]
                                             
    # Venkatesh bias correction approach                        
    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy, sample_size=ntrials, unbiased = True)
    
    # ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]                                        
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) ) 
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )


        
    return imxy, uix, uiy, imx, imy, ri, si, HM, HXY, HMXY,imxy_ind, HM_ind, HXY_ind, HMXY_ind            

def No_bias_correction_routine(input_data, T, ntrials, seed, cov_method='numpy', precomp_cov=None):
    # Determine M from the data dimensions (assumes 3*M columns)
    M = input_data.shape[1] // 3
    dm = dx = dy = M

    if precomp_cov is None:
        if cov_method == 'numpy':
            cov_hat = np.cov(input_data.T)
        elif cov_method == 'shrink':
            lw = LedoitWolf()
            lw.fit(input_data)  # A: rows are samples, columns are features
            cov_hat = lw.covariance_
        elif cov_method == 'lasso':
            imputer = SimpleImputer(strategy='mean')
            X_clean = imputer.fit_transform(input_data)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X_clean)
            model = GraphicalLassoCV()
            model.fit(X_scaled)
            cov_hat = model.covariance_
        elif cov_method == 'osa':
            cov_hat, _ = oas(input_data)
        else:
            cov_hat = precomp_cov
    # cov_hat = np.cov(input_data.T)

    # lw = LedoitWolf()
    # lw.fit(input_data)  # A: rows are samples, columns are features
    # cov_hat = lw.covariance_

    # model = GraphicalLassoCV()
    # model.fit(input_data)
    # cov_hat = model.covariance_

    # cov_hat, _ = oas(input_data)


    Sigma_MY_combined = np.block([
        [cov_hat[:dm, :dm], cov_hat[:dm, dm + dx:]],
        [cov_hat[dm + dx:, :dm], cov_hat[dm + dx:, dm + dx:]]
    ])

    HM = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)))
    HX = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:dm + dx, dm:dm + dx])[1] + (dx) * (1 + np.log(2 * np.pi)))
    HY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm + dx:, dm + dx:])[1] + (dy) * (1 + np.log(2 * np.pi)))
    HXY = 0.5 / np.log(2) * (la.slogdet(cov_hat[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)))
    HMXY = 0.5 / np.log(2) * (la.slogdet(cov_hat)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)))
    HMX = 0.5 / np.log(2) * (la.slogdet(cov_hat[:dm + dx, :dm + dx])[1] + (dm + dx) * (1 + np.log(2 * np.pi)))
    HMY = 0.5 / np.log(2) * (la.slogdet(Sigma_MY_combined)[1] + (dm + dy) * (1 + np.log(2 * np.pi)))
    # Compute biased PID estimates
    try:
        # Venkatesh bias correction approach
        ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy, sample_size=ntrials ,unbiased = False)

        # ret = exact_gauss_tilde_pid(cov_hat, dm, dx, dy)
        imxy, uix, uiy, ri, si = ret[2], *ret[-4:]
    except:  # Mainly to catch LinAlgWarning's
        imxy, uix, uiy, ri, si = [np.nan,] * 5
    
    # uix_uiy = uix + uiy
    imx = uix + ri
    imy = uiy + ri
                        
    cov_hat_ind = np.copy(cov_hat)                                                                
    sig_m = cov_hat[:dm, :dm]
    sig_xy = cov_hat[dm:, dm:]
    sig_xy_m = cov_hat[dm:, :dm]
    sig_xy = sig_xy_m @ solve(sig_m, sig_xy_m.T)
    cov_hat_ind[dm:dm+dx, dm+dx:] = sig_xy[:dx, dx:]
    cov_hat_ind[dm+dx:, dm:dm+dx] = sig_xy[dx:, :dx]
                                             
    # Venkatesh bias correction approach                        
    ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy, sample_size=ntrials, unbiased = False)
    
    # ret_ind = exact_gauss_tilde_pid(cov_hat_ind, dm, dx, dy)
    imxy_ind, uix_ind, uiy_ind, ri_ind, si_ind = ret_ind[2], *ret_ind[-4:]                                        
    HM_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[:dm, :dm])[1] + (dm) * (1 + np.log(2 * np.pi)) ) 
    HXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind[dm:, dm:])[1] + (dx + dy) * (1 + np.log(2 * np.pi)) )
    HMXY_ind = 0.5 / np.log(2) * (la.slogdet(cov_hat_ind)[1] + (dm + dx + dy) * (1 + np.log(2 * np.pi)) )

    # imx = HM + HX - HMX
    # imy = HM + HY - HMY
    # imxy = HM + HXY - HMXY
    return imxy, uix, uiy, imx, imy, ri, si, HM, HXY, HMXY,imxy_ind, HM_ind, HXY_ind, HMXY_ind            