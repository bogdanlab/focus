import logging

import numpy as np
import pandas as pd
import pyfocus
import scipy.linalg as lin

from itertools import chain, combinations
from numpy.linalg import multi_dot as mdot


__all__ = ["fine_map"]

def create_output(foo):
    pass


def estimate_cor(wcollection, resid_vars, ldmat, intercept=False):
    # this needs to be made much more robust but okay for now
    wmat = wcollection.get_all_weights()
    wcov = mdot([wmat.T, ldmat, wmat])
    scale = np.diag(1 / (wcov.std(axis=0) * np.sqrt(resid_vars)))

    wcor = mdot([scale, wcov, scale])

    if intercept:
        inter = mdot([scale, wcov, ldmat])
        return wcor, inter
    else:
        return wcor, None


def assoc_test(weights, gwas, ldmat, heterogeneity=False):
    # this needs to be made much more robust but okay for now
    p = ldmat.shape[0]

    assoc = np.dot(weights.w, gwas.Z)
    if heterogeneity:
        resid = assoc - gwas.Z
        resid_var = mdot([resid, np.linalg.inv(ldmat), resid]) / p
    else:
        resid_var = 1

    se = np.sqrt(resid_var * mdot([weights.w, ldmat, weights.w]))

    # z-score is imputed score standardized by sqrt of estimated variance
    zscore = assoc / se

    return zscore, resid_var


def get_resid(zscores, swld, wcor):
    m, m = wcor.shape
    m, p = swld.shape

    # create mean factor
    intercept = swld.dot(np.ones(p))

    # estimate under the null for variance components, i.e. V = SW LD SW
    wcor_inv = np.linalg.svd(wcor)
    rank = np.linalg.matrix_rank(wcor)

    numer = mdot([intercept.T, wcor_inv, zscores])
    denom = mdot([intercept.T, wcor_inv, intercept])
    alpha = numer / denom
    resid = zscores - intercept * alpha

    # really this should be divided by the rank of wcor - 1
    s2 = mdot([resid, wcor_inv, resid]) / (rank - 1)
    inter_se = np.sqrt(s2 / denom)
    inter_z = alpha / inter_se

    return resid, inter_z


def bayes_factor(zscores, idx_set, wcor, prior_chisq, prb, use_log=True):
    idx_set = np.array(idx_set)

    m = len(zscores)

    # only need genes in the causal configuration using FINEMAP BF trick
    nc = len(idx_set)
    cur_chi2 = prior_chisq / nc

    cur_wcor = wcor[idx_set].T[idx_set].T
    cur_zscores = zscores[idx_set]

    # compute SVD for robust estimation
    if nc > 1:
        cur_U, cur_EIG, _ = lin.svd(cur_wcor)
        scaled_chisq = (cur_zscores.T.dot(cur_U)) ** 2
    else:
        cur_U, cur_EIG = 1, cur_wcor
        scaled_chisq = cur_zscores ** 2

    # log BF + log prior
    cur_bf = 0.5 * -np.sum(np.log(1 + cur_chi2 * cur_EIG)) + \
        0.5 *  np.sum((cur_chi2 / (1 + cur_chi2 * cur_EIG)) * scaled_chisq) + \
        nc * np.log(prb) + (m - nc) * np.log(1 - prb)

    if use_log:
        return cur_bf
    else:
        return np.exp(cur_bf)


def fine_map(gwas, wcollection, ref_geno, intercept=False, heterogeneity=False, ridge=0.1, prior_prob=1e-3, prior_var=40):
    log = logging.getLogger(pyfocus.LOG)
    log.info("Starting fine-mapping at region {}".format(ref_geno))

    # align gwas with weights / ref
    merged_snps = ref_geno.overlap_gwas(gwas)
    ref_snps = merged_snps.loc[~pd.isna(merged_snps.i)]

    # compute linkage-disequilibrium estimate
    log.debug("Estimating LD for {} SNPs".format(len(ref_snps)))
    ldmat = ref_geno.estimate_ld(ref_snps, adjust=ridge)

    zscores = []
    resid_vars = []

    # run local TWAS
    # we need to check overlap somewhere in here
    for weights in wcollection:
        result, resid_var = assoc_test(weights, gwas, ldmat, heterogeneity)
        zscores.append(result)
        resid_vars.append(resid_var)

    # perform fine-mapping
    zscores = np.array(zscores)
    log.debug("Estimating local TWAS correlation structure")
    wcor, swld = estimate_cor(wcollection, ldmat, resid_vars, intercept)

    m = len(zscores)
    rm = range(m)
    pips = np.zeros(m)

    if intercept:
        zscores, inter_z = get_resid(zscores, swld, wcor)
    else:
        inter_z = None

    k = m
    null_res = m * np.log(1 - prior_prob)
    marginal = null_res
    # enumerate all subsets
    for subset in chain.from_iterable(combinations(rm, n) for n in range(1, k + 1)):

        local = bayes_factor(zscores, subset, wcor, prior_var, prior_prob)

        # keep track for marginal likelihood
        marginal = np.logaddexp(local, marginal)

        # marginalize the posterior for marginal-posterior on causals
        for idx in subset:
            if pips[idx] == 0:
                pips[idx] = local
            else:
                pips[idx] = np.logaddexp(pips[idx], local)

    pips = np.exp(pips - marginal)
    null_res = np.exp(null_res - marginal)

    # clean up and return results
    df = create_output(None)
    log.info("Completed fine-mapping at region {}".format(ref_geno))

    return df
