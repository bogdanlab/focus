import logging

import numpy as np
import pandas as pd
import pyfocus as pf
import scipy.linalg as lin

from itertools import chain, combinations
from numpy.linalg import multi_dot as mdot


__all__ = ["fine_map"]


def create_output(foo):
    return foo


def align_data(gwas, ref_geno, wcollection, ridge=0.1):
    """
    Align and merge gwas, LD reference, and eQTL weight data to the same reference alleles.

    :param gwas: pyfocus.GWAS object containing a risk region
    :param ref_geno:  pyfocus.LDRefPanel object containing reference genotypes at risk region
    :param wcollection: pandas.DataFrame object containing overlapping eQTL weights
    :param ridge: ridge adjustment for LD estimation (default = 0.1)

    :return: tuple of aligned GWAS, eQTL weight-marix W, gene-names list, LD-matrix V
    """
    log = logging.getLogger(pf.LOG)

    # align gwas with ref snps
    merged_snps = ref_geno.overlap_gwas(gwas)
    ref_snps = merged_snps.loc[~pd.isna(merged_snps.i)]

    # flip Zscores to match reference panel
    ref_snps[pf.GWAS.ZCOL] = pf.flip_values(ref_snps[pf.GWAS.ZCOL].values,
                                            ref_snps[pf.GWAS.A1COL],
                                            ref_snps[pf.GWAS.A2COL],
                                            ref_snps[pf.LDRefPanel.A1COL],
                                            ref_snps[pf.LDRefPanel.A2COL])

    # collapse the gene models into a single weight matrix
    genes = []
    final_df = None
    for eid, model in wcollection.groupby("ens_gene_id"):
        log.debug("Aligning weights for gene {}".format(eid))

        # merge local model with the reference panel
        # effect_allele alt_allele effect
        m_merged = pd.merge(ref_snps, model, how="inner", left_on=pf.GWAS.SNPCOL, right_on="rsid")

        # make sure effects are for same ref allele as GWAS + reference panel
        m_merged["effect"] = pf.flip_values(m_merged["effect"].values,
                                            m_merged["effect_allele"],
                                            m_merged["alt_allele"],
                                            m_merged[pf.LDRefPanel.A1COL],
                                            m_merged[pf.LDRefPanel.A2COL])

        # skip genes whose overlapping weights are all 0s
        if np.isclose(np.var(m_merged["effect"]), 0):
            continue

        # keep track of name
        genes.append(eid)

        # perform a union (outer merge) to build the aligned/flipped weight (possibly jagged) matrix
        if final_df is None:
            final_df = m_merged[[pf.GWAS.SNPCOL, "effect"]]
        else:
            final_df = pd.merge(final_df, m_merged[[pf.GWAS.SNPCOL, "effect"]], how="outer", on="SNP")

    # break out early
    if len(genes) == 0:
        return None

    # final align back with GWAS + reference panel
    ref_snps = pd.merge(ref_snps, final_df, how="inner", on=pf.GWAS.SNPCOL)

    # compute linkage-disequilibrium estimate
    log.debug("Estimating LD for {} SNPs".format(len(ref_snps)))
    ldmat = ref_geno.estimate_ld(ref_snps, adjust=ridge)

    # subset down to just actual GWAS data
    ref_snps = ref_snps[pf.GWAS.REQ_COLS]

    # need to replace NA with 0 due to jaggedness across genes
    wmat = final_df.loc[:, final_df.columns != "SNP"].values
    wmat[np.isnan(wmat)] = 0.0

    return ref_snps, wmat, genes, ldmat


def estimate_cor(wmat, ldmat, intercept=False):
    """
    Estimate the sample correlation structure for predicted expression.

    :param wmat: numpy.ndarray eQTL weight matrix for a risk region
    :param ldmat: numpy.ndarray LD matrix for a risk region
    :param intercept: bool to return the intercept variable or not

    :return: tuple (pred_expr correlation, intercept variable; None if intercept=False)
    """
    wcov = mdot([wmat.T, ldmat, wmat])
    scale = np.diag(1 / np.sqrt(np.diag(wcov))
    wcor = mdot([scale, wcov, scale])

    if intercept:
        inter = mdot([scale, wcov, ldmat])
        return wcor, inter
    else:
        return wcor, None


def assoc_test(weights, gwas, ldmat, heterogeneity=False):
    """
    TWAS association test.

    :param weights: numpy.ndarray of eQTL weights
    :param gwas: pyfocus.GWAS object
    :param ldmat: numpy.ndarray LD matrix
    :param heterogeneity:  bool estimate variance from multiplicative random effect

    :return: tuple (beta, se)
    """
    p = ldmat.shape[0]

    assoc = np.dot(weights, gwas.Z)
    if heterogeneity:
        resid = assoc - gwas.Z
        resid_var = mdot([resid, np.linalg.inv(ldmat), resid]) / p
    else:
        resid_var = 1

    se = np.sqrt(resid_var * mdot([weights, ldmat, weights]))

    return assoc, se


def get_resid(zscores, swld, wcor):
    """
    Regress out the average pleiotropic signal tagged by TWAS at the region

    :param zscores: numpy.ndarray TWAS zscores
    :param swld: numpy.ndarray intercept variable
    :param wcor: numpy.ndarray predicted expression correlation

    :return: tuple (residual TWAS zscores, intercept z-score)
    """
    m, m = wcor.shape
    m, p = swld.shape

    # create mean factor
    intercept = swld.dot(np.ones(p))

    # estimate under the null for variance components, i.e. V = SW LD SW
    wcor_inv, rank = lin.pinvh(wcor, return_rank=True)

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
    """
    Compute the Bayes Factor for the evidence that a set of genes explain the observed association signal under the
    correlation structure.


    :param zscores: numpy.ndarray TWAS zscores
    :param idx_set: list the indices representing the causal gene-set
    :param wcor: numpy.ndarray predicted expression correlation
    :param prior_chisq: float prior effect-size variance scaled by GWAS sample size
    :param prb:  float prior probability for a gene to be causal
    :param use_log: bool whether to compute the log Bayes Factor

    :return: float the Bayes Factor (log Bayes Factor if use_log = True)
    """
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


def fine_map(gwas, wcollection, ref_geno, intercept=False, heterogeneity=False, max_genes=None, ridge=0.1, prior_prob=1e-3, prior_chisq=40):
    """
    Perform a TWAS and fine-map the results.

    :param gwas: pyfocus.GWAS object for the risk region
    :param wcollection: pandas.DataFrame containing overlapping eQTL weight information for the risk region
    :param ref_geno: pyfocus.LDRefPanel object for the risk region
    :param intercept: bool flag to estimate the average TWAS signal due to tagged pleiotropy
    :param heterogeneity: bool flag to compute sample variance in TWAS test assuming multiplicative random effect
    :param max_genes: int or None the maximum number of genes to include in any given causal configuration. None if all genes
    :param ridge: float ridge adjustment for LD estimation (default = 0.1)
    :param prior_prob: float prior probability for a gene to be causal
    :param prior_chisq: float prior effect-size variance scaled by GWAS sample size

    :return: pandas.DataFrame containing the TWAS statistics and fine-mapping results
    """
    log = logging.getLogger(pf.LOG)
    log.info("Starting fine-mapping at region {}".format(ref_geno))

    # align all GWAS, LD reference, and overlapping molecular weights
    parameters = align_data(gwas, ref_geno, wcollection, ridge=ridge)

    if parameters is None:
        # TODO: log and return empty result
        pass
    else:
        gwas, wmat, genes, ldmat = parameters

    zscores = []
    # run local TWAS
    # we need to check overlap somewhere in here
    for idx, weights in enumerate(wmat.T):
        log.debug("Computing TWAS association statistic for gene {}".format(genes[idx]))
        beta, se = assoc_test(weights, gwas, ldmat, heterogeneity)
        zscores.append(beta / se)

    # perform fine-mapping
    zscores = np.array(zscores)

    log.debug("Estimating local TWAS correlation structure")
    wcor, swld = estimate_cor(wmat, ldmat, intercept)

    m = len(zscores)
    rm = range(m)
    pips = np.zeros(m)

    if intercept:
        # should really be done at the SNP level first ala Barfield et al 2018
        log.debug("Regressing out average tagged pleiotropic associations")
        zscores, inter_z = get_resid(zscores, swld, wcor)
    else:
        inter_z = None

    k = m if max_genes is None else max_genes
    null_res = m * np.log(1 - prior_prob)
    marginal = null_res
    # enumerate all subsets
    for subset in chain.from_iterable(combinations(rm, n) for n in range(1, k + 1)):

        local = bayes_factor(zscores, subset, wcor, prior_chisq, prior_prob)

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
