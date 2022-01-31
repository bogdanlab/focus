import logging

import numpy as np
import pandas as pd
import scipy.linalg as lin

from itertools import chain, combinations
from numpy.linalg import multi_dot as mdot
import pyfocus as pf

# __all__ = ["fine_map", "me_fine_map", "num_convert"]
__all__ = ["fine_map", "num_convert"]


# def add_credible_set(df, credible_set=0.9):
#     """
#     Compute the credible gene set and add it to the dataframe.
#
#     :param df: pandas.DataFrame containing TWAS summary results
#     :param credible_set: float sensitivity to compute the credble set
#
#     :return: pandas.DataFrame containing TWAS summary results augmented with `in_cred_set` flag
#     """
#     df = df.sort_values(by=["pip"], ascending=False)
#
#     # add credible set flag
#     psum = np.sum(df.pip.values)
#     csum = np.cumsum(df.pip.values / psum)
#     in_cred_set = np.zeros(len(csum), dtype=int)
#     for idx, npip in enumerate(csum):
#         in_cred_set[idx] = True
#         if npip >= credible_set:
#             break
#
#     df["in_cred_set"] = in_cred_set
#
#     return df

def add_credible_set(df, credible_set=0.9):
    """
    Compute the credible gene set and add it to the dataframe.

    :param df: pandas.DataFrame containing TWAS summary results
    :param credible_set: float sensitivity to compute the credble set

    :return: pandas.DataFrame containing TWAS summary results augmented with `in_cred_set` flag
    """
    n_pips = np.sum(df.columns.str.contains("pips"))

    for i in range(n_pips):
        if i == (n_pips - 1) and i > 0:
            name_pip = "pips_me"
            name_set = "in_cred_set_me"
        else:
            name_pip = f"pips_pop{i+1}"
            name_set = f"in_cred_set_pop{i+1}"

        df = df.sort_values(by=[name_pip], ascending=False)

        # add credible set flag
        psum = np.sum(df[name_pip].values)
        csum = np.cumsum(df[name_pip].values / psum)
        in_cred_set = np.zeros(len(csum), dtype=int)
        for idx, npip in enumerate(csum):
            in_cred_set[idx] = True
            if npip >= credible_set:
                break

        df[name_set] = in_cred_set

    return df

# def create_output(meta_data, attr, zscores, pips, null_res, region):
#     """
#     Creates TWAS pandas.DataFrame output.
#
#     :param meta_data: pandas.DataFrame Metadata about the gene models
#     :param attr: pandas.DataFrame Prediction performance metadata about gene models
#     :param zscores: numpy.ndarray of TWAS zscores
#     :param pips: numpy.ndarray of posterior inclusion probabilities (PIPs)
#     :param null_res: float posterior probability of the null
#     :param region: str region identifier
#
#     :return: pandas.DataFrame TWAS summary results
#     """
#
#     # merge attributes
#     df = pd.merge(meta_data, attr, left_on="model_id", right_index=True)
#     df["twas_z"] = zscores
#     df["pip"] = pips
#     df["in_cred_set"] = 0
#     df["region"] = region
#
#     # sort by tx start site and we're good to go
#     df = df.sort_values(by="tx_start")
#
#     # chrom is listed twice (once from weight and once from molfeature)
#     idxs = np.where(df.columns == "chrom")[0]
#     if len(idxs) > 1:
#         df = df.iloc[:, [j for j, c in enumerate(df.columns) if j != idxs[0]]]
#
#     # drop model-id
#     df = df.drop("model_id", axis=1)
#
#     # add null model result
#     null_dict = dict()
#     for c in df.columns:
#         null_dict[c] = None
#
#     null_dict["ens_gene_id"] = "NULL.MODEL"
#     null_dict["mol_name"] = "NULL"
#     null_dict["type"] = "NULL"
#     null_dict["pip"] = null_res
#     null_dict["region"] = region
#     null_dict["twas_z"] = 0
#     null_dict["chrom"] = df["chrom"].values[0]
#     df = df.append(null_dict, ignore_index=True)
#
#     return df

def create_output(meta_data, attr, null_res, region):
    """
    Creates TWAS pandas.DataFrame output.

    :param meta_data: pandas.DataFrame Metadata about the gene models
    :param attr: pandas.DataFrame Prediction performance metadata about gene models
    :param null_res: float posterior probability of the null
    :param region: str region identifier

    :return: pandas.DataFrame TWAS summary results
    """
    # merge attributes
    n_pop = len(meta_data)
    join_col = ["ens_gene_id", "tissue", "model_id"]
    pick_col = ["ens_gene_id", "tissue", "model_id", "inference", "cv.R2", "cv.R2.pval", "inter_z", "twas_z", "pips"]
    for i in range(n_pop):
        df_tmp = pd.merge(meta_data[i], attr[i], left_on="model_id", right_index=True)
        if i != 0:
            df_tmp = df_tmp[pick_col]
        df_tmp = df_tmp.rename(columns = {"cv.R2": f"cv.R2_pop{i+1}",
                                          "cv.R2.pval": f"cv.R2.pval_pop{i+1}",
                                          "inference": f"inference_pop{i+1}",
                                          "twas_z": f"twas_z_pop{i+1}",
                                          "inter_z": f"inter_z_pop{i+1}",
                                          "pips": f"pips_pop{i+1}"})
        df_tmp[f"ldregion_pop{i+1}"] = region[i]
        if i == 0:
            df = df_tmp
        else:
            if len(df) == len(df_tmp):
                df = pd.merge(df, df_tmp, how = "left", on = join_col)
            else:
                raise ValueError(f"Cannot column binds output due to different rows for populations.")
                return None

    # sort by tx start site and we're good to go
    df = df.sort_values(by="tx_start")

    # chrom is listed twice (once from weight and once from molfeature)
    idxs = np.where(df.columns.str.contains("chrom"))[0]
    if len(idxs) > 1:
        df = df.iloc[:, [j for j, c in enumerate(df.columns) if j != idxs[0]]]

    df.columns.values[np.where(df.columns.str.contains("chrom"))] = "chrom"

    # drop model-id
    # idx = [c for c in df.columns if "model_id" not in c]
    # df = df[idx]
    df.pop("model_id")

    # add null model result
    null_dict = dict()
    for c in df.columns:
        null_dict[c] = None

    null_dict["ens_gene_id"] = "NULL.MODEL"
    null_dict["mol_name"] = "NULL"
    null_dict["type"] = "NULL"
    null_dict["chrom"] = df["chrom"].values[0]

    for i in range(n_pop):
        null_dict[f"twas_z_pop{i+1}"] = 0
        null_dict[f"pips_pop{i+1}"] = null_res[i]
        null_dict[f"ldregion_pop{i+1}"] = region[i]
        if i == (n_pop - 1) and i > 0:
            null_dict["pips_me"] = null_res[i+1]

    df = df.append(null_dict, ignore_index=True)

    return df

# def align_data(gwas, ref_geno, wcollection, ridge=0.1):
#     """
#     Align and merge gwas, LD reference, and eQTL weight data to the same reference alleles.
#     :param gwas: pyfocus.GWAS object containing a risk region
#     :param ref_geno:  pyfocus.LDRefPanel object containing reference genotypes at risk region
#     :param wcollection: pandas.DataFrame object containing overlapping eQTL weights
#     :param ridge: ridge adjustment for LD estimation (default = 0.1)
#     :return: tuple of aligned GWAS, eQTL weight-matrix W, gene-names list, LD-matrix V
#     """
#     log = logging.getLogger(pf.LOG)
#
#     # align gwas with ref snps
#     merged_snps = ref_geno.overlap_gwas(gwas)
#     if len(merged_snps) == 0:
#         log.info("No overlap between LD reference and GWAS")
#         return None
#
#     ref_snps = merged_snps.loc[~pd.isna(merged_snps.i)]
#
#     # filter out mis-matched SNPs
#     matched = pf.check_valid_alleles(ref_snps[pf.GWAS.A1COL],
#                                      ref_snps[pf.GWAS.A2COL],
#                                      ref_snps[pf.LDRefPanel.A1COL],
#                                      ref_snps[pf.LDRefPanel.A2COL])
#     n_miss = sum(np.logical_not(matched))
#     if n_miss > 0:
#         log.debug("Pruned {} SNPs due to invalid allele pairs between GWAS/RefPanel.".format(n_miss))
#
#     ref_snps = ref_snps.loc[matched]
#
#     # flip Zscores to match reference panel
#     ref_snps[pf.GWAS.ZCOL] = pf.flip_alleles(ref_snps[pf.GWAS.ZCOL].values,
#                                              ref_snps[pf.GWAS.A1COL],
#                                              ref_snps[pf.GWAS.A2COL],
#                                              ref_snps[pf.LDRefPanel.A1COL],
#                                              ref_snps[pf.LDRefPanel.A2COL])
#
#     # collapse the gene models into a single weight matrix
#     idxs = []
#     final_df = None
#     for eid, model in wcollection.groupby(["ens_gene_id", "tissue", "inference", "ref_name"]):
#         log.debug("Aligning weights for gene {}".format(eid))
#
#         # merge local model with the reference panel
#         # effect_allele alt_allele effect
#         m_merged = pd.merge(ref_snps, model, how="inner", left_on=pf.GWAS.SNPCOL, right_on="snp")
#
#         m_matched = pf.check_valid_alleles(m_merged["effect_allele"],
#                                            m_merged["alt_allele"],
#                                            m_merged[pf.LDRefPanel.A1COL],
#                                            m_merged[pf.LDRefPanel.A2COL])
#
#         n_miss = sum(np.logical_not(m_matched))
#         if n_miss > 0:
#             log.debug("Gene {} pruned {} SNPs due to invalid allele pairs between weight-db/GWAS.".format(eid, n_miss))
#
#         m_merged = m_merged.loc[m_matched]
#
#         # make sure effects are for same ref allele as GWAS + reference panel
#         m_merged["effect"] = pf.flip_alleles(m_merged["effect"].values,
#                                              m_merged["effect_allele"],
#                                              m_merged["alt_allele"],
#                                              m_merged[pf.LDRefPanel.A1COL],
#                                              m_merged[pf.LDRefPanel.A2COL])
#
#         # skip genes whose overlapping weights are all 0s
#         if len(m_merged) > 1 and all(np.isclose(m_merged["effect"], 0)):
#             log.debug("Gene {} has only zero-weights. This will break variance estimate. Skipping.".format(eid))
#             continue
#
#         # skip genes that do not have weights at referenced SNPs
#         if all(pd.isnull(m_merged["effect"])):
#             log.debug("Gene {} has no overlapping weights. Skipping.".format(eid))
#             continue
#
#         # keep model_id around to grab other attributes (pred-R2, etc) later on
#         cur_idx = model.index[0]
#         idxs.append(cur_idx)
#
#         # perform a union (outer merge) to build the aligned/flipped weight (possibly jagged) matrix
#         if final_df is None:
#             final_df = m_merged[[pf.GWAS.SNPCOL, "effect"]]
#             final_df = final_df.rename(index=str, columns={"effect": "model_{}".format(cur_idx)})
#         else:
#             final_df = pd.merge(final_df, m_merged[[pf.GWAS.SNPCOL, "effect"]], how="outer", on="SNP")
#             final_df = final_df.rename(index=str, columns={"effect": "model_{}".format(cur_idx)})
#
#     # break out early
#     if len(idxs) == 0:
#         log.info("No weights overlapped GWAS data")
#         return None
#
#     # final align back with GWAS + reference panel
#     ref_snps = pd.merge(ref_snps, final_df, how="inner", on=pf.GWAS.SNPCOL)
#
#     # compute linkage-disequilibrium estimate
#     log.debug("Estimating LD for {} SNPs".format(len(ref_snps)))
#     ldmat = ref_geno.estimate_ld(ref_snps, adjust=ridge)
#
#     # subset down to just actual GWAS data
#     gwas = ref_snps[pf.GWAS.REQ_COLS]
#
#     # need to replace NA with 0 due to jaggedness across genes
#     wmat = ref_snps.filter(like="model").values
#     wmat[np.isnan(wmat)] = 0.0
#
#     # Meta-data on the current model
#     # what other things should we include in here?
#     meta_data = wcollection.loc[idxs,
#                                 ["ens_gene_id", "ens_tx_id", "mol_name", "tissue", "ref_name", "type", "chrom", "tx_start",
#                                 "tx_stop", "inference", "model_id"]
#     ]
#
#     # re-rorder by tx_start
#     ranks = np.argsort(meta_data["tx_start"].values)
#     wmat = wmat.T[ranks].T
#     meta_data = meta_data.iloc[ranks]
#
#     return gwas, wmat, meta_data, ldmat

def align_data(gwas, ref_geno, wcollection, min_r2pred=0.7, max_impute=0.5, ridge=0.1):
    """
    Align and merge gwas, LD reference, and eQTL weight data to the same reference alleles.

    :param gwas: pyfocus.GWAS object containing a risk region
    :param ref_geno: pyfocus.LDRefPanel object containing reference genotypes at risk region
    :param wcollection: pandas.DataFrame object containing overlapping eQTL weights
    :param min_r2pred: minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores (default = 0.7)
    :param max_impute: maximum fraction of SNPs allowed to be missing per gene, and will be imputed using LD (default = 0.5)
    :param ridge: ridge adjustment for LD estimation (default = 0.1)

    :return: tuple of aligned GWAS, eQTL weight-matrix W, gene-names list, LD-matrix V
    """
    log = logging.getLogger(pf.LOG)
    # align gwas with ref snps
    merged_snps = ref_geno.overlap_gwas(gwas, enable_impg = True)

    if pd.isna(merged_snps[pf.GWAS.SNPCOL]).all():
        log.warning("No overlap between LD reference and GWAS. Skipping.")
        return None

    # to make sure no NA in LD Ref
    ref_snps = merged_snps.loc[~pd.isna(merged_snps.i)]

    # filter out mis-matched SNPs
    matched = pf.check_valid_alleles(ref_snps[pf.GWAS.A1COL],
                                     ref_snps[pf.GWAS.A2COL],
                                     ref_snps[pf.LDRefPanel.A1COL],
                                     ref_snps[pf.LDRefPanel.A2COL],
                                     enable_impg = True)

    n_miss = sum(np.logical_not(matched))
    if n_miss > 0:
        log.debug(f"Pruned {n_miss} SNPs due to invalid allele pairs between GWAS/RefPanel.")

    ref_snps = ref_snps.loc[matched]

    # flip Zscores to match reference panel
    ref_snps[pf.GWAS.ZCOL] = pf.flip_alleles(ref_snps[pf.GWAS.ZCOL].values,
                                             ref_snps[pf.GWAS.A1COL],
                                             ref_snps[pf.GWAS.A2COL],
                                             ref_snps[pf.LDRefPanel.A1COL],
                                             ref_snps[pf.LDRefPanel.A2COL],
                                             enable_impg = True)

    # Re-assign in case there are NA GWAS value
    ref_snps[pf.GWAS.SNPCOL] = ref_snps[pf.LDRefPanel.SNPCOL]
    ref_snps[pf.GWAS.A1COL] = ref_snps[pf.LDRefPanel.A1COL]
    ref_snps[pf.GWAS.A2COL] = ref_snps[pf.LDRefPanel.A2COL]

    # collapse the gene models into a single weight matrix
    idxs = []
    final_df = None
    for eid, model in wcollection.groupby(["ens_gene_id", "tissue", "inference", "ref_name"]):
        log.debug(f"Aligning weights for gene {eid}.")

        # merge local model with the reference panel
        # effect_allele alt_allele effect
        m_merged = pd.merge(ref_snps, model, how="inner", left_on=pf.GWAS.SNPCOL, right_on="snp")

        if len(m_merged) == 0:
            log.debug(f"Gene {eid} has no overlapping weights. Skipping.")
            continue

        # if the percentage of untyped SNP GWAS SS is more than default max_impute, skip.
        if np.sum(pd.isna(m_merged[pf.GWAS.ZCOL])) / len(m_merged) > max_impute:
            log.debug(f"Gene {eid} has {len(pd.isna(m_merged[pf.GWAS.ZCOL]))} out of {len(m_merged)} non-overlapping GWAS Z scores, which is less than the default max_impute {max_impute}. Skipping.")
            continue

        # here, it has nothing to do with GWAS columns, so emable_impg set to False
        m_matched = pf.check_valid_alleles(m_merged["effect_allele"],
                                           m_merged["alt_allele"],
                                           m_merged[pf.LDRefPanel.A1COL],
                                           m_merged[pf.LDRefPanel.A2COL],
                                           enable_impg = False)

        n_miss = sum(np.logical_not(m_matched))
        if n_miss > 0:
            log.debug(f"Gene {eid} pruned {n_miss} SNPs due to invalid allele pairs between weight-db/GWAS.")

        m_merged = m_merged.loc[m_matched]

        # make sure effects are for same ref allele as GWAS + reference panel
        # here, it has nothing to do with GWAS columns, so emable_impg set to False
        m_merged["effect"] = pf.flip_alleles(m_merged["effect"].values,
                                             m_merged["effect_allele"],
                                             m_merged["alt_allele"],
                                             m_merged[pf.LDRefPanel.A1COL],
                                             m_merged[pf.LDRefPanel.A2COL],
                                             enable_impg = False)

        # skip genes whose overlapping weights are all 0s
        if len(m_merged) > 0 and all(np.isclose(m_merged["effect"], 0)):
            log.debug(f"Gene {eid} has only zero-weights. This will break variance estimate. Skipping.")
            continue

        # skip genes that do not have weights at referenced SNPs
        if all(pd.isnull(m_merged["effect"])):
            log.debug(f"Gene {eid} has no overlapping weights. Skipping.")
            continue

        # keep model_id around to grab other attributes (pred-R2, etc) later on
        cur_idx = model.index[0]
        idxs.append(cur_idx)

        # perform a union (outer merge) to build the aligned/flipped weight (possibly jagged) matrix
        if final_df is None:
            final_df = m_merged[[pf.GWAS.SNPCOL, "effect"]]
            final_df = final_df.rename(index=str, columns={"effect": f"model_{cur_idx}"})
        else:
            final_df = pd.merge(final_df, m_merged[[pf.GWAS.SNPCOL, "effect"]], how="outer", on="SNP")
            final_df = final_df.rename(index=str, columns={"effect": f"model_{cur_idx}"})

    # break out early
    if len(idxs) == 0:
        log.warning("No weights overlapped GWAS data. Skipping.")
        return None

    # final align back with GWAS + reference panel
    ref_snps = pd.merge(ref_snps, final_df, how="inner", on=pf.GWAS.SNPCOL)

    # compute linkage-disequilibrium estimate
    log.debug(f"Estimating LD for {len(ref_snps)} SNPs.")
    ldmat = ref_geno.estimate_ld(ref_snps, adjust=ridge)

    # Run ImpG-Summary when GWAS statistics are not available for SNPs available in LD Ref
    miss_idx = ref_snps[pd.isna(ref_snps[pf.GWAS.ZCOL])].index
    nmiss_idx = ref_snps[pd.notna(ref_snps[pf.GWAS.ZCOL])].index

    if len(miss_idx) > 0:
        log.info(f"Found {len(miss_idx)} GWAS SS are missing, but in the LD reference panel. Imputing them using ImpG-Summary.")

        tmp = lin.inv(ldmat[nmiss_idx].T[nmiss_idx].T + 0.1 * np.eye(len(nmiss_idx)))
        impu_wgt = np.dot(ldmat[miss_idx].T[nmiss_idx].T, tmp)
        impu_z = np.dot(impu_wgt, ref_snps[pf.GWAS.ZCOL][nmiss_idx])
        r2_pred = np.diagonal(mdot([impu_wgt, ldmat[nmiss_idx].T[nmiss_idx].T, impu_wgt.T]))
        ref_snps.loc[miss_idx, pf.GWAS.ZCOL] = impu_z / np.sqrt(r2_pred)

        # Remove genes that have too many GWAS imputation or have SNPs whose GWAS SS could not be imputed
        all_r2pred = np.array([1] * len(ref_snps), np.float)
        all_r2pred[miss_idx] = r2_pred
        # track the index for removal genes
        rem_idx = []
        for i in range(len(idxs)):
            notna_snps = ref_snps[f"model_{idxs[i]}"].loc[pd.notna(ref_snps[f"model_{idxs[i]}"])]
            # we want to only look at those SNPs that have actual weights
            mean_r2 = np.mean(all_r2pred[notna_snps.index])
            gene_name = wcollection.loc[idxs[0]]["mol_name"]
            if mean_r2 < min_r2pred:
                log.warning(f"{gene_name} has mean GWAS Z-score imputation r2 of {mean_r2}, which is less than the default {min_r2pred}. Skipping.")
                rem_idx.append(idxs[i])
                continue
            notna_gwas = ref_snps.iloc[notna_snps.index][pf.GWAS.ZCOL]
            if np.sum(pd.isna(notna_gwas)) > 0:
                log.warning(f"{gene_name} has missing GWAS Z-scores that could not be imputed. Skipping.")
                rem_idx.append(idxs[i])

        # need to remove it from idxs (for meta data) and ref_snps (for wmat)
        idxs = [j for j in idxs if j not in rem_idx]
        for i in rem_idx:
            ref_snps = ref_snps.drop(f"model_{i}", axis = 1)

    # subset down to just actual GWAS data
    gwas = ref_snps[pf.GWAS.REQ_COLS]

    # need to replace NA with 0 due to jaggedness across genes
    wmat = ref_snps.filter(like="model").values
    wmat[np.isnan(wmat)] = 0.0

    # Meta-data on the current model
    # what other things should we include in here?
    meta_data = wcollection.loc[idxs,
                                ["ens_gene_id", "ens_tx_id", "mol_name", "tissue", "ref_name", "type", "chrom", "tx_start",
                                "tx_stop", "inference", "model_id"]]
    # remove duplicated columns
    meta_data = meta_data.loc[:, ~meta_data.columns.duplicated()]

    # re-rorder by tx_start
    ranks = np.argsort(meta_data["tx_start"].values)
    wmat = wmat.T[ranks].T
    meta_data = meta_data.iloc[ranks].reset_index(drop = True)

    return gwas, wmat, meta_data, ldmat


def estimate_cor(wmat, ldmat, intercept=False):
    """
    Estimate the sample correlation structure for predicted expression.

    :param wmat: numpy.ndarray eQTL weight matrix for a risk region
    :param ldmat: numpy.ndarray LD matrix for a risk region
    :param intercept: bool to return the intercept variable or not

    :return: tuple (pred_expr correlation, intercept variable; None if intercept=False)
    """
    wcov = mdot([wmat.T, ldmat, wmat])
    scale = np.diag(1 / np.sqrt(np.diag(wcov)))
    wcor = mdot([scale, wcov, scale])

    if intercept:
        inter = mdot([scale, wmat.T, ldmat])
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
        resid_var = mdot([resid, lin.pinvh(ldmat), resid]) / p
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

    s2 = mdot([resid, wcor_inv, resid]) / (rank - 1)
    inter_se = np.sqrt(s2 / denom)
    inter_z = alpha / inter_se

    return resid, inter_z


# def bayes_factor(zscores, idx_set, wcor, prior_chisq, prb, use_log=True):
#     """
#     Compute the Bayes Factor for the evidence that a set of genes explain the observed association signal under the
#     correlation structure.
#
#     :param zscores: numpy.ndarray TWAS zscores
#     :param idx_set: list the indices representing the causal gene-set
#     :param wcor: numpy.ndarray predicted expression correlation
#     :param prior_chisq: float prior effect-size variance scaled by GWAS sample size
#     :param prb:  float prior probability for a gene to be causal
#     :param use_log: bool whether to compute the log Bayes Factor
#
#     :return: float the Bayes Factor (log Bayes Factor if use_log = True)
#     """
#     idx_set = np.array(idx_set)
#
#     m = len(zscores)
#
#     # only need genes in the causal configuration using FINEMAP BF trick
#     nc = len(idx_set)
#     cur_chi2 = prior_chisq / nc
#
#     cur_wcor = wcor[idx_set].T[idx_set].T
#     cur_zscores = zscores[idx_set]
#
#     # compute SVD for robust estimation
#     if nc > 1:
#         cur_U, cur_EIG, _ = lin.svd(cur_wcor)
#         scaled_chisq = (cur_zscores.T.dot(cur_U)) ** 2
#     else:
#         cur_U, cur_EIG = 1, cur_wcor
#         scaled_chisq = cur_zscores ** 2
#
#     # log BF + log prior
#     cur_bf = 0.5 * -np.sum(np.log(1 + cur_chi2 * cur_EIG)) + \
#         0.5 * np.sum((cur_chi2 / (1 + cur_chi2 * cur_EIG)) * scaled_chisq) + \
#         nc * np.log(prb) + (m - nc) * np.log(1 - prb)
#
#     if use_log:
#         return cur_bf
#     else:
#         return np.exp(cur_bf)

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

    # log BF
    cur_bf = 0.5 * -np.sum(np.log(1 + cur_chi2 * cur_EIG)) + \
        0.5 * np.sum((cur_chi2 / (1 + cur_chi2 * cur_EIG)) * scaled_chisq)

    # log prior
    cur_prior = nc * np.log(prb) + (m - nc) * np.log(1 - prb)

    if use_log:
        return cur_bf, cur_prior
    else:
        return np.exp(cur_bf), np.exp(cur_prior)

# def fine_map(gwas, wcollection, ref_geno, intercept=False, heterogeneity=False, max_genes=3, ridge=0.1, prior_prob=1e-3,prior_chisq=40, credible_level=0.9, plot=False):
#     """
#     Perform a TWAS and fine-map the results.
#
#     :param gwas: pyfocus.GWAS object for the risk region
#     :param wcollection: pandas.DataFrame containing overlapping eQTL weight information for the risk region
#     :param ref_geno: pyfocus.LDRefPanel object for the risk region
#     :param intercept: bool flag to estimate the average TWAS signal due to tagged pleiotropy
#     :param heterogeneity: bool flag to compute sample variance in TWAS test assuming multiplicative random effect
#     :param max_genes: int or None the maximum number of genes to include in any given causal configuration. None if all genes
#     :param ridge: float ridge adjustment for LD estimation (default = 0.1)
#     :param prior_prob: float prior probability for a gene to be causal
#     :param prior_chisq: float prior effect-size variance scaled by GWAS sample size
#     :param credible_level: float the credible-level to compute credible gene sets (default = 0.9)
#     :param plot: bool whether or not to generate visualizations/plots at the risk region
#
#     :return: pandas.DataFrame containing the TWAS statistics and fine-mapping results if plot=False.
#         (pandas.DataFrame, list of plot-objects) if plot=True
#     """
#     log = logging.getLogger(pf.LOG)
#     log.info("Starting fine-mapping at region {}".format(ref_geno))
#
#     # align all GWAS, LD reference, and overlapping molecular weights
#     parameters = align_data(gwas, ref_geno, wcollection, ridge=ridge)
#
#     if parameters is None:
#         # break; logging of specific reason should be in align_data
#         return None
#     else:
#         gwas, wmat, meta_data, ldmat = parameters
#
#     zscores = []
#     # run local TWAS
#     for idx, weights in enumerate(wmat.T):
#         log.debug("Computing TWAS association statistic for gene {}".format(meta_data.iloc[idx]["ens_gene_id"]))
#         beta, se = assoc_test(weights, gwas, ldmat, heterogeneity)
#         zscores.append(beta / se)
#
#     # perform fine-mapping
#     zscores = np.array(zscores)
#
#     log.debug("Estimating local TWAS correlation structure")
#     wcor, swld = estimate_cor(wmat, ldmat, intercept)
#
#     m = len(zscores)
#     rm = range(m)
#     pips = np.zeros(m)
#
#     if intercept:
#         # should really be done at the SNP level first ala Barfield et al 2018
#         log.debug("Regressing out average tagged pleiotropic associations")
#         zscores, inter_z = get_resid(zscores, swld, wcor)
#     else:
#         inter_z = None
#
#     k = m if max_genes > m else max_genes
#     null_res = m * np.log(1 - prior_prob)
#     marginal = null_res
#     # enumerate all subsets
#     for subset in chain.from_iterable(combinations(rm, n) for n in range(1, k + 1)):
#
#         local = bayes_factor(zscores, subset, wcor, prior_chisq, prior_prob)
#
#         # keep track for marginal likelihood
#         marginal = np.logaddexp(local, marginal)
#
#         # marginalize the posterior for marginal-posterior on causals
#         for idx in subset:
#             if pips[idx] == 0:
#                 pips[idx] = local
#             else:
#                 pips[idx] = np.logaddexp(pips[idx], local)
#
#     pips = np.exp(pips - marginal)
#     null_res = np.exp(null_res - marginal)
#
#     # Query the db to grab model attributes
#     # We might want to filter to only certain attributes at some point
#     session = pf.get_session()
#     attr = pd.read_sql(session.query(pf.ModelAttribute)
#                        .filter(pf.ModelAttribute.model_id.in_(meta_data.model_id.values.astype(object)))  # why doesn't inte64 work!?!
#                        .statement, con=session.connection())
#
#     # convert from long to wide format
#     attr = attr.pivot("model_id", "attr_name", "value")
#
#     # clean up and return results
#     region = str(ref_geno).replace(" ", "")
#
#     # dont sort here to make plotting easier
#     df = create_output(meta_data, attr, zscores, pips, null_res, region)
#
#     log.info("Completed fine-mapping at region {}".format(ref_geno))
#     if plot:
#         log.info("Creating FOCUS plots at region {}".format(ref_geno))
#         plot_arr = pf.focus_plot(wcor, df)
#
#         # sort here and create credible set
#         df = add_credible_set(df, credible_set=credible_level)
#         return df, plot_arr
#
#     else:
#         # sort here and create credible set
#         df = add_credible_set(df, credible_set=credible_level)
#         return df

def calculate_pips(meta_data, wmat, ldmat, max_genes, prior_prob, intercept):
    """
    Calculate Posterior Inclusion Probability (PIPs)

    :param meta_data: pandas.DataFrame Metadata about the gene models (including Z scores)
    :param wmat: numpy.ndarray eQTL weight matrix for a risk region
    :param ldmat: numpy.ndarray LD matrix for a risk region
    :param max_genes: int or None the maximum number of genes to include in any given causal configuration. None if all genes
    :param prior_prob: float prior probability for a gene to be causal
    :param intercept: bool to return the intercept variable or not

    :return: meta data with PIPs, null PIPs, and Expression correlation matrix
    """
    log = logging.getLogger(pf.LOG)
    n_pop = len(meta_data)

    # Calculate prior chisq
    prior_chisq = [None] * n_pop
    for i in range(n_pop):
        wcor, swld = estimate_cor(wmat[i], ldmat[i], intercept = False)
        solution, resid, rank_wcor, svs = lin.lstsq(wcor, np.asarray(meta_data[i]["twas_z"]))
        prior_chisq[i] = np.dot(np.asarray(meta_data[i]["twas_z"]), solution)

    # perform fine-mapping
    log.debug("Estimating local TWAS correlation structure.")
    wcor = [None] * n_pop
    swld = [None] * n_pop
    for i in range(n_pop):
        wcor_tmp, swld_tmp = estimate_cor(wmat[i], ldmat[i], intercept)
        wcor[i] = wcor_tmp
        swld[i] = swld_tmp

    null_res = [None] if n_pop == 1 else [None] * (n_pop + 1)
    m = len(meta_data[0])
    rm = range(m)
    k = m if max_genes > m else max_genes
    for i in range(1 if n_pop == 1 else n_pop+1):
        pips_tmp = np.zeros(m)
        null_res_tmp = m * np.log(1 - prior_prob)
        marginal = null_res_tmp
        # enumerate all subsets
        for subset in chain.from_iterable(combinations(rm, n) for n in range(1, k + 1)):
            local = 0
            # if it's the case, we need to do Multi-ancestry Fine-mapping
            if i == (len(null_res)-1) and i > 0:
                for j in range(n_pop):
                    local_bf, local_prior = bayes_factor(np.asarray(meta_data[j]["twas_z"]), subset, wcor[j], prior_chisq[j], prior_prob)
                    if j == 0:
                        local += local_bf + local_prior
                    else:
                        local += local_bf
            else:
                local_bf, local_prior = bayes_factor(np.asarray(meta_data[i]["twas_z"]), subset, wcor[i], prior_chisq[i], prior_prob)
                local = local_bf + local_prior

            # keep track for marginal likelihood
            marginal = np.logaddexp(local, marginal)

            # marginalize the posterior for marginal-posterior on causals
            for idx in subset:
                if pips_tmp[idx] == 0:
                    pips_tmp[idx] = local
                else:
                    pips_tmp[idx] = np.logaddexp(pips_tmp[idx], local)

        pips_tmp = np.exp(pips_tmp - marginal)
        null_res_tmp = np.exp(null_res_tmp - marginal)

        if i == (len(null_res)-1) and i > 0:
            # store the me pips into first pop's metadata
            meta_data[0].insert(len(meta_data[0].columns), "pips_me", pips_tmp)
        else:
            meta_data[i].insert(len(meta_data[i].columns), "pips", pips_tmp)
        null_res[i] = null_res_tmp

    return meta_data, null_res, wcor

def num_convert(i):
    """
    A hashmap to convert numeric number to 1st, 2nd.

    :param i: int natural number

    :return: str 1st, 2nd, and so on
    """

    nth = {1: "1st",
    2: "2nd",
    3: "3rd",
    4: "4th",
    5: "5th",
    6: "6th"}

    return nth[i]

def rearrange_columns(df, prior_prob):
    """
    Re-arrange FOCUS output table columns for better user experience.

    :param df: pandas.DataFrame containing all the necessary FOCUS columns
    :param prioor_prob: float use this number to back the number of genes in the block 

    :return: pandas.DataFrame containing all the necessary FOCUS columns in order of non-pop-specific parameters, pop-specific non-pips parameters, pop-specific pips parameters, and me-pips parameters.
    """
    n_pop = np.sum(df.columns.str.contains("pips"))
    n_pop = 1 if n_pop == 1 else n_pop - 1
    df["block_genes"] = 1 / prior_prob
    # move non-pop parameters upfront
    for i in range(n_pop):
        tmp = df.columns.str.contains(f"pop{i+1}")
        tmp_name = df.columns.values[~tmp].tolist() + df.columns.values[tmp].tolist()
        df = df[tmp_name]

    # move twas_z, pips, and credible set in the endswith
    for i in range(n_pop):
        tmp1 = df.pop(f"twas_z_pop{i+1}")
        tmp2 = df.pop(f"pips_pop{i+1}")
        tmp3 = df.pop(f"in_cred_set_pop{i+1}")
        df.insert(len(df.columns), f"twas_z_pop{i+1}", tmp1)
        df.insert(len(df.columns), f"pips_pop{i+1}", tmp2)
        df.insert(len(df.columns), f"in_cred_set_pop{i+1}", tmp3)

        # rearrange pips_me to the last rows
    if n_pop > 1:
        tmp1 = df.pop("pips_me")
        tmp2 = df.pop("in_cred_set_me")
        df.insert(len(df.columns), "pips_me", tmp1)
        df.insert(len(df.columns), "in_cred_set_me", tmp2)

    return df

def fine_map(gwas, wcollection, ref_geno, block, intercept=False, heterogeneity=False, max_genes=3, ridge=0.1, prior_prob=1e-3, credible_level=0.9, plot=False, max_impute=0.5, min_r2pred=0.7, tissue_pr_gene = False):
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
    :param credible_level: float the credible-level to compute credible gene sets (default = 0.9)
    :param plot: bool whether or not to generate visualizations/plots at the risk region
    :param min_r2pred: minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores (default = 0.7)
    :param max_impute: maximum fraction of SNPs allowed to be missing per gene, and will be imputed using LD (default = 0.5)
    :param tissue_pr_gene: boolean variable to indicate whether genes are prioritized by tissues.

    :return: pandas.DataFrame containing the TWAS statistics and fine-mapping results if plot=False.
        (pandas.DataFrame, list of plot-objects) if plot=True
    """
    log = logging.getLogger(pf.LOG)
    log.info(f"Fine-mapping starts at region {block}.")

    # align all GWAS, LD reference, and overlapping molecular weights
    n_pop = len(gwas)
    gwas_copy = gwas
    gwas = [None] * n_pop
    wmat = [None] * n_pop
    meta_data = [None] * n_pop
    ldmat = [None] * n_pop

    for i in range(n_pop):

        if n_pop == 1:
            log.info(f"Aligning GWAS, LD, and eQTL weights for the single population. Region {block} will skip if following errors occur.")
        else:
            log.info(f"Aligning GWAS, LD, and eQTL weights for {num_convert(i+1)} population. Region {block} will skip if following errors occur.")


        parameters_tmp = align_data(gwas_copy[i], ref_geno[i], wcollection[i], ridge=ridge, max_impute=max_impute, min_r2pred=min_r2pred)
        if parameters_tmp is None:
            # break; logging of specific reason should be in align_data
            return None
        else:
            gwas[i], wmat[i], meta_data[i], ldmat[i] = parameters_tmp

    # Get common genes (genes-tissue pair) across populations
    # if genes are not prioritized by tissue, use ens_gene_id, tissue, and model_id to be identifier
    # if genes are prioritized by tissue, use ens_gene_id, and model_id to be identifier
    # import pdb; pdb.set_trace()
    gene_identifier = ["ens_gene_id", "model_id", "ref_name"] if tissue_pr_gene else ["ens_gene_id", "model_id", "tissue", "ref_name"]
    gene_set = meta_data[0][gene_identifier]
    for i in range(n_pop - 1):
        gene_set = pd.merge(gene_set, meta_data[i+1][gene_identifier], how = "inner")

    if len(gene_set) == 0:
        log.warning(f"No common genes found at region {block}. Skipping.")
        return
    else:
        log.info(f"Find {len(gene_set)} common genes to be fine-mapped at region {block}.")

    # foo columns to get the index
    gene_set.insert(len(gene_set.columns), "foo", "foo")
    for i in range(n_pop):
        df_tmp = pd.merge(meta_data[i], gene_set, how = "left", on = gene_identifier)
        idx = df_tmp[pd.notna(df_tmp.foo)].index
        wmat[i] = wmat[i].T[idx].T
        meta_data[i] = meta_data[i].iloc[idx, :]

    # run local TWAS
    for i in range(n_pop):

        if n_pop == 1:
            log.info(f"Running TWAS for the single population.")
        else:
            log.info(f"Running TWAS for {num_convert(i + 1)} population.")

        zscores_tmp = []
        for idx, weights in enumerate(wmat[i].T):
            log.debug(f"Computing TWAS association statistic for gene {meta_data[i].iloc[idx]['ens_gene_id']}.")
            beta, se = assoc_test(weights, gwas[i], ldmat[i], heterogeneity)
            zscores_tmp.append(beta / se)
        meta_data[i].insert(len(meta_data[i].columns), "twas_z", zscores_tmp)

    for i in range(n_pop):
        if intercept:
            # should really be done at the SNP level first ala Barfield et al 2018
            log.debug(f"Regressing out average tagged pleiotropic associations for {num_convert(i + 1)} population")
            twas_z_tmp, inter_z_tmp = get_resid(zscores[i], swld[i], wcor[i])
            meta_data[i].insert(len(meta_data[i].columns), "twas_z", twas_z_tmp)
            meta_data[i].insert(len(meta_data[i].columns), "inter_z", inter_z_tmp)
        else:
            meta_data[i].insert(len(meta_data[i].columns), "inter_z", None)

    # calculate pips for single pop, and me.
    log.info(f"Calculating PIPs.")
    meta_data, null_res, wcor = calculate_pips(meta_data, wmat, ldmat, max_genes, prior_prob, intercept)

    # Query the db to grab model attributes
    # We might want to filter to only certain attributes at some point

    session = [None] * n_pop
    attr = [None] * n_pop
    region = [None] * n_pop
    for i in range(n_pop):
        session_tmp = pf.get_session(idx = i)
        attr_tmp = pd.read_sql(session_tmp.query(pf.ModelAttribute)
                        .filter(pf.ModelAttribute.model_id.in_(meta_data[i].model_id.values.astype(object)))  # why doesn't inte64 work!?!
                        .statement, con=session_tmp.connection())
        # convert from long to wide format
        attr_tmp = attr_tmp.pivot("model_id", "attr_name", "value")
        attr[i] = attr_tmp
        # clean up and return results
        region_tmp = str(ref_geno[i]).replace(" ", "")
        region[i] = region_tmp

    # dont sort here to make plotting easier

    df = create_output(meta_data, attr, null_res, region)
    # Output the partition blocks in which focus performs
    df.insert(0, "block", block)

    log.info(f"Completed fine-mapping at region {block}.")
    if plot:
        log.info(f"Creating FOCUS plots at region {block}.")
        plot_arr = pf.focus_plot(wcor, df)

        # sort here and create credible set
        df = add_credible_set(df, credible_set=credible_level)
        df = rearrange_columns(df, prior_prob)
        return df, plot_arr

    else:
        # sort here and create credible set
        df = add_credible_set(df, credible_set=credible_level)
        df = rearrange_columns(df, prior_prob)
        return df
