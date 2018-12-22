import logging

import numpy as np
import pyfocus


__all__ = ["train_model", "METHODS"]

LASSO = "LASSO"
GBLUP = "GBLUP"
HORSESHOE = "HORSESHOE"
METHODS = [LASSO, GBLUP, HORSESHOE]


def _lrt_pvalue(logl_H0, logl_H1, dof=1):
    from scipy.stats import chi2
    from numpy import asarray, clip, inf

    epsilon = np.finfo(float).tiny

    lrs = clip(-2 * logl_H0 + 2 * asarray(logl_H1), epsilon, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, epsilon, 1 - epsilon)


def _fit_cis_herit(y, K_cis, X, compute_lrt=False):
    from glimix_core.lmm import LMM
    from numpy_sugar.linalg import economic_qs_linear
    from scipy.stats import norm
    from scipy.linalg import lstsq

    log = logging.getLogger(pyfocus.LOG)
    log.info("Computing cis-heritability")

    K_cis = economic_qs_linear(K_cis)
    lmm = LMM(y, X, K_cis)
    lmm.fit(verbose=False)

    fixed_betas = lmm.beta
    logl_1 = lmm.lml()

    cis_scale = lmm.v0
    noise_scale = lmm.v1
    fe_scale = lmm.fixed_effects_variance

    if compute_lrt:
        n, p = X.shape
        # reduced model is just OLS regression for fixed-effects
        fixed_betas_0, yresid, ranks, svals = lstsq(X, y)
        s2e = np.mean(yresid ** 2)

        logl_0 = np.sum(norm.logpdf(loc=np.dot(X, fixed_betas_0), scale=np.sqrt(s2e)))
        pval = _lrt_pvalue(logl_0, logl_1)
    else:
        pval = None

    return fe_scale, cis_scale, noise_scale, logl_1, fixed_betas, pval


def _train_lasso(y, Z, X, include_ses=False):
    from sklearn.linear_model import Lasso
    from scipy.linalg import lstsq

    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing LASSO model")

    # we only want to penalize SNP effects and not covariate effects...
    fixed_betas, yresid, ranks, svals = lstsq(X, y)

    # we need to estimate the h2g for tuning parameter
    lasso = Lasso()
    lasso.fit(X, y)
    betas = lasso.coef_

    attrs = dict()
    return betas, attrs


def _train_gblup(y, Z, X, include_ses=False):
    import limix
    from numpy.linalg import multi_dot as mdot

    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing GBLUP model")

    attrs = dict()

    # estimate heritability using LIMIX
    K_cis = np.dot(Z, Z.T)
    K_cis = limix.qc.normalise_covariance(K_cis)
    fe_var, s2u, s2e, logl, fixed_betas, pval = _fit_cis_herit(y, K_cis, X, compute_lrt=True)

    attrs["h2g"] = s2u / (fe_var + s2u + s2e)
    attrs["h2g.pvalue"] = pval

    # Total variance
    V = s2u * K_cis + s2e * np.ones(len(y))
    Vinv = np.linalg.inv(V)

    # estimate BLUP for genomic prediction
    betas = mdot([Z * s2u, Vinv, (y - np.dot(X, fixed_betas))])

    return betas, attrs


def _train_horseshoe(y, Z, X, include_ses=False):
    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing Horseshoe model")

    # need to install STAN to do this
    betas = None
    attrs = dict()
    return betas, attrs


def train_model(y, X, G, method="GBLUP", include_ses=False):
    log = logging.getLogger(pyfocus.LOG)
    log.info("Initializing model inference")

    if method == LASSO:
        betas, attrs = _train_lasso(y, G, X, include_ses)
    elif method == GBLUP:
        betas, attrs = _train_gblup(y, G, X, include_ses)
    elif method == HORSESHOE:
        betas, attrs = _train_gblup(y, G, X, include_ses)
    else:
        raise ValueError("Unknown inference method {}".format(method))

    return betas, attrs
