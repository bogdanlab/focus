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
    from numpy import clip, inf

    epsilon = np.finfo(float).tiny

    lrs = clip(2 * (logl_H1 - logl_H0), epsilon, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, epsilon, 1 - epsilon)


def _fit_cis_herit(y, K_cis, X, compute_lrt=True):
    from glimix_core.lmm import LMM
    from numpy_sugar.linalg import economic_qs_linear
    from scipy.stats import norm
    from scipy.linalg import lstsq

    log = logging.getLogger(pyfocus.LOG)

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
        fixed_betas_0, sosqs, ranks, svals = lstsq(X, y)
        s2e = sosqs / len(y)  # LMM also uses MLE estimation, so don't adjust for bias right now

        logl_0 = np.sum(norm.logpdf(y, loc=np.dot(X, fixed_betas_0), scale=np.sqrt(s2e)))
        pval = _lrt_pvalue(logl_0, logl_1)
        log.debug("Estimated cis-h2g = {} (P = {})".format(cis_scale / (cis_scale + noise_scale + fe_scale), pval))
    else:
        pval = None
        log.debug("Estimated cis-h2g = {}".format(cis_scale / (cis_scale + noise_scale + fe_scale)))

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
    fe_var, s2u, s2e, logl, fixed_betas, pval = _fit_cis_herit(y, K_cis, X)
    yresid = y - np.dot(X, fixed_betas)

    attrs["h2g"] = s2u / (fe_var + s2u + s2e)
    attrs["h2g.logl"] = logl
    attrs["h2g.pvalue"] = pval

    # Total variance
    n, p = Z.shape
    if n < p:
        V = s2u * K_cis + s2e * np.ones(len(y))
        Vinv = np.linalg.pinv(V)

        # estimate BLUP for genomic prediction
        betas = mdot([Z.T * s2u, Vinv, (y - np.dot(X, fixed_betas))])
        if include_ses:
            # I'm not sure how to pull out marker effects in fast leave-one-out using classic GBLUP, yet
            ses = None
        else:
            ses = None

    else:
        # ridge solution (i.e. rrBLUP) should be faster when p < n
        ZtZpDinv = np.linalg.pinv(np.dot(Z.T, Z) + np.eye(p) * (s2e / s2u))
        betas = mdot([ZtZpDinv, Z.T, yresid])
        if include_ses:
            # jack-knife standard-errors over the fast leave-one-out estimates using rrBLUP
            h = np.array([mdot([Z[i], ZtZpDinv, Z[i]]) for i in range(n)])
            e = yresid - np.dot(Z, betas)
            beta_jk = [betas - np.dot(ZtZpDinv, Z[i] * e[i] / (1 - h[i])) for i in range(n)]
            ses = np.sqrt(np.mean(beta_jk, axis=0) * (n - 1))
        else:
            ses = None

    return betas, ses, attrs


def _train_horseshoe(y, Z, X, include_ses=False):
    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing Horseshoe model")

    # need to install STAN to do this
    betas = None
    attrs = dict()
    return betas, attrs


def train_model(y, X, G, method="GBLUP", include_ses=False):

    if method == LASSO:
        betas, ses, attrs = _train_lasso(y, G, X, include_ses)
    elif method == GBLUP:
        betas, ses, attrs = _train_gblup(y, G, X, include_ses)
    elif method == HORSESHOE:
        betas, ses, attrs = _train_gblup(y, G, X, include_ses)
    else:
        raise ValueError("Unknown inference method {}".format(method))

    return betas, ses, attrs
