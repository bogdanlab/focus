import logging

import numpy as np
import pyfocus


__all__ = ["train_model"]

LASSO = "LASSO"
GBLUP = "GBLUP"
HORSESHOE = "HORSESHOE"
METHODS = [LASSO, GBLUP, HORSESHOE]


def _fit_cis_herit(y, K_cis, X):
    import limix

    log = logging.getLogger(pyfocus.LOG)
    log.info("Computing cis-heritability")

    glmm = limix.glmm.GLMMComposer(len(y))
    glmm.y = y
    glmm.fixed_effects.append(X)
    glmm.fixed_effects.append_offset()
    glmm.covariance_matrices.append(K_cis)
    glmm.covariance_matrices.append_iid_noise()
    glmm.fit(verbose=False)

    cis_scale = glmm.covariance_matrices[0].scale
    noise_scale = glmm.covariance_matrices[1].scale
    logl = glmm.lml()

    return cis_scale, noise_scale, logl


def _train_lasso(y, Z, X, include_ses=False):
    from sklearn.linear_model import Lasso, LinearRegression
    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing LASSO model")

    # we only want to penalize SNP effects and not covariate effects...
    ols = LinearRegression()
    ols.fit(X, y)
    yresid = y - ols.predict(Z)

    # we need to estimate the h2g for tuning parameter
    lasso = Lasso()
    lasso.fit(X, y)
    betas = lasso.coef_

    attrs = dict()
    return betas, attrs


def _train_gblup(y, Z, X, include_ses=False):
    import limix
    from sklearn.linear_model import LinearRegression
    from numpy.linalg import multi_dot as mdot
    from scipy.linalg import cholesky

    log = logging.getLogger(pyfocus.LOG)
    log.debug("Initializing GBLUP model")

    attrs = dict()

    # estimate heritability using LIMIX
    K_cis = np.dot(Z, Z.T)
    K_cis = limix.qc.normalise_covariance(K_cis) # do we need this?
    s2u, s2e, logl = _fit_cis_herit(y, K_cis, X)
    V = s2u * K_cis + s2e * np.ones(len(y))
    Vinv = np.linalg.inv(V)

    # estimate BLUE estimates using estimated variance parameters
    L = cholesky(Vinv, lower=True)
    Xstar = np.dot(X, L)
    ystar = np.dot(L, y)
    ols = LinearRegression()
    ols.fit(Xstar, ystar)
    fixed_betas = ols.coef_

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


def train_model(y, Z, X, method="LASSO", include_ses=False):
    log = logging.getLogger(pyfocus.LOG)
    log.info("Initializing model inference")

    if method == LASSO:
        betas, attrs = _train_lasso(y, Z, X, include_ses)
    elif method == GBLUP:
        betas, attrs = _train_gblup(y, Z, X, include_ses)
    elif method == HORSESHOE:
        betas, attrs = _train_gblup(y, Z, X, include_ses)
    else:
        raise ValueError("Unknown inference method {}".format(method))

    return betas, attrs
