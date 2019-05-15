import logging

import numpy as np
import pyfocus


__all__ = ["train_model", "METHODS"]

LASSO = "LASSO"
GBLUP = "GBLUP"
ENET = "ENET"
METHODS = [LASSO, GBLUP, ENET]


def _lrt_pvalue(logl_H0, logl_H1, dof=1):
    from scipy.stats import chi2
    from numpy import clip, inf

    epsilon = np.finfo(float).tiny

    lrs = clip(2 * (logl_H1 - logl_H0), epsilon, inf)
    pv = chi2(df=dof).sf(lrs)
    return clip(pv, epsilon, 1 - epsilon)


def _fit_cis_herit(y, K_cis, X=None, compute_lrt=True):
    log = logging.getLogger(pyfocus.LOG)

    try:
        from glimix_core.lmm import LMM
        from numpy_sugar.linalg import economic_qs_linear
    except ImportError as ie:
        log.error("Training submodule requires glimix-core>=2.0.0 and numpy-sugar to be installed.")
        raise

    from scipy.stats import norm
    from scipy.linalg import lstsq

    if X is None:
        X = np.ones((len(y), 1))

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


def _train_lasso(y, Z, X, include_ses=False, p_threshold=0.01):
    log = logging.getLogger(pyfocus.LOG)
    try:
        from limix.qc import normalise_covariance
        from sklearn.linear_model import Lasso
    except ImportError as ie:
        log.error("Training submodule requires limix>=2.0.0 and sklearn to be installed.")
        raise
    from scipy.linalg import lstsq

    log.debug("Initializing LASSO model")

    n = len(y)
    attrs = dict()

    K_cis = np.dot(Z, Z.T)
    K_cis = normalise_covariance(K_cis)
    fe_var, s2u, s2e, logl, fixed_betas, pval = _fit_cis_herit(y, K_cis, X)
    if pval > p_threshold:
        log.info("h2g pvalue {} greater than threshold {}. Skipping".format(pval, p_threshold))
        return None

    h2g = s2u / (s2u + s2e + fe_var)

    attrs["h2g"] = h2g
    attrs["h2g.logl"] = logl
    attrs["h2g.pvalue"] = pval

    # we only want to penalize SNP effects and not covariate effects...
    fixed_betas, sum_resid, ranks, svals = lstsq(X, y)
    yresid = y - np.dot(X, fixed_betas)

    # PLINK-style LASSO
    lambda_max = np.linalg.norm(Z.T.dot(yresid), np.inf) / float(n)

    def _gen_e():
        e = np.random.normal(size=n)
        return np.linalg.norm(Z.T.dot(e), np.inf)

    min_tmp = np.median([_gen_e() for _ in range(1000)])
    sige = np.sqrt(1.0 - h2g + (1.0 / float(n)))
    lambda_min = (sige / n) * min_tmp

    # 100 values spaced logarithmically from lambda-min to lambda-max
    alphas = np.exp(np.linspace(np.log(lambda_min), np.log(lambda_max), 100))

    # fit LASSO solution using coordinate descent, updating with consecutively smaller penalties
    lasso = Lasso(fit_intercept=True, warm_start=True)
    for penalty in reversed(alphas):
        lasso.set_params(alpha=penalty)
        lasso.fit(Z, yresid)

    betas = lasso.coef_

    attrs["r2"] = lasso.score(Z, yresid)
    attrs["resid.var"] = sum((yresid - lasso.predict(Z)) ** 2) / (n - 1)

    if include_ses:
        # TODO: bootstrap?
        ses = None
    else:
        ses = None

    return betas, ses, attrs


def _train_enet(y, Z, X, include_ses=False, p_threshold=0.01):
    log = logging.getLogger(pyfocus.LOG)
    try:
        from limix.qc import normalise_covariance
        from sklearn.linear_model import ElasticNetCV
    except ImportError as ie:
        log.error("Training submodule requires limix>=2.0.0 and sklearn to be installed.")
        raise
    from scipy.linalg import lstsq

    log.debug("Initializing ElasticNet model")

    n = len(y)
    attrs = dict()

    K_cis = np.dot(Z, Z.T)
    K_cis = normalise_covariance(K_cis)
    fe_var, s2u, s2e, logl, fixed_betas, pval = _fit_cis_herit(y, K_cis, X)
    if pval > p_threshold:
        log.info("h2g pvalue {} greater than threshold {}. Skipping".format(pval, p_threshold))
        return None

    h2g = s2u / (s2u + s2e + fe_var)

    attrs["h2g"] = h2g
    attrs["h2g.logl"] = logl
    attrs["h2g.pvalue"] = pval

    # we only want to penalize SNP effects and not covariate effects...
    fixed_betas, sum_resid, ranks, svals = lstsq(X, y)
    yresid = y - np.dot(X, fixed_betas)

    enet = ElasticNetCV(l1_ratio=0.5, fit_intercept=True, cv=5)
    enet.fit(Z, yresid)
    betas = enet.coef_

    attrs["r2"] = enet.score(Z, yresid)
    attrs["resid.var"] = sum((yresid - enet.predict(Z)) ** 2) / (n - 1)

    if include_ses:
        # TODO: bootstrap?
        ses = None
    else:
        ses = None

    return betas, ses, attrs


def _train_gblup(y, Z, X, include_ses=False, p_threshold=0.01):
    log = logging.getLogger(pyfocus.LOG)

    try:
        from limix.qc import normalise_covariance
    except ImportError as ie:
        log.error("Training submodule requires limix>=2.0.0 and sklearn to be installed.")
        raise
    from numpy.linalg import multi_dot as mdot
    from scipy.linalg import pinvh

    log.debug("Initializing GBLUP model")

    attrs = dict()

    # estimate heritability using limix
    K_cis = np.dot(Z, Z.T)
    K_cis = normalise_covariance(K_cis)
    fe_var, s2u, s2e, logl, fixed_betas, pval = _fit_cis_herit(y, K_cis, X)
    yresid = y - np.dot(X, fixed_betas)

    if pval > p_threshold:
        log.info("h2g pvalue {} greater than threshold {}. Skipping".format(pval, p_threshold))
        return None

    attrs["h2g"] = s2u / (fe_var + s2u + s2e)
    attrs["h2g.logl"] = logl
    attrs["h2g.pvalue"] = pval

    # Total variance
    n, p = Z.shape

    # ridge solution (i.e. rrBLUP)
    # this will be slower than normal GBLUP when p > n but is a little bit more flexible
    ZtZpDinv = pinvh(np.dot(Z.T, Z) + np.eye(p) * (s2e / s2u))
    betas = mdot([ZtZpDinv, Z.T, yresid])

    if include_ses:
        # TODO: come back to this with matrix operations rather than list comprehensions
        # jack-knife standard-errors over the fast leave-one-out estimates using rrBLUP
        """
        h = np.array([mdot([Z[i], ZtZpDinv, Z[i]]) for i in range(n)])
        e = yresid - np.dot(Z, betas)
        beta_jk = [betas - np.dot(ZtZpDinv, Z[i] * e[i]) / (1 - h[i]) for i in range(n)]
        ses = np.sqrt(np.mean(beta_jk, axis=0) * (n - 1))
        """
        ses = None
    else:
        ses = None

    return betas, ses, attrs


def train_model(y, X, G, method="GBLUP", include_ses=False, p_threshold=0.01):

    if method == LASSO:
        res = _train_lasso(y, G, X, include_ses, p_threshold)
    elif method == GBLUP:
        res = _train_gblup(y, G, X, include_ses, p_threshold)
    elif method == ENET:
        res = _train_enet(y, G, X, include_ses, p_threshold)
    else:
        raise ValueError("Unknown inference method {}".format(method))

    return res
