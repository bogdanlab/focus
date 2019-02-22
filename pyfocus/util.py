
__all__ = ["inv_norm"]


def inv_norm(pheno):
    import scipy.stats as stats

    if pheno is None:
        raise ValueError("Expected non-null numpy vector")
    if pheno.ndim > 1:
        raise ValueError("Expected numpy vector")

    k = 3 / 8.0
    n = len(pheno)
    ranks = stats.rankdata(pheno)

    return stats.norm.isf((ranks - k) / (n - 2 * k + 1))
