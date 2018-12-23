
__all__ = ["clean_chrom", "inv_norm"]


def clean_chrom(df, chrom, column_name):
    import pandas as pd

    if pd.api.types.is_string_dtype(df[column_name]) or pd.api.types.is_categorical_dtype(df[column_name]):
        ret_val = str(chrom)
    else:
        ret_val = int(chrom)

    return ret_val


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
