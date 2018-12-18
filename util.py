import pandas as pd

__all__ = ["clean_chrom"]


def clean_chrom(df, chrom, column_name):
    if pd.api.types.is_string_dtype(df[column_name]) or pd.api.types.is_categorical_dtype(df[column_name]):
        ret_val = str(chrom)
    else:
        ret_val = int(chrom)

    return ret_val

