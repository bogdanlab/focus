import os

__all__ = ["get_compression", "write_output"]


def get_compression(fh):
    # This function from LDSC regression
    # (c) 2014 Brendan Bulik-Sullivan and Hilary Finucane
    """Which sort of compression should we use with read_csv?"""
    if isinstance(fh, file):
        _, ext = os.path.splitext(fh.name)
    elif isinstance(fh, str):
        _, ext = os.path.splitext(fh)
    else:
        raise ValueError("get_compression: argument must be file handle or path")

    if ext.endswith('gz'):
        compression = 'gzip'
    elif ext.endswith('bz2'):
        compression = 'bz2'
    else:
        compression = None

    return compression


def write_output(imputed_gwas, output, append=False):
    imputed_gwas.to_csv(output, sep="\t", mode="a" if append else "w", header=not append, index=False)
    return
