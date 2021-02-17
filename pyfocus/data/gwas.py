import pandas as pd
import scipy.stats as stats

import pyfocus

__all__ = ["GWAS", "GWASSeries"]


class GWASSeries(pd.Series):
    @property
    def _constructor(self):
        return pyfocus.GWASSeries

    @property
    def _constructor_expanddim(self):
        return pyfocus.GWAS


class GWAS(pd.DataFrame):
    """
    Thin wrapper for a pandas DataFrame object containing GWAS summary data.
    Assumes the GWAS data have gone through LDSC munge-sumstat
    """

    CHRCOL = "CHR"
    SNPCOL = "SNP"
    BPCOL = "BP"
    A1COL = "A1"
    A2COL = "A2"
    ZCOL = "Z"

    PCOL = "P"
    NCOL = "N"

    REQ_COLS = [CHRCOL, SNPCOL, BPCOL, A1COL, A2COL, ZCOL]

    def __init__(self, *args, **kwargs):
        super(GWAS, self).__init__(*args, **kwargs)
        return

    @property
    def _constructor(self):
        return pyfocus.GWAS

    @property
    def _constructor_sliced(self):
        return pyfocus.GWASSeries

    def subset_by_pos(self, chrom, start=None, stop=None):

        if start is not None and stop is not None:
            snps = self.loc[(self[GWAS.CHRCOL] == chrom) & (self[GWAS.BPCOL] >= start) & (self[GWAS.BPCOL] <= stop)]
        elif start is not None and stop is None:
            snps = self.loc[(self[GWAS.CHRCOL] == chrom) & (self[GWAS.BPCOL] >= start)]
        elif start is None and stop is not None:
            snps = self.loc[(self[GWAS.CHRCOL] == chrom) & (self[GWAS.BPCOL] <= stop)]
        else:
            snps = self.loc[(self[GWAS.CHRCOL] == chrom)]

        return GWAS(snps)


    @classmethod
    def parse_gwas(cls, stream):
        dtype_dict = {'CHR': "category", 'SNP': str, 'Z': float, 'N': float, 'A1': str, 'A2': str}
        try:
            df = pd.read_csv(stream, dtype=dtype_dict, delim_whitespace=True, compression='infer')
        except Exception as e:
            raise Exception(f"Parsing GWAS failed for population at {stream}." + str(e))

        for column in GWAS.REQ_COLS:
            if column not in df:
                raise ValueError(f"{column}-column not found in summary statistics at {stream}.")

        if GWAS.PCOL not in df:
            df[GWAS.PCOL] = stats.chi2.sf(df[GWAS.ZCOL].values ** 2, 1)

        return cls(df)
