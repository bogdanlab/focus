import logging
import warnings

import pyfocus
import numpy as np
import pandas as pd
import pkg_resources

import scipy.linalg as lin

from pandas_plink import read_plink


class IndBlocks(object):
    """
    Class to wrap/iterate approximately independent regions
    """

    CHRCOL = "CHR"
    STARTCOL = "START"
    STOPCOL = "STOP"

    REQ_COLS = [CHRCOL, STARTCOL, STOPCOL]

    def __init__(self, regions=None):
        if regions is None:
            local_ld_blocks = pkg_resources.resource_filename(__name__, 'ld_blocks/grch37.eur.loci.bed')
            self._regions = pd.read_csv(local_ld_blocks)
        else:
            # type checking in python == dumb
            if type(regions) is str or type(regions) is file:
                self._regions = pd.read_csv(regions)
            elif type(regions) is pd.core.frame.DataFrame:
                self._regions = regions

            for column in IndBlocks.REQ_COLS:
                if column not in self._regions:
                    raise ValueError("{}-column not found in regions".format(column))

            self._regions = self._regions[IndBlocks.REQ_COLS]

        return

    def subset_by_pos(self, chrom, start, stop):

        if chrom is None:
            raise ValueError("chrom argument cannot be `None` in subset_by_pos")

        chrom = pyfocus.clean_chrom(self._regions, chrom, IndBlocks.CHRCOL)
        df = self._regions.loc[self._regions[IndBlocks.CHRCOL] == chrom]

        if stop is None:
            stop = max(df[IndBlocks.STOPCOL])
        if start is None:
            start = min(df[IndBlocks.STARTCOL])

        # grab the intervals that overap the provided start and stop positions on same chromosome
        locs = df.loc[df.apply(lambda x: min(x[IndBlocks.STOPCOL], stop) - max(x[IndBlocks.STARTCOL], start)) > 0]

        return IndBlocks(locs)

    def __iter__(self):
        for row in self._regions.itertuples(name=None):
            yield row

        return


class RefPanel(object):
    CHRCOL = "chrom"
    SNPCOL = "snp"
    BPCOL = "pos"
    A1COL = "a1"
    A2COL = "a0"

    # chrom snp  cm pos a0 a1 i

    def __init__(self, snp_info, sample_info, geno):
        self._snp_info = snp_info
        if pd.api.types.is_categorical_dtype(self._snp_info[RefPanel.A1COL]):
            self._snp_info.loc[:, RefPanel.A1COL] = self._snp_info[RefPanel.A1COL].astype('str')
            self._snp_info.loc[:, RefPanel.A2COL] = self._snp_info[RefPanel.A2COL].astype('str')
        self._sample_info = sample_info
        self._geno = geno
        return

    def __len__(self):
        return len(self._snp_info)

    def __str__(self):
        start_chr = self._snp_info[RefPanel.CHRCOL].iloc[0]
        stop_chr = self._snp_info[RefPanel.CHRCOL].iloc[-1]

        start_bp = self._snp_info[RefPanel.BPCOL].iloc[0]
        stop_bp = self._snp_info[RefPanel.BPCOL].iloc[-1]
        return "{}:{} - {}:{}".format(start_chr, int(start_bp), stop_chr, int(stop_bp))

    def get_partitions(self, window_size, chrom=None, start=None, stop=None):
        """
        Lazily iterate over location partitions
        """
        log = logging.getLogger(pyfocus.LOG)

        chroms = self._snp_info[RefPanel.CHRCOL].unique()

        if chrom is not None:
            chrom = pyfocus.clean_chrom(self._snp_info, chrom, RefPanel.CHRCOL)
            if chrom not in chroms:
                msg = "User supplied chromosome {} is not found in data".format(chrom)
                log.error(msg)
                return

        for chrm in chroms:
            if chrom is not None and chrom != chrm:
                continue

            snps = self._snp_info.loc[self._snp_info[RefPanel.CHRCOL] == chrm]

            min_pos_indata = snps[RefPanel.BPCOL].min()
            max_pos_indata = snps[RefPanel.BPCOL].max()

            # check against user arguments
            if start is not None:
                min_pos = int(start)
                if min_pos < min_pos_indata and min_pos < max_pos_indata:
                    msg = "User supplied start {} is less than min start found in data {}. Switching to data start"
                    msg = msg.format(min_pos, min_pos_indata)
                    log.warning(msg)
                    min_pos = min_pos_indata
            else:
                min_pos = min_pos_indata

            if stop is not None:
                max_pos = int(stop)
                if max_pos > max_pos_indata and max_pos > min_pos_indata:
                    msg = "User supplied stop {} is greater than max stop found in data {}. Switching to data stop"
                    msg = msg.format(max_pos, max_pos_indata)
                    log.warning(msg)
                    max_pos = max_pos_indata
            else:
                max_pos = max_pos_indata

            nwin = int(np.ceil((max_pos - min_pos + 1) / window_size))
            yield [chrm, min_pos, min(min_pos + window_size, max_pos)]

            last_stop = min_pos + window_size
            for i in range(1, nwin):
                start = last_stop + 1
                stop = min(start + window_size, max_pos)
                yield [chrm, start, stop]
                last_stop = stop

        return

    def subset_by_pos(self, chrom, start=None, stop=None, filter_ambig=True):
        ambig = ["AT", "TA", "CG", "GC"]
        df = self._snp_info

        chrom = pyfocus.clean_chrom(df, chrom, RefPanel.CHRCOL)
        if start is not None and stop is not None:
            snps = df.loc[(df[RefPanel.CHRCOL] == chrom) & (df[RefPanel.BPCOL] >= start) & (df[RefPanel.BPCOL] <= stop)]
        elif start is not None and stop is None:
            snps = df.loc[(df[RefPanel.CHRCOL] == chrom) & (df[RefPanel.BPCOL] >= start)]
        elif start is None and stop is not None:
            snps = df.loc[(df[RefPanel.CHRCOL] == chrom) & (df[RefPanel.BPCOL] <= stop)]
        else:
            snps = df.loc[(df[RefPanel.CHRCOL] == chrom)]

        if filter_ambig:
            alleles = snps[RefPanel.A1COL] + snps[RefPanel.A2COL]
            non_ambig = alleles.apply(lambda y: y.upper() not in ambig)
            snps = snps[non_ambig]

        return RefPanel(snps, self._sample_info, self._geno)

    def overlap_gwas(self, gwas):
        df = self._snp_info
        merged_snps = pd.merge(gwas, df, how="outer", left_on=pyfocus.GWAS.SNPCOL, right_on=pyfocus.RefPanel.SNPCOL)
        merged_snps.drop_duplicates(subset=pyfocus.RefPanel.SNPCOL, inplace=True)
        return merged_snps

    def get_geno(self, snps=None):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            if snps is None:
                return self._geno.compute().T
            else:
                return self._geno[snps.i.values, :].compute().T

    @property
    def sample_size(self):
        return len(self._sample_info)

    def estimate_ld(self, snps=None, adjust=0.1, return_eigvals=False):
        G = self.get_geno(snps)
        n, p = G.shape
        col_mean = np.nanmean(G, axis=0)

        # impute missing with column mean
        inds = np.where(np.isnan(G))
        G[inds] = np.take(col_mean, inds[1])

        if return_eigvals:
            G = (G - np.mean(G, axis=0)) / np.std(G, axis=0)
            _, S, V = lin.svd(G, full_matrices=True)

            # adjust 
            D = np.full(p, adjust)
            D[:len(S)] = D[:len(S)] + (S**2 / n)

            return (np.dot(V.T * D, V), D)
        else:
            return np.corrcoef(G.T) + np.eye(p) * adjust

    @classmethod
    def parse_plink(cls, path):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', 'FutureWarning')
            bim, fam, bed = read_plink(path, verbose=False)
        return RefPanel(bim, fam, bed)
