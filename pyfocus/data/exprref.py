import pandas as pd
import pyfocus as pf
import numpy as np

from pandas_plink import read_plink


class ExprRef(object):

    def __init__(self, geno, sample_info, var_info, pheno=None, covar=None, gene_info=None):
        self._geno = geno                # IID x SNP
        self._sample_info = sample_info  # info on individuals
        self._var_info = var_info        # info on genotype variants

        self._pheno = pheno
        self._covar = covar
        self._gene_info = gene_info

        self._has_aligned = False
        return

    def __iter__(self):
        if not self._has_aligned:
            self.align()

        if self._gene_info is None:
            raise ValueError("Gene Info must be set before iterating through expression reference data")
        elif self._pheno is None:
            raise ValueError("Phenotype must be set before iterating through expression reference data")

        # right now the covariates are a pandas data-frame
        # since data are aligned convert into matrix suitable for regression
        X = self._prepare_covar()

        bound = 5e5
        for idx, row in self._gene_info.iterrows():
            tx = row.txid
            chrom = row.chrom
            txstart = row.txstart

            # find overlapping SNPs
            snps = self._var_info.loc[(self._var_info.chrom == chrom) & (np.abs(self._var_info.pos - txstart) <= bound)]

            # pull genotype data
            G = self._geno[snps.i.values, :].compute().T

            # impute missing with column mean
            col_mean = np.nanmean(G, axis=0)
            inds = np.where(np.isnan(G))
            G[inds] = np.take(col_mean, inds[1])

            # center and scale
            G -= np.mean(G, axis=0)
            G /= np.std(G, axis=0)

            # inverse-normal transform the expression data
            y = pf.inv_norm(self._pheno[tx].values)

            yield y, X, G, snps, row

        return

    def parse_pheno(self, path):
        self._pheno = pd.read_table(path, delim_whitespace=True)
        return

    def parse_covar(self, path):
        self._covar = pd.read_table(path, delim_whitespace=True)
        return

    def _prepare_covar(self):
        # prepare the covariates data-frame into a matrix usable for regression
        # assumes that the data have already been aligned
        sub_df = self._covar  # we need to first remove FID/IID if they exist
        sub_df = sub_df.loc[:, (sub_df.columns != "FID") & (sub_df.columns != "IID")]

        # one-hot encoding of categorical variables
        covar_mat = pd.get_dummies(sub_df, drop_first=True).values

        return covar_mat

    def parse_gene_info(self, path):
        # TODO: validate meta-data. Must have columns we need for training
        self._gene_info = pd.read_table(path, delim_whitespace=True)
        self._gene_info.chrom = self._gene_info.chrom.astype(str)
        return

    def align(self):
        # TODO: implement data-alignment
        self._has_aligned = True
        return

    @classmethod
    def from_plink(cls, path):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', 'FutureWarning')
            bim, fam, bed = read_plink(path, verbose=False)
            bim.chrom = bim.chrom.astype(str)

            return cls(bed, fam, bim)

    # @classmethod
    # def from_vcf(cls, path):
    #     pass
