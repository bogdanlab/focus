import numpy as np
#import cyvcf2

from pandas_plink import read_plink


class ExprRef(object):

    def __init__(self, G=None, covariates=None, pheno_map=None):
        self._G = G             # IID x SNP
        self._X = covariates    # IID COVAR_1 COVAR_2 ... COVAR_k
        self._Y = pheno_map     # IID PHENO_1 PHENO_2 ... PHENO_m
        return

    @property
    def covariates(self):
        return self._X

    @covariates.setter
    def covariates(self, value):
        print("setter of x called")
        self._x = value

    @classmethod
    def from_plink(cls, path):
        with np.warnings.catch_warnings():
            np.warnings.filterwarnings('ignore', 'FutureWarning')
            bim, fam, bed = read_plink(path, verbose=False)

            # 0. Family ID ('FID')
            # 1. Within-family ID ('IID'; cannot be '0')
            # 2. Within-family ID of father ('0' if father isn't in dataset)
            # 3. Within-family ID of mother ('0' if mother isn't in dataset)
            # 4. Sex code ('1' = male, '2' = female, '0' = unknown)
            # 5. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            sample_info = fam.iloc[:, [0, 5]]

            return cls(bed, None, sample_info)

    @classmethod
    def from_vcf(cls, path):
        pass
