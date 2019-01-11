import itertools as it
import logging

import pyfocus

__all__ = ["flip_values", "clean_chrom", "inv_norm"]

# Base-handling code is from LDSC...
# complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# bases
BASES = COMPLEMENT.keys()
# true iff strand ambiguous
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                    for x in it.product(BASES, BASES)
                    if x[0] != x[1]}

# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {x for x in map(lambda y: ''.join(y), it.product(BASES, BASES))
              if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}

# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {x for x in map(lambda y: ''.join(y), it.product(VALID_SNPS, VALID_SNPS))
                 # strand and ref match
                 if ((x[0] == x[2]) and (x[1] == x[3])) or
                 # ref match, strand flip
                 ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                 # ref flip, strand match
                 ((x[0] == x[3]) and (x[1] == x[2])) or
                 ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip

# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
                ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                # strand flip
                ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                for x in MATCH_ALLELES}


def flip_values(values, data_a1, data_a2, ref_a1, ref_a2):
    alleles = data_a1 + data_a2 + ref_a1 + ref_a2
    flip_flags = alleles.apply(lambda y: FLIP_ALLELES[y])
    values *= (-1) ** flip_flags

    return values


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
