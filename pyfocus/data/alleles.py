import itertools as it
import logging

import pyfocus

__all__ = ["check_valid_snp", "check_valid_alleles", "flip_alleles"]

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


def check_valid_snp(a1, a2):
    """
    Checks that the alleles a1, a2 are unambiguous.

    :param a1: pandas column for the A1 allele
    :param a2: pandas column for the A2 allele

    :return: pandas boolean column indicating if SNPs are unambiguous or not
    """
    a = a1 + a2
    return a.isin(VALID_SNPS)


def check_valid_alleles(gwas_a1, gwas_a2, ref_a1, ref_a2):
    """
    Check that A1/A2 alleles in GWAS match A1/A2 alleles in LDRefPanel.
    Allows for A1/A2 flip between GWAS/LDRefPanel.

    :param gwas_a1: pandas column for the A1 allele in GWAS data
    :param gwas_a2: pandas column for the A2 allele in GWAS data
    :param ref_a1: pandas column for the A1 allele in LDRefPanel data
    :param ref_a2: pandas column for the A2 allele in LDRefPanel data

    :return: pandas boolean column indicating if alleles match or not
    """
    alleles = gwas_a1 + gwas_a2 + ref_a1 + ref_a2
    return alleles.apply(lambda y: y in MATCH_ALLELES)


def flip_alleles(zscores, gwas_a1, gwas_a2, ref_a1, ref_a2):
    """
    Flips zscores in the GWAS data so that the sign matches with the reference allele in LDRefPanel genotypes.

    :param zscores: numpy.ndarray of zscores
    :param gwas_a1: pandas column for the A1 allele in GWAS data
    :param gwas_a2: pandas column for the A2 allele in GWAS data
    :param ref_a1: pandas column for the A1 allele in LDRefPanel data
    :param ref_a2: pandas column for the A2 allele in LDRefPanel data

    :return: numpy.ndarray of sign-aligned zscores

    :raises: ValueError if incompatible alleles are present
    """
    log = logging.getLogger(pyfocus.LOG)

    alleles = gwas_a1 + gwas_a2 + ref_a1 + ref_a2

    try:
        flip_flags = alleles.apply(lambda y: FLIP_ALLELES[y])
        flipped_zscores = zscores * (-1) ** flip_flags
        log.debug("Flipped {} alleles to match reference".format(sum(flip_flags)))
    except KeyError:
        raise ValueError('Incompatible alleles for z-score flipping. Filter invalid SNPs first.')

    return flipped_zscores
