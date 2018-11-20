#!/usr/bin/env python

# This is a modified munge_sumstat.py from LDSC project
# FIZI needed a similar tool but with a few extra output columns
# Credit to Brendan Bulik-Sullivan and Hilary Finucane

from __future__ import division
import argparse
import bz2
import gzip
import logging
import os
import re
import sys

import pyfocus
import numpy as np
import pandas as pd

from scipy.stats import chi2

np.seterr(invalid='ignore')

null_values = {
    'LOG_ODDS': 0,
    'BETA': 0,
    'OR': 1,
    'Z': 0
}

default_cnames = {
    # Chromosome
    'CHR': 'CHR',
    'CHROM': 'CHR',

    # BP
    'BP': 'BP',
    'POS': 'BP',

    # RS NUMBER
    'SNP': 'SNP',
    'MARKERNAME': 'SNP',
    'SNPID': 'SNP',
    'RS': 'SNP',
    'RSID': 'SNP',
    'RS_NUMBER': 'SNP',
    'RS_NUMBERS': 'SNP',
    # NUMBER OF STUDIES
    'NSTUDY': 'NSTUDY',
    'N_STUDY': 'NSTUDY',
    'NSTUDIES': 'NSTUDY',
    'N_STUDIES': 'NSTUDY',
    # P-VALUE
    'P': 'P',
    'PVALUE': 'P',
    'P_VALUE':  'P',
    'PVAL': 'P',
    'P_VAL': 'P',
    'GC_PVALUE': 'P',
    # ALLELE 1
    'A1': 'A1',
    'ALLELE1': 'A1',
    'ALLELE_1': 'A1',
    'EFFECT_ALLELE': 'A1',
    'REFERENCE_ALLELE': 'A1',
    'INC_ALLELE': 'A1',
    'EA': 'A1',
    # ALLELE 2
    'A2': 'A2',
    'ALLELE2': 'A2',
    'ALLELE_2': 'A2',
    'OTHER_ALLELE': 'A2',
    'NON_EFFECT_ALLELE': 'A2',
    'DEC_ALLELE': 'A2',
    'NEA': 'A2',
    # N
    'N': 'N',
    'NCASE': 'N_CAS',
    'CASES_N': 'N_CAS',
    'N_CASE': 'N_CAS',
    'N_CASES': 'N_CAS',
    'N_CONTROLS': 'N_CON',
    'N_CAS': 'N_CAS',
    'N_CON': 'N_CON',
    'NCONTROL': 'N_CON',
    'CONTROLS_N': 'N_CON',
    'N_CONTROL': 'N_CON',
    'WEIGHT': 'N',  # metal does this. possibly risky.
    # SIGNED STATISTICS
    'ZSCORE': 'Z',
    'Z-SCORE': 'Z',
    'GC_ZSCORE': 'Z',
    'Z': 'Z',
    'OR': 'OR',
    'B': 'BETA',
    'BETA': 'BETA',
    'LOG_ODDS': 'LOG_ODDS',
    'EFFECTS': 'BETA',
    'EFFECT': 'BETA',
    'SIGNED_SUMSTAT': 'SIGNED_SUMSTAT',
    # INFO
    'INFO': 'INFO',
    # MAF
    'EAF': 'FRQ',
    'FRQ': 'FRQ',
    'MAF': 'FRQ',
    'FRQ_U': 'FRQ',
    'F_U': 'FRQ',
}

describe_cname = {
    'CHR': 'Chromsome',
    'BP': 'Base position',
    'SNP': 'Variant ID (e.g., rs number)',
    'P': 'p-Value',
    'A1': 'Allele 1, interpreted as ref allele for signed sumstat',
    'A2': 'Allele 2, interpreted as non-ref allele for signed sumstat',
    'N': 'Sample size',
    'N_CAS': 'Number of cases',
    'N_CON': 'Number of controls',
    'Z': 'Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'OR': 'Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)',
    'BETA': '[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'LOG_ODDS': 'Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)',
    'INFO': 'INFO score (imputation quality; higher --> better imputation)',
    'FRQ': 'Allele frequency',
    'SIGNED_SUMSTAT': 'Directional summary statistic as specified by --signed-sumstats',
    'NSTUDY': 'Number of studies in which the SNP was genotyped'
}

numeric_cols = ['P', 'N', 'N_CAS', 'N_CON', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ', 'SIGNED_SUMSTAT', 'NSTUDY']


def parse_pos(pos, option):
    """
    Parse a specified genomic position.
    Should be digits followed optionally by case-insensitive Mb or Kb modifiers.
    """
    match = re.match("^(([0-9]*[.])?[0-9]+)(mb|kb)?$", pos, flags=re.IGNORECASE)
    if match:
        pos_tmp = float(match.group(1))  # position
        pos_mod = match.group(3)  # modifier
        if pos_mod:
            pos_mod = pos_mod.upper()
            if pos_mod == "MB":
                pos_tmp *= 1000000
            elif pos_mod == "KB":
                pos_tmp *= 1000

        position = pos_tmp
    else:
        raise ValueError("Option {} {} is an invalid genomic position".format(option, pos))

    return position


def get_command_string(args):
    """
    Format pyfocus call and options into a string for logging/printing
    """
    base = "focus.py " + " ".join(args[:2]) + os.linesep
    rest = args[2:]
    rest_strs = []
    for cmd in rest:
        if "--" in cmd:
            if cmd in ["--quiet", "--verbose"]:
                rest_strs.append("\t{}".format(cmd) + os.linesep)
            else:
                rest_strs.append("\t{}".format(cmd))
        else:
            rest_strs.append(" " + cmd + os.linesep)

    return base + "".join(rest_strs) + os.linesep


def read_header(fh):
    """Read the first line of a file and returns a list with the column names."""
    (openfunc, compression) = get_compression(fh)
    return [x.rstrip('\n') for x in openfunc(fh).readline().split()]


def get_cname_map(flag, default, ignore):
    """
    Figure out which column names to use.

    Priority is
    (1) ignore everything in ignore
    (2) use everything in flags that is not in ignore
    (3) use everything in default that is not in ignore or in flags

    The keys of flag are cleaned. The entries of ignore are not cleaned. The keys of defualt
    are cleaned. But all equality is modulo clean_header().

    """
    clean_ignore = [clean_header(x) for x in ignore]
    cname_map = {x: flag[x] for x in flag if x not in clean_ignore}
    cname_map.update(
        {x: default[x] for x in default if x not in clean_ignore + flag.keys()})
    return cname_map


def get_compression(fh):
    """
    Read filename suffixes and figure out whether it is gzipped,bzip2'ed or not compressed
    """
    if fh.endswith('gz'):
        compression = 'gzip'
        openfunc = gzip.open
    elif fh.endswith('bz2'):
        compression = 'bz2'
        openfunc = bz2.BZ2File
    else:
        openfunc = open
        compression = None

    return openfunc, compression


def clean_header(header):
    """
    For cleaning file headers.
    - convert to uppercase
    - replace dashes '-' with underscores '_'
    - replace dots '.' (as in R) with underscores '_'
    - remove newlines ('\n')
    """
    return header.upper().replace('-', '_').replace('.', '_').replace('\n', '')


def filter_pvals(pvals, log):
    """Remove out-of-bounds P-values"""

    ii = (pvals > 0) & (pvals <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        msg = '{N} SNPs had P outside of (0,1]. The P column may be mislabeled'
        log.warning(msg.format(N=bad_p))

    return ii


def filter_info(info, log, args):
    """Remove INFO < args.info_min (default 0.9) and complain about out-of-bounds INFO."""

    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= args.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = (((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()).any(axis=1))
        ii = (info.sum(axis=1) >= args.info_min * (len(info.columns)))
    else:
        raise ValueError('Expected pd.DataFrame or pd.Series.')

    bad_info = jj.sum()
    if bad_info > 0:
        msg = '{N} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled'
        log.warning(msg.format(N=bad_info))

    return ii


def filter_frq(frq, log, args):
    """
    Filter on MAF. Remove MAF < args.maf_min and out-of-bounds MAF.
    """
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        msg = '{N} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled'
        log.warning(msg.format(N=bad_frq))

    frq = np.minimum(frq, 1 - frq)
    ii = frq > args.maf_min
    return ii & ~jj


def filter_alleles(a):
    """Remove alleles that do not describe strand-unambiguous SNPs"""
    return a.isin(pyfocus.VALID_SNPS)


def parse_dat(dat_gen, convert_colname, merge_alleles, log, args):
    """Parse and filter a sumstats file chunk-wise"""

    tot_snps = 0
    dat_list = []
    msg = 'Reading sumstats from {F} into memory {N} SNPs at a time'
    log.info(msg.format(F=args.sumstats, N=int(args.chunksize)))
    drops = {'NA': 0, 'P': 0, 'INFO': 0,
             'FRQ': 0, 'A': 0, 'SNP': 0, 'MERGE': 0}
    for block_num, dat in enumerate(dat_gen):
        log.info('Reading SNP chunk {}'.format(block_num + 1))
        tot_snps += len(dat)
        old = len(dat)
        
        for c in dat.columns:
            # sometimes column types change when streaming the data
            if c in numeric_cols and not np.issubdtype(dat[c].dtype, np.number):
                log.warning('Column {} expected to be numeric. Attempting to convert'.format(c))
                dat[c] = pd.to_numeric(dat[c], errors="coerce")

        dat = dat.dropna(axis=0, how="any", subset=filter(
            lambda x: x != 'INFO', dat.columns)).reset_index(drop=True)
        drops['NA'] += old - len(dat)
        dat.columns = map(lambda x: convert_colname[x], dat.columns)

        ii = np.ones(len(dat), dtype=bool)
        if args.merge_alleles:
            old = ii.sum()
            ii = dat.SNP.isin(merge_alleles.SNP)
            drops['MERGE'] += old - ii.sum()
            if ii.sum() == 0:
                continue

            dat = dat[ii].reset_index(drop=True)
            ii = np.ones(len(dat), dtype=bool)

        if 'INFO' in dat.columns:
            old = ii.sum()
            ii &= filter_info(dat['INFO'], log, args)
            new = ii.sum()
            drops['INFO'] += old - new

        if 'FRQ' in dat.columns:
            old = ii.sum()
            ii &= filter_frq(dat['FRQ'], log, args)
            new = ii.sum()
            drops['FRQ'] += old - new

        old = ii.sum()
        if args.keep_maf:
            dat.drop(
                [x for x in ['INFO'] if x in dat.columns], inplace=True, axis=1)
        else:
            dat.drop(
                [x for x in ['INFO', 'FRQ'] if x in dat.columns], inplace=True, axis=1)

        ii &= filter_pvals(dat.P, log)
        new = ii.sum()
        drops['P'] += old - new
        old = new
        if not args.no_alleles:
            dat.A1 = dat.A1.str.upper()
            dat.A2 = dat.A2.str.upper()
            ii &= filter_alleles(dat.A1 + dat.A2)
            new = ii.sum()
            drops['A'] += old - new

        if ii.sum() == 0:
            continue

        dat_list.append(dat[ii].reset_index(drop=True))

    log.info('Done reading SNP chunks')
    dat = pd.concat(dat_list, axis=0).reset_index(drop=True)
    log.info('Read {N} SNPs from --sumstats file'.format(N=tot_snps))
    if args.merge_alleles:
        log.info('Removed {N} SNPs not in --merge-alleles'.format(N=drops['MERGE']))

    log.info('Removed {N} SNPs with missing values'.format(N=drops['NA']))
    log.info('Removed {N} SNPs with INFO <= {I}'.format(N=drops['INFO'], I=args.info_min))
    log.info('Removed {N} SNPs with MAF <= {M}'.format(N=drops['FRQ'], M=args.maf_min))
    log.info('Removed {N} SNPs with out-of-bounds p-values'.format(N=drops['P']))
    log.info('Removed {N} variants that were not SNPs or were strand-ambiguous'.format(N=drops['A']))
    log.info('{N} SNPs remain'.format(N=len(dat)))

    return dat


def process_n(dat, args, log):
    """Determine sample size from --N* flags or N* columns. Filter out low N SNPs.s"""

    if all(i in dat.columns for i in ['N_CAS', 'N_CON']):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat['N'] = N * P / P[N == N.max()].mean()
        dat.drop(['N_CAS', 'N_CON'], inplace=True, axis=1)
        # NB no filtering on N done here -- that is done in the next code block

    if 'N' in dat.columns:
        n_min = args.n_min if args.n_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        dat = dat[dat.N >= n_min].reset_index(drop=True)
        new = len(dat)
        log.info('Removed {M} SNPs with N < {MIN} ({N} SNPs remain)'.format(M=old - new, N=new, MIN=n_min))

    elif 'NSTUDY' in dat.columns and 'N' not in dat.columns:
        nstudy_min = args.nstudy_min if args.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        dat = dat[dat.NSTUDY >= nstudy_min].drop(
            ['NSTUDY'], axis=1).reset_index(drop=True)
        new = len(dat)
        log.info('Removed {M} SNPs with NSTUDY < {MIN} ({N} SNPs remain)'.format(M=old - new, N=new, MIN=nstudy_min))

    if 'N' not in dat.columns:
        if args.N:
            dat['N'] = args.N
            log.info('Using N = {N}'.format(N=args.N))
        elif args.N_cas and args.N_con:
            dat['N'] = args.N_cas + args.N_con
            if args.daner is None:
                msg = 'Using N_cas = {N1}; N_con = {N2}'
                log.info(msg.format(N1=args.N_cas, N2=args.N_con))
        else:
            raise ValueError(
                'Cannot determine N. This message indicates a bug. N should have been checked earlier in the program.')

    return dat


def p_to_z(pvals):
    """Convert P-value and N to standardized beta."""
    return np.sqrt(chi2.isf(pvals, 1))


def pass_median_check(m, expected_median, tolerance):
    """Check that median(x) is within tolerance of expected_median."""
    return np.abs(m - expected_median) <= tolerance


def parse_flag_cnames(args):
    """
    Parse flags that specify how to interpret nonstandard column names.

    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    """

    cname_options = [
        [args.nstudy, 'NSTUDY', '--nstudy'],
        [args.snp, 'SNP', '--snp'],
        [args.N_col, 'N', '--N'],
        [args.N_cas_col, 'N_CAS', '--N-cas-col'],
        [args.N_con_col, 'N_CON', '--N-con-col'],
        [args.a1, 'A1', '--a1'],
        [args.a2, 'A2', '--a2'],
        [args.p, 'P', '--P'],
        [args.frq, 'FRQ', '--nstudy'],
        [args.info, 'INFO', '--info']
    ]
    flag_cnames = {clean_header(x[0]): x[1]
                   for x in cname_options if x[0] is not None}
    if args.info_list:
        try:
            flag_cnames.update(
                {clean_header(x): 'INFO' for x in args.info_list.split(',')})
        except ValueError:
            raise ValueError('The argument to --info-list should be a comma-separated list of column names')

    null_value = None
    if args.signed_sumstats:
        try:
            cname, null_value = args.signed_sumstats.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError:
            raise ValueError('The argument to --signed-sumstats should be column header comma number')

    return [flag_cnames, null_value]


def allele_merge(dat, alleles, log):
    """
    WARNING: dat now contains a bunch of NA's~
    Note: dat now has the same SNPs in the same order as --merge alleles.
    """
    dat = pd.merge(dat, alleles, how='left', on='SNP', sort=False).reset_index(drop=True)
    ii = dat.A1.notnull()
    a1234 = dat.A1[ii] + dat.A2[ii] + dat.MA[ii]
    match = a1234.apply(lambda y: y in pyfocus.MATCH_ALLELES)
    jj = pd.Series(np.zeros(len(dat), dtype=bool))
    jj[ii] = match        # This breaks the dtype sometimes
    jj = jj.astype(bool)  # Enforce bool dtype
    old = ii.sum()
    n_mismatch = (~match).sum()
    if n_mismatch < old:
        N = old - n_mismatch
        log.info(
            'Removed {M} SNPs whose alleles did not match --merge-alleles ({N} SNPs remain).'.format(M=n_mismatch, N=N))
    else:
        raise ValueError('All SNPs have alleles that do not match --merge-alleles.')

    dat.loc[~jj, [i for i in dat.columns if i != 'SNP']] = float('nan')
    dat.drop(['MA'], axis=1, inplace=True)
    return dat


def munge(args):
    log = logging.getLogger(pyfocus.LOG)
    try:
        if args.sumstats is None:
            raise ValueError('The --sumstats flag is required')
        if args.no_alleles and args.merge_alleles:
            raise ValueError(
                '--no-alleles and --merge-alleles are not compatible')
        if args.daner and args.daner_n:
            raise ValueError('--daner and --daner-n are not compatible. Use --daner for sample ' +
                             'size from FRQ_A/FRQ_U headers, use --daner-n for values from Nca/Nco columns')

        file_cnames = read_header(args.sumstats)  # note keys not cleaned
        flag_cnames, signed_sumstat_null = parse_flag_cnames(args)
        if args.ignore:
            ignore_cnames = [clean_header(x) for x in args.ignore.split(',')]
        else:
            ignore_cnames = []

        # remove LOG_ODDS, BETA, Z, OR from the default list
        if args.signed_sumstats is not None or args.a1_inc:
            mod_default_cnames = {x: default_cnames[x] for x in default_cnames if default_cnames[x] not in null_values}
        else:
            mod_default_cnames = default_cnames

        cname_map = get_cname_map(flag_cnames, mod_default_cnames, ignore_cnames)
        if args.daner:
            frq_u = filter(lambda x: x.startswith('FRQ_U_'), file_cnames)[0]
            frq_a = filter(lambda x: x.startswith('FRQ_A_'), file_cnames)[0]
            N_cas = float(frq_a[6:])
            N_con = float(frq_u[6:])
            log.info(
                'Inferred that N_cas = {N1}, N_con = {N2} from the FRQ_[A/U] columns.'.format(N1=N_cas, N2=N_con))
            args.N_cas = N_cas
            args.N_con = N_con
            # drop any N, N_cas, N_con or FRQ columns
            for c in ['N', 'N_CAS', 'N_CON', 'FRQ']:  # c isn't used!?!?
                for d in [x for x in cname_map if cname_map[x] == 'c']:
                    del cname_map[d]

            cname_map[frq_u] = 'FRQ'

        if args.daner_n:
            frq_u = filter(lambda x: x.startswith('FRQ_U_'), file_cnames)[0]
            cname_map[frq_u] = 'FRQ'
            try:
                dan_cas = clean_header(file_cnames[file_cnames.index('Nca')])
            except ValueError:
                raise ValueError('Could not find Nca column expected for daner-n format')
        
            try:
                dan_con = clean_header(file_cnames[file_cnames.index('Nco')])
            except ValueError:
                raise ValueError('Could not find Nco column expected for daner-n format')
    
            cname_map[dan_cas] = 'N_CAS'
            cname_map[dan_con] = 'N_CON'
    
            cname_translation = {x: cname_map[clean_header(x)] for x in file_cnames if
                                 clean_header(x) in cname_map}  # note keys not cleaned
            cname_description = {
                x: describe_cname[cname_translation[x]] for x in cname_translation}
            if args.signed_sumstats is None and not args.a1_inc:
                sign_cnames = [
                    x for x in cname_translation if cname_translation[x] in null_values]
                if len(sign_cnames) > 1:
                    raise ValueError(
                        'Too many signed sumstat columns. Specify which to ignore with the --ignore flag.')
                if len(sign_cnames) == 0:
                    raise ValueError(
                        'Could not find a signed summary statistic column.')
    
                sign_cname = sign_cnames[0]
                signed_sumstat_null = null_values[cname_translation[sign_cname]]
                cname_translation[sign_cname] = 'SIGNED_SUMSTAT'
            else:
                sign_cname = 'SIGNED_SUMSTATS'
    
            # check that we have all the columns we need
            if not args.a1_inc:
                req_cols = ['SNP', 'P', 'SIGNED_SUMSTAT']
            else:
                req_cols = ['SNP', 'P']
    
            for c in req_cols:
                if c not in cname_translation.values():
                    raise ValueError('Could not find {C} column.'.format(C=c))

        # check aren't any duplicated column names in mapping
        for field in cname_translation:
            numk = file_cnames.count(field)
            if numk > 1:
                raise ValueError('Found {num} columns named {C}'.format(C=field, num=str(numk)))
    
            # check multiple different column names don't map to same data field
            for head in cname_translation.values():
                numc = cname_translation.values().count(head)
            if numc > 1:
                    raise ValueError('Found {num} different {C} columns'.format(C=head, num=str(numc)))
    
            if (not args.N) and (not (args.N_cas and args.N_con)) and ('N' not in cname_translation.values()) and\
                    (any(x not in cname_translation.values() for x in ['N_CAS', 'N_CON'])):
                raise ValueError('Could not determine N.')
            if ('N' in cname_translation.values() or all(x in cname_translation.values() for x in ['N_CAS', 'N_CON']))\
                    and 'NSTUDY' in cname_translation.values():
                nstudy = [
                    x for x in cname_translation if cname_translation[x] == 'NSTUDY']
                for x in nstudy:
                    del cname_translation[x]
            if not args.no_alleles and not all(x in cname_translation.values() for x in ['A1', 'A2']):
                raise ValueError('Could not find A1/A2 columns.')
    
            log.info('Interpreting column names as follows:')
            for x in cname_description:
                log.info(x + ': ' + cname_description[x])
    
            if args.merge_alleles:
                log.info(
                    'Reading list of SNPs for allele merge from {F}'.format(F=args.merge_alleles))
                (openfunc, compression) = get_compression(args.merge_alleles)
                merge_alleles = pd.read_csv(args.merge_alleles, compression=compression, header=0,
                                            delim_whitespace=True, na_values='.')
                if any(x not in merge_alleles.columns for x in ["SNP", "A1", "A2"]):
                    raise ValueError(
                        '--merge-alleles must have columns SNP, A1, A2.')
    
                log.info(
                    'Read {N} SNPs for allele merge.'.format(N=len(merge_alleles)))
                merge_alleles['MA'] = (
                    merge_alleles.A1 + merge_alleles.A2).apply(lambda y: y.upper())
                merge_alleles.drop(
                    [x for x in merge_alleles.columns if x not in ['SNP', 'MA']], axis=1, inplace=True)
            else:
                merge_alleles = None
    
            (openfunc, compression) = get_compression(args.sumstats)
    
            # figure out which columns are going to involve sign information, so we can ensure
            # they're read as floats
            signed_sumstat_cols = [k for k, v in cname_translation.items() if v == 'SIGNED_SUMSTAT']
            dat_gen = pd.read_csv(args.sumstats, delim_whitespace=True, header=0,
                                  compression=compression, usecols=cname_translation.keys(),
                                  na_values=['.', 'NA'], iterator=True, chunksize=args.chunksize,
                                  dtype={c: np.float64 for c in signed_sumstat_cols})
    
            dat = parse_dat(dat_gen, cname_translation, merge_alleles, log, args)
            if len(dat) == 0:
                raise ValueError('After applying filters, no SNPs remain.')
    
            old = len(dat)
            dat = dat.drop_duplicates(subset='SNP').reset_index(drop=True)
            new = len(dat)
            log.info('Removed {M} SNPs with duplicated rs numbers ({N} SNPs remain).'.format(M=old - new, N=new))
            # filtering on N cannot be done chunkwise
            dat = process_n(dat, args, log)
            dat.P = p_to_z(dat.P)
            dat.rename(columns={'P': 'Z'}, inplace=True)
            if not args.a1_inc:
                m = np.median(dat.SIGNED_SUMSTAT)
                if not pass_median_check(m, signed_sumstat_null, 0.1):
                    msg = 'Median value of {F} is {V} (should be close to {M}). This column may be mislabeled.'
                    raise ValueError(msg.format(F=sign_cname, M=signed_sumstat_null, V=round(m, 2)))
                else:
                    msg = 'Median value of {F} was {C}, which seems sensible.'.format(C=m, F=sign_cname)
                    log.info(msg)
    
                dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
                dat.drop('SIGNED_SUMSTAT', inplace=True, axis=1)
    
            # do this last so we don't have to worry about NA values in the rest of
            # the program
            if args.merge_alleles:
                dat = allele_merge(dat, merge_alleles, log)
    
            out_fname = args.output + '.sumstats.gz'
            print_colnames = [
                c for c in dat.columns if c in ['CHR', 'BP', 'SNP', 'N', 'Z', 'A1', 'A2']]
            if args.keep_maf and 'FRQ' in dat.columns:
                print_colnames.append('FRQ')
            msg = 'Writing summary statistics for {M} SNPs ({N} with nonmissing beta) to {F}.'
            log.info(
                msg.format(M=len(dat), F=out_fname, N=dat.N.notnull().sum()))
    
            dat.to_csv(
                out_fname, sep="\t", index=False, columns=print_colnames, float_format='%.3f', compression="gzip")
    
            CHISQ = (dat.Z ** 2)
            mean_chisq = CHISQ.mean()
    
            log.info('METADATA - Mean chi^2 = ' + str(round(mean_chisq, 3)))
            if mean_chisq < 1.02:
                log.warning("METADATA: Mean chi^2 may be too small")
    
            log.info('METADATA - Lambda GC = ' + str(round(CHISQ.median() / 0.4549, 3)))
            log.info('METADATA - Max chi^2 = ' + str(round(CHISQ.max(), 3)))
            log.info(
                'METADATA - {N} Genome-wide significant SNPs (some may have been removed by filtering)'.format(N=(CHISQ > 29).sum()))

    except Exception as err:
        log.error(err.message)
    finally:
        log.info('Conversion finished')

    return


def run_twas(args):
    log = logging.getLogger(pyfocus.LOG)
    try:
        # perform sanity arguments checking before continuing
        chrom = None
        start_bp = None
        stop_bp = None
        if any(x is not None for x in [args.chr, args.start, args.stop]):
            if args.start is not None and args.chr is None:
                raise ValueError("Option --start cannot be set unless --chr is specified")
            if args.stop is not None and args.chr is None:
                raise ValueError("Option --stop cannot be set unless --chr is specified")

            chrom = args.chr

            # parse start/stop positions and make sure they're ordered (if exist)
            if args.start is not None:
                start_bp = parse_pos(args.start, "--start")
            else:
                start_bp = None

            if args.stop is not None:
                stop_bp = parse_pos(args.stop, "--stop")
            else:
                stop_bp = None

            if args.start is not None and args.stop is not None:
                if start_bp >= stop_bp:
                    raise ValueError("Specified --start position must be before --stop position")

        # load GWAS summary data
        log.info("Preparing GWAS summary file")
        gwas = pyfocus.GWAS.parse_gwas(args.gwas)

        # load reference genotype data
        log.info("Preparing reference SNP data")
        ref = pyfocus.RefPanel.parse_plink(args.ref)

        log.info("Preparing eQTL weights")
        wcollection = pyfocus.WeightsCollection(args.weights)

        with open("{}.focus.txt".format(args.output), "w") as output:
            if args.locations is not None:
                log.info("Preparing user-defined locations")
                partitions = pyfocus.IndBlocks(args.locations)
            else:
                partitions = pyfocus.IndBlocks()

            partitions = partitions.subset_by_pos(chrom, start_bp, stop_bp)

            written = False
            for region in partitions:
                chrom, start, stop = region
                local_gwas = gwas.subset_by_pos(chrom, start, stop)

                # only fine-map regions that contain GWAS signal
                if len(local_gwas) == 0:
                    log.warning("No SNPs found at {}:{} - {}:{}. Skipping".format(chrom, int(start), chrom, int(stop)))
                    continue

                if not local_gwas.contains_hit(args.p_threshold):
                    log.warning("No SNPs passing threshold {} found at {}:{} - {}:{}. Skipping".format(args.p_threshold,
                                chrom, int(start), chrom, int(stop)))
                    continue

                # grab local weight data
                local_weights = wcollection.subset_by_pos(chrom, start, stop)

                # grab local reference genotype data
                local_ref = ref.subset_by_pos(chrom, start, stop)

                # perform fine-mapping
                result = pyfocus.fine_map(local_gwas, local_weights, local_ref)

                # fine-map can break and return early if there are issues so check for non-none result
                if result is not None:
                    pyfocus.write_output(result, output, append=written)
                    written = True

    except Exception as err:
        log.error(err.message)
    finally:
        log.info("Finished twas & fine-mapping")

    return 0


def build_weights(args):
    pass


def build_munge_parser(subp):
    munp = subp.add_parser("munge",
                           description="Munge summary statistics input to conform to FIZI requirements")

    munp.add_argument('sumstats', type=argparse.FileType("r"),
                      help="Input filename.")
    munp.add_argument('--N', default=None, type=float,
                      help="Sample size If this option is not set, will try to infer the sample "
                           "size from the input file. If the input file contains a sample size "
                           "column, and this flag is set, the argument to this flag has priority.")
    munp.add_argument('--N-cas', default=None, type=float,
                      help="Number of cases. If this option is not set, will try to infer the number "
                           "of cases from the input file. If the input file contains a number of cases "
                           "column, and this flag is set, the argument to this flag has priority.")
    munp.add_argument('--N-con', default=None, type=float,
                      help="Number of controls. If this option is not set, will try to infer the number "
                           "of controls from the input file. If the input file contains a number of controls "
                           "column, and this flag is set, the argument to this flag has priority.")
    munp.add_argument('--info-min', default=0.9, type=float,
                      help="Minimum INFO score.")
    munp.add_argument('--maf-min', default=0.01, type=float,
                      help="Minimum MAF.")
    munp.add_argument('--daner', default=False, action='store_true',
                      help="Use this flag to parse Stephan Ripke's daner* file format.")
    munp.add_argument('--daner-n', default=False, action='store_true',
                      help="Use this flag to parse more recent daner* formatted files, which "
                           "include sample size column 'Nca' and 'Nco'.")
    munp.add_argument('--no-alleles', default=False, action="store_true",
                      help="Don't require alleles. Useful if only unsigned summary statistics are available "
                           "and the goal is h2 / partitioned h2 estimation rather than rg estimation.")
    munp.add_argument('--merge-alleles', default=None, type=str,
                      help="Same as --merge, except the file should have three columns: SNP, A1, A2, "
                           "and all alleles will be matched to the --merge-alleles file alleles.")
    munp.add_argument('--n-min', default=None, type=float,
                      help='Minimum N (sample size). Default is (90th percentile N) / 2.')
    munp.add_argument('--chunksize', default=5e6, type=int,
                      help='Chunksize.')

    # optional args to specify column names
    munp.add_argument('--snp', default=None, type=str,
                      help='Name of SNP column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--N-col', default=None, type=str,
                      help='Name of N column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--N-cas-col', default=None, type=str,
                      help='Name of N column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--N-con-col', default=None, type=str,
                      help='Name of N column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--a1', default=None, type=str,
                      help='Name of A1 column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--a2', default=None, type=str,
                      help='Name of A2 column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--p', default=None, type=str,
                      help='Name of p-value column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--frq', default=None, type=str,
                      help='Name of FRQ or MAF column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--signed-sumstats', default=None, type=str,
                      help='Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.')
    munp.add_argument('--info', default=None, type=str,
                      help='Name of INFO column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--info-list', default=None, type=str,
                      help='Comma-separated list of INFO columns. Will filter on the mean. NB: case insensitive.')
    munp.add_argument('--nstudy', default=None, type=str,
                      help='Name of NSTUDY column (if not a name that pyfocus understands). NB: case insensitive.')
    munp.add_argument('--nstudy-min', default=None, type=float,
                      help='Minimum # of studies. Default is to remove everything below the max, unless there is an N column,'
                           ' in which case do nothing.')
    munp.add_argument('--ignore', default=None, type=str,
                      help='Comma-separated list of column names to ignore.')
    munp.add_argument('--a1-inc', default=False, action='store_true',
                      help='A1 is the increasing allele.')
    munp.add_argument('--keep-maf', default=False, action='store_true',
                      help='Keep the MAF column (if one exists).')

    # misc options
    munp.add_argument("-q", "--quiet", default=False, action="store_true",
                      help="Do not print anything to stdout.")
    munp.add_argument("--verbose", default=False, action="store_true",
                      help="Verbose logging. Includes debug info.")
    munp.add_argument("-o", "--output", default="FIZI",
                      help="Prefix for output data.")

    return munp


def build_focus_parser(subp):
    # add imputation parser
    fmp = subp.add_parser("finemap", description="Perform a TWAS and fine-map regional results.")

    # main arguments
    fmp.add_argument("gwas", type=argparse.FileType("r"),
                      help="GWAS summary data. Supports gzip and bz2 compression.")
    fmp.add_argument("ref",
                      help="Path to reference panel PLINK data.")
    fmp.add_argument("weights",
                     help="Path to eQTL weights.")

    # location options
    fmp.add_argument("--chr", default=None,
                      help="Perform imputation for specific chromosome.")
    fmp.add_argument("--start", default=None,
                      help="Perform imputation starting at specific location (in base pairs). Accepts kb/mb modifiers. Requires --chr to be specified.")
    fmp.add_argument("--stop", default=None,
                      help="Perform imputation until at specific location (in base pairs). Accepts kb/mb modifiers. Requires --chr to be specified.")
    fmp.add_argument("--locations", default=None, type=argparse.FileType("r"),
                     help="Path to a BED file containing windows (e.g., CHR START STOP) to impute.")

    fmp.add_argument("--min-prop", default=0.9, type=float,
                      help="Minimum required proportion of gwas/reference panel overlap to perform twas.")
    fmp.add_argument("--p-threshold", default=5e-8, type=float,
                     help="GWAS p-value threshold required for fine-mapping.")

    return fmp


def build_weights_parser(subp):
    wgtp = subp.add_parser("build", description="Compute weights for downstream TWAS and fine-mapping.")
    return wgtp


def main(argsv):
    # setup main parser
    argp = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subp = argp.add_subparsers(help="Subcommands: munge to clean up summary statistics. finemap to perform run twas & finemap.")

    # add subparser for the  munge summary statistics program
    munp = build_munge_parser(subp)
    munp.set_defaults(func=munge)

    fmp = build_focus_parser(subp)
    fmp.set_defaults(func=run_twas)

    wgtp = build_weights_parser(subp)
    wgtp.set_defaults(func=build_weights)

    # parse arguments
    args = argp.parse_args(argsv)

    cmd_str = get_command_string(argsv)

    masthead =  "===================================" + os.linesep
    masthead += "              FOCUS v{}            ".format(pyfocus.VERSION) + os.linesep
    masthead += "===================================" + os.linesep

    # setup logging
    log_format = "[%(asctime)s - %(levelname)s] %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"
    log = logging.getLogger(pyfocus.LOG)
    log.setLevel(logging.INFO)
    fmt = logging.Formatter(fmt=log_format, datefmt=date_format)

    # write to stdout unless quiet is set
    if not args.quiet:
        sys.stdout.write(masthead)
        sys.stdout.write(cmd_str)
        sys.stdout.write("Starting log..." + os.linesep)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setFormatter(fmt)
        log.addHandler(stdout_handler)

    # setup log file, but write PLINK-style command first
    disk_log_stream = open("{}.log".format(args.output), "w")
    disk_log_stream.write(masthead)
    disk_log_stream.write(cmd_str)
    disk_log_stream.write("Starting log..." + os.linesep)

    disk_handler = logging.StreamHandler(disk_log_stream)
    disk_handler.setFormatter(fmt)
    log.addHandler(disk_handler)

    # launch either munge or impute
    args.func(args)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
