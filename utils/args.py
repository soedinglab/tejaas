#!/usr/bin/env python

import argparse
import logging
from utils.logs import MyLogger
from utils import project


class Error(Exception):
    pass

class LargerEndSnpError(Error):
    pass

def snprange(mstring):
    try:
        incsnps = mstring.split(":")
        startsnp  = int(incsnps[0].strip())
        endsnp    = int(incsnps[1].strip())
        mlist = [startsnp, endsnp]
        if endsnp < startsnp:
            raise LargerEndSnpError
    except IndexError:
        raise argparse.ArgumentTypeError('Value has to be a range of SNPs separated by colon')
    except LargerEndSnpError:
        raise argparse.ArgumentTypeError('End SNP position must be higher than start SNP position')
    return mlist


def method_strings(mstring):
    try:
        assert (mstring == 'jpa' or mstring == 'jpa-rr' or mstring == 'rr' or mstring == 'rr-sparse')
    except AssertionError:
        raise argparse.ArgumentTypeError('Please specify a correct method')
    return mstring


def null_strings(mstring):
    try:
        assert (mstring == 'perm' or mstring == 'maf')
    except AssertionError:
        raise argparse.ArgumentTypeError('Please specify a correct null model')
    return mstring

class Args():

    def __init__(self, comm, rank):

        self.logger = MyLogger(__name__)
        self.rank = rank
        self.comm = comm

        args = None
        if self.rank == 0:
            args = self.parse_args()
        args = self.comm.bcast(args, root = 0)

        self.vcf_file    = args.vcf_filename
        self.oxf_file    = args.oxf_filename
        # self.oxf_columns = args.oxf_columns
        self.selected_donors = args.selected_donors
        self.isdosage    = args.isdosage
        self.gxtrim      = args.gxtrim
        self.forcetrans  = args.forcetrans
        self.forcecis    = args.forcecis
        self.cismasking  = args.cismasking

        self.shuffle        = args.shuffle
        self.shuffle_file   = args.shuffle_file
        if self.shuffle_file is not None:
            self.shuffle = True

        self.masking(self.forcecis, self.forcetrans, self.cismasking)

        if args.chrom is not None:
            self.chrom   = int(args.chrom)
        else:
            self.chrom   = None
        self.fam_file    = args.fam_filename
        self.gx_file     = args.gx_filename
        self.gtf_file    = args.gtf_filename
        self.outprefix   = args.outprefix
        if args.incsnps is not None:
            self.startsnp = args.incsnps[0] - 1
            self.endsnp   = args.incsnps[1]
        else:
            self.startsnp  = 0
            self.endsnp    = 1e15 # an unusually high number to ensure all SNPs are read.
        self.psnpcut   = args.psnpthres
        self.pgenecut  = args.pgenethres
        self.maf_file  = args.maf_filename
        self.sigmabeta = args.sigmabeta
        self.jpacut    = args.jpathres
        self.jpafile   = args.jpa_filename
        self.ntransmax = args.optim_ntrans
        self.nullmodel = args.nullmodel
        self.window    = args.window
        self.qnullfile = args.qnullfile

        self.jpa, self.rr, self.pms, self.onlyjpa = project.method_selector(args.method)

        self.simulate  = args.simulate
        self.simparams = args.simparams
        self.maketest  = args.maketest

        self.gxsim = False
        if self.gx_file is None:
            self.gxsim = True

        if self.rank == 0:
            self.logger.info('Method: {:s}'.format(args.method))
            if self.rr:
                self.logger.info('Null Model: {:s}'.format(args.nullmodel))
                self.logger.info('Sigma_beta: {:g}'.format(args.sigmabeta))
            if self.gxsim:
                self.logger.warn('No gene expression file provided. Simulating gene expression')

    def masking(self, forcecis, forcetrans, cismasking):
        if sum([forcecis, forcetrans]) > 1:
            raise argparse.ArgumentTypeError("Mutually exclusive options, choose one of --force-trans or --force-cis")

    def parse_args(self):

        self.logger.info('Running TEJAAS v{:s}'.format(project.version()))

        parser = argparse.ArgumentParser(description='Trans-Eqtls from Joint Association AnalysiS (TEJAAS)')
    
        parser.add_argument('--vcf',
                            type=str,
                            dest='vcf_filename',
                            metavar='FILE',
                            help='input VCF file')

        parser.add_argument('--oxf',
                            type=str,
                            dest='oxf_filename',
                            metavar='FILE',
                            help='input Oxford file')

        # parser.add_argument('--oxf-columns',
        #                     type=int,
        #                     default=6,
        #                     dest='oxf_columns',
        #                     metavar='FILE',
        #                     help='number of columns before genotype data (rsid, bpos, chr, refAllele, altAllele, maf)')

        parser.add_argument('--dosage',
                            dest='isdosage',
                            action='store_true',
                            help='Read dosages')

        parser.add_argument('--no-dosage',
                            dest='isdosage',
                            action='store_false',
                            help='Do not read dosages')

        parser.set_defaults(isdosage=False)

        parser.add_argument('--trim',
                            dest='gxtrim',
                            action='store_true',
                            help='Whether to trim version number from gene Ensembl IDs')

        parser.add_argument('--shuffle',
                            dest='shuffle',
                            action='store_true',
                            help='Shuffle the genotypes randomly')

        parser.add_argument('--shuffle-with',
                            type=str,
                            dest='shuffle_file',
                            metavar='FILE',
                            help='Shuffle the genotypes using the supplied donor IDs file')

        parser.add_argument('--selected-donors',
                            type=str,
                            dest='selected_donors',
                            metavar='FILE',
                            help='List of donors ids to select')

        parser.add_argument('--force-trans',
                            dest='forcetrans',
                            action='store_true',
                            help='Select only trans SNPs-Genes pairs not in the same Chromosome')

        # parser.set_defaults(forcetrans=False)

        parser.add_argument('--force-cis',
                            dest='forcecis',
                            action='store_true',
                            help='Select only cis SNPs-Genes pairs in the same Chromosome')

        # parser.set_defaults(forcecis=False)

        parser.add_argument('--cismask',
                            dest='cismasking',
                            action='store_true',
                            help='Generate cismasks for the expression matrix for each SNP')

        # parser.set_defaults(cismasking=False)

        parser.add_argument('--chrom',
                            dest='chrom',
                            metavar='NUMBER',
                            help="Chromosome number being processed")

        parser.add_argument('--fam',
                            type=str,
                            dest='fam_filename',
                            metavar='FILE',
                            help='input fam file')
    
        parser.add_argument('--gx',
                            type=str,
                            dest='gx_filename',
                            metavar='FILE',
                            help='input expression file')

        parser.add_argument('--gtf',
                            type=str,
                            dest='gtf_filename',
                            metavar='FILE',
                            help='input gtf file')
    
        parser.add_argument('--method',
                            default='jpa-rr',
                            type=method_strings,
                            dest='method',
                            metavar='STR',
                            help='which method to run: jpa / rr / optim / jpa-rr / rr-sparse')
    
        parser.add_argument('--outprefix',
                            type=str,
                            dest='outprefix',
                            default='out',
                            metavar='STR',
                            help='prefix for all output files')
    
        parser.add_argument('--include-SNPs',
                            type=snprange,
                            dest='incsnps',
                            metavar='START:END',
                            help='colon-separated index of SNPs to be included')
    
        parser.add_argument('--psnpthres',
                            default=0.05,
                            type=float,
                            dest='psnpthres',
                            metavar='PVAL',
                            help='target genes will be reported for trans-eQTLs, which are below this threshold p-value for RR/JPA statistics')
    
        parser.add_argument('--pgenethres',
                            default=0.05,
                            type=float,
                            dest='pgenethres',
                            metavar='PVAL',
                            help='target genes whose linear regression association with trans-eQTLs are below this threshold p-value will be reported')
    
        parser.add_argument('--null',
                            default='perm',
                            type=null_strings,
                            dest='nullmodel',
                            metavar='STR',
                            help='which null model to use: perm / maf. The later requires separate maf file')
    
        parser.add_argument('--maf-file',
                            type=str,
                            dest='maf_filename',
                            metavar='FILE',
                            help='file name of the MAF, see Documentation for filetype')
    
        parser.add_argument('--prior-sigma',
                            default=0.005,
                            type=float,
                            dest='sigmabeta',
                            metavar='FLOAT',
                            help='standard deviation of the normal prior for reverse multiple linear regression')
    
        parser.add_argument('--jpathres',
                            default=20,
                            type=float,
                            dest='jpathres',
                            metavar='FLOAT',
                            help='RR statistics are calculated for SNPs above this threshold of JPA statistics')
    
        parser.add_argument('--jpafile',
                            type=str,
                            dest='jpa_filename',
                            metavar='FILE',
                            help='name of jpa output file; required for selecting SNPs for optimization of sigmabeta')
    
        parser.add_argument('--optim-ntrans',
                            type=int,
                            dest='optim_ntrans',
                            metavar='INT',
                            help='number of trans-eQTLs assumed by TEJAAS for optimization of sigma_beta')

        parser.add_argument('--simulate',
                            dest='simulate',
                            action='store_true',
                            help='perform simulation, modify simparams to change default values')

        parser.add_argument('--simparams',
                            nargs='*',
                            default=['0.1', '0.9', '1000', '0', '0.005'],
                            dest='simparams',
                            help='fmin, fmax, nsnp')

        parser.add_argument('--test',
                            dest='maketest',
                            action='store_true',
                            help='whether to do test run')

        parser.add_argument('--window',
                            type = int,
                            default = 1e6,
                            dest = 'window',
                            help = 'Window (number of base pairs) used for masking cis genes')

        parser.add_argument('--nullfile',
                            type = str,
                            dest = 'qnullfile',
                            help = 'Filename for storing / reading null JPA scores')

        res = parser.parse_args()
        return res
