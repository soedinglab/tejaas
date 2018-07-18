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
        raise argparse.ArgumentTypeError('End SNP position must be less than start SNP position')
    return mlist


def method_strings(mstring):
    try:
        assert (mstring == 'jpa' or mstring == 'jpa-rr' or mstring == 'rr')
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

        self.vcf_file  = args.vcf_filename
        self.oxf_file  = args.oxf_filename
        self.fam_file  = args.fam_filename
        self.gx_file   = args.gx_filename
        self.gtf_file  = args.gtf_filename
        self.method    = args.method
        self.outprefix = args.outprefix
        if args.incsnps is not None:
            self.startsnp = args.incsnps[0]
            self.endsnp   = args.incsnps[1]
        else:
            self.startsnp  = 0
            self.endsnp    = 1
        self.psnpcut   = args.psnpthres
        self.pgenecut  = args.pgenethres
        self.maf_file  = args.maf_filename
        self.sigmabeta = args.sigmabeta
        self.jpacut    = args.jpathres
        self.jpafile   = args.jpa_filename
        self.ntransmax = args.optim_ntrans
        self.nullmodel = args.nullmodel

        self.jpa, self.rr = project.method_selector(args.method)

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
                self.logger.info('NO GENE EXPRESSION FILE PROVIDED. SIMULATING GENE EXPRESSION.')


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
                            help='which method to run: jpa / rr / optim / jpa-rr')
    
        parser.add_argument('--outprefix',
                            type=str,
                            dest='outprefix',
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
    
        parser.add_argument('--maf',
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

        res = parser.parse_args()
        return res
