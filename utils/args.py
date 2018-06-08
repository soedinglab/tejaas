#!/usr/bin/env python

'''
    Get user inputs,
    and check for options.
'''

import argparse

class Args():

    def __init__(self):
        args = self.parse_args()
        self.vcf_file  = args.vcf_filename
        self.fam_file  = args.fam_filename
        self.gx_file   = args.gx_filename
        self.method    = args.method
        self.outprefix = args.outprefix

        incsnps = args.incsnps.split(":")
        self.startsnp  = int(incsnps[0].strip())
        self.endsnp    = int(incsnps[1].strip())
        self.psnpcut   = args.psnpthres
        self.pgenecut  = args.pgenethres
        self.maf_file  = args.maf_filename
        self.sigmabeta = args.sigmabeta
        self.jpacut    = args.jpathres
        self.jpafile   = args.jpa_filename
        self.ntransmax = args.optim_ntrans


    def parse_args(self):

        parser = argparse.ArgumentParser(description='Trans-Eqtls from Joint Association AnalysiS (TEJAAS)')
    
        parser.add_argument('--vcf',
                            type=str,
                            dest='vcf_filename',
                            metavar='FILE',
                            help='input VCF file')
    
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
    
        parser.add_argument('--method',
                            default='jpa-rr',
                            type=str,
                            dest='method',
                            metavar='STR',
                            help='which method to run: jpa / rr / optim / jpa-rr')
    
        parser.add_argument('--outprefix',
                            type=str,
                            dest='outprefix',
                            metavar='STR',
                            help='prefix for all output files')
    
        parser.add_argument('--include-SNPs',
                            type=str,
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
                            type=str,
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
                            help='RR statistics are calculated for SNPs above this threshold of JPA statistics')
    
        parser.add_argument('--jpafile',
                            type=str,
                            dest='jpa_filename',
                            help='name of jpa output file; required for selecting SNPs for optimization of sigmabeta')
    
        parser.add_argument('--optim-ntrans',
                            type=int,
                            dest='optim_ntrans',
                            help='number of trans-eQTLs assumed by TEJAAS for optimization of sigma_beta')
    
        res = parser.parse_args()
        return res
