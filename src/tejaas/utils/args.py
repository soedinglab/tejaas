#!/usr/bin/env python
'''
Collect and streamline all arguments and format them.
Sanity check on input values.
'''

import argparse
import os
import errno
from tejaas.utils.logs import MyLogger
from tejaas.utils import project


class InputError(Exception):
    ''' Raise when incorrect options are used in the argument '''
    pass


class LargerEndSnpError(InputError):
    ''' Raise when end SNP index is larger than start SNP index '''
    pass


def snprange(mstring):
    ''' 
    snprange: a colon-separated user input style for specifying start:end of SNPs to read from VCF file.
    This function is a parser for the colon-separated string.
    '''
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


def biotype_fmt(stringlist):
    allowed_types = ['protein_coding', 'lncRNA']
    try:
        assert (all([x in allowed_types for x in stringlist]))
    except AssertionError:
        raise argparse.ArgumentTypeError('Please specify a correct biotype')
    return stringlist


def method_strings(mstring):
    '''
    Check if the specified method name is valid.
    '''
    try:
        assert (mstring == 'jpa' or mstring == 'rr' or mstring == 'jpa-rr')
    except AssertionError:
        raise argparse.ArgumentTypeError('Please specify a correct method')
    return mstring


def null_strings(mstring):
    '''
    Check if the specified null model is valid.
    '''
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
        self.isdosage    = args.isdosage
        self.fam_file    = args.fam_filename
        if args.chrom is not None:
            self.chrom   = int(args.chrom)
        else:
            self.chrom   = None
        self.gx_file     = args.gx_filename
        self.gxcorr_file = args.gxcorr_filename
        self.gx_datafmt  = args.gx_datafmt
        self.gtf_file    = args.gtf_filename
        self.gxtrim      = args.gxtrim
        self.biotype     = args.biotype
        self.outprefix   = args.outprefix
        if args.incsnps is not None:
            self.startsnp = args.incsnps[0] - 1
            self.endsnp   = args.incsnps[1]
        else:
            self.startsnp  = 0
            self.endsnp    = 1e15 # an unusually high number to ensure all SNPs are read.

        self.jpa, self.rr = project.method_selector(args.method)
        self.nullmodel   = args.nullmodel
        self.cismasking  = args.cismasking
        self.window      = args.window
        self.sigmabeta   = args.sigmabeta
        self.knn_nbr     = args.knn
        self.knncorr     = True
        if args.knn == 0:
            self.knncorr = False

        self.shuffle        = args.shuffle
        self.shuffle_file   = args.shuffle_file
        if self.shuffle_file is not None:
            self.shuffle = True

        self.psnpcut   = args.psnpthres
        self.pgenecut  = args.pgenethres
        self.maf_file  = args.maf_filename
        self.jpanull_file = args.qnullfile
        self.jpa_calc_null = project.need_new_jpanull_file(self.jpa, self.jpanull_file)
        self.jpanull_iter = args.qnull_iter
        self.seed      = args.seed

        self.maketest  = args.maketest

        if not self.maketest: self.check_inputs()
        self.crossmapfile = args.crossmapfile
        self.usefdr = False
        self.target_fdr = args.target_fdr
        if self.target_fdr is not None:
            self.usefdr = True


        if self.rank == 0:
            self.logger.info('Method: {:s}'.format(args.method))
            if self.rr:
                self.logger.info('Null Model: {:s}'.format(args.nullmodel))
                self.logger.info('Sigma_beta: {:g}'.format(args.sigmabeta))


    def parse_args(self):

        self.logger.info('Running TEJAAS v{:s}'.format(project.version()))

        parser = argparse.ArgumentParser(description='Tejaas: Discover trans-eQTLs!')
    
        parser.add_argument('--vcf',
                            type=str,
                            dest='vcf_filename',
                            metavar='FILE',
                            help='Input VCF file in vcf.gz format')

        parser.add_argument('--oxf',
                            type=str,
                            dest='oxf_filename',
                            metavar='FILE',
                            help='Input Oxford file')

        parser.add_argument('--dosage',
                            dest='isdosage',
                            action='store_true',
                            help='Read dosages')

        parser.add_argument('--fam',
                            type=str,
                            dest='fam_filename',
                            metavar='FILE',
                            help='Input fam file')
    
        parser.add_argument('--chrom',
                            dest='chrom',
                            metavar='NUMBER',
                            help="Chromosome number of the genotype file")

        parser.add_argument('--include-SNPs',
                            type=snprange,
                            dest='incsnps',
                            metavar='START:END',
                            help='Colon-separated index of SNPs to be included')
    
        parser.add_argument('--gx',
                            type=str,
                            dest='gx_filename',
                            metavar='FILE',
                            help='input expression file for finding trans-eQTLs')

        parser.add_argument('--gxcorr',
                            type=str,
                            dest='gxcorr_filename',
                            metavar='FILE',
                            help='input expression file for finding target genes')

        parser.add_argument('--gxfmt',
                            type=str,
                            dest='gx_datafmt',
                            metavar='GX_FORMAT',
                            default='gtex',
                            help='Format of input gene expression file. Supported: gtex, cardiogencis and geuvadis')

        parser.add_argument('--biotype',
                            nargs='*',
                            type=biotype_fmt,
                            dest='biotype',
                            metavar='BIOTYPE_OPTIONS',
                            default=['protein_coding', 'lncRNA'],
                            help='List of biotypes to be selected from the GENCODE annotation file. Supported options: protein_coding, lncRNA')

        parser.add_argument('--gtf',
                            type=str,
                            dest='gtf_filename',
                            metavar='FILE',
                            help='input gtf file')
    
        parser.add_argument('--trim',
                            dest='gxtrim',
                            action='store_true',
                            help='Trim version number from GENCODE Ensembl IDs')

        parser.add_argument('--outprefix',
                            type=str,
                            dest='outprefix',
                            default='out',
                            metavar='STR',
                            help='prefix for all output files')

        parser.add_argument('--method',
                            default='rr',
                            type=method_strings,
                            dest='method',
                            metavar='STR',
                            help='which method to run: jpa / rr')
    
        parser.add_argument('--null',
                            default='perm',
                            type=null_strings,
                            dest='nullmodel',
                            metavar='STR',
                            help='which null model to use: perm / maf')
    
        parser.add_argument('--cismask',
                            dest='cismasking',
                            action='store_true',
                            help='Generate cismasks for the expression matrix for each SNP')

        parser.add_argument('--window',
                            type = int,
                            default = 1e6,
                            dest = 'window',
                            help = 'Window (number of base pairs) used for masking cis genes')

        parser.add_argument('--prior-sigma',
                            default=0.1,
                            type=float,
                            dest='sigmabeta',
                            metavar='FLOAT',
                            help='standard deviation of the normal prior for reverse multiple linear regression')

        parser.add_argument('--knn',
                            type = int,
                            dest = 'knn',
                            help = 'Number of neighbours for KNN (use 0 if you do not want KNN correction)',
                            default = 0)

        parser.add_argument('--psnpthres',
                            default=0.0001,
                            type=float,
                            dest='psnpthres',
                            metavar='PVAL',
                            help='target genes will be reported for trans-eQTLs, which are below this threshold p-value for RR/JPA statistics')
    
        parser.add_argument('--pgenethres',
                            default=0.001,
                            type=float,
                            dest='pgenethres',
                            metavar='PVAL',
                            help='target genes whose linear regression association with trans-eQTLs are below this threshold p-value will be reported')
    
        parser.add_argument('--jpanull',
                            type = str,
                            dest = 'qnullfile',
                            help = 'Filename for storing / reading null JPA scores')

        parser.add_argument('--jpanull-iter',
                            default = 100000,
                            type = int,
                            dest = 'qnull_iter',
                            help = 'Number of iterations for creating null JPA scores')

        parser.add_argument('--seed',
                            default = None,
                            type = int,
                            dest = 'seed',
                            help = 'Seed the random generator for numpy, used for development purpose')

        parser.add_argument('--maf-file',
                            type=str,
                            dest='maf_filename',
                            metavar='FILE',
                            help='file name of the MAF, see Documentation for filetype')
    
        parser.add_argument('--shuffle',
                            dest='shuffle',
                            action='store_true',
                            help='Shuffle the genotypes randomly')

        parser.add_argument('--shuffle-with',
                            type=str,
                            dest='shuffle_file',
                            metavar='FILE',
                            help='Shuffle the genotypes using the supplied donor IDs file')

        parser.add_argument('--test',
                            dest='maketest',
                            action='store_true',
                            help='whether to do test run')

        parser.add_argument('--crossmap',
                            type = str,
                            default=None,
                            dest = 'crossmapfile',
                            help = 'Crossmapability file (Saha, Battle 2018) ')

        parser.add_argument('--fdrgenethres',
                            default=None,
                            type=float,
                            dest='target_fdr',
                            metavar='FDR',
                            help='enable FDR correction up to a certain cutoff for target gene discovery')

        res = parser.parse_args()
        return res


    def check_inputs(self):
        '''
        Perform sanity checks on the input options.
        '''
        if self.rank == 0:
            '''
            Check if any genotype file is specified.
            '''
            try:
                assert (self.vcf_file is not None) or (self.oxf_file is not None)
            except AssertionError:
                print ('Input error: Specify either --vcf or --oxf. See --help for details.')
                raise

            if (self.oxf_file is not None):
                try:
                    assert (self.fam_file is not None)
                except AssertionError:
                    print ('Input error: Specify the sample file with --fam. See --help for details.')
                    raise

            '''
            Check if gene expression and GTF files are specified
            '''
            try:
                assert (self.gx_file is not None)
            except AssertionError:
                print ('Input error: Specify gene expression file. See --help for details')
                raise
            try:
                assert (self.gtf_file is not None)
            except AssertionError:
                print ('Input error: Specify GENCODE file. See --help for details')
                raise

            '''
            Check if files exist.
            '''
            for filepath in [self.vcf_file, self.oxf_file, self.fam_file, self.gx_file, self.gtf_file, self.gxcorr_file]:
                if (filepath is not None):
                    try:
                        assert os.path.isfile(filepath)
                    except AssertionError:
                        print ('File {:s} does not exist.'.format(filepath))
                        raise

            '''
            Check if output directory is writable / can be created
            '''
            outdir = os.path.dirname(os.path.realpath(self.outprefix))
            try:
                if not os.path.exists(outdir): os.makedirs(outdir)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print ('Unable to create output directory: {:s}'.format(outdir))
                    raise

            try:
                assert os.path.isdir(outdir) and os.access(outdir, os.W_OK | os.X_OK)
                #filepath = "{:s}.write_tester".format(self.outprefix)
                #filehandle = open( filepath, 'w' )
                #filehandle.close()
                #os.remove(filepath)
            except AssertionError:
                print('Unable to create files in {:s}'.format(outdir))
                raise
