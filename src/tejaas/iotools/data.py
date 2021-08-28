import numpy as np
import scipy.stats as ss
import random
import os
from collections import defaultdict
from sklearn.decomposition import PCA

from tejaas.utils import cismasking
from tejaas.utils import knn

from tejaas.iotools import readgtf
from tejaas.iotools.readOxford import ReadOxford
from tejaas.iotools.readvcf import ReadVCF
from tejaas.iotools.readRPKM import ReadRPKM
from tejaas.utils.containers import GeneInfo, CisMask
from tejaas.utils.logs import MyLogger

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

from functools import wraps
import time

def timeit(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        args[0].logger.debug('{:s} took: {:.6f} seconds'.format(f.__name__, te-ts))
        return result
    return wrap

class Data():

    def __init__(self, args):
        self.logger = MyLogger(__name__)
        self.args = args
        self._gtcent = None
        self._gtnorm = None
        self._snpinfo = None
        self._geneinfo = None
        self._expr = None
        self._cismaskcomp = None
        self._cismasklist = None
        self._tgene_gtnorm = None
        self._tgene_gtcent = None
        self._tgene_expr = None

    @property
    def geno_centered(self):
        return self._gtcent

    @property
    def geno_normed(self):
        return self._gtnorm

    @property
    def snpinfo(self):
        return self._snpinfo

    @property
    def geneinfo(self):
        return self._geneinfo

    @property
    def cismasks_comp(self):
        return self._cismaskcomp

    @property
    def cismasks_list(self):
        return self._cismasklist

    @property
    def expression(self):
        return self._expr

    @property
    def tgene_geno_normed(self):
        return self._tgene_gtnorm

    @property
    def tgene_geno_centered(self):
        return self._tgene_gtcent

    @property
    def tgene_expression(self):
        return self._tgene_expr


    def select_donors(self, vcf_donors, expr_donors):
        ''' Make sure that donors are in the same order for both expression and genotype
        '''
        common_donors = [x for x in vcf_donors if x in expr_donors]
        vcfmask = np.array([vcf_donors.index(x) for x in common_donors])
        exprmask = np.array([expr_donors.index(x) for x in common_donors])
        return vcfmask, exprmask


    def select_genes(self, info, names):
        ''' Select genes which would be analyzed. 
            Make sure the indices are not mixed up
        '''
        allowed = [x.ensembl_id for x in info]
        common  = [x for x in names if x in allowed]
        genes = [x for x in info if x.ensembl_id in common]
        indices = [names.index(x.ensembl_id) for x in genes]
        return genes, np.array(indices)


    def match_gx_indices(self, ref_gx, ref_donors, ref_gnames, gx, donors, gnames):
        '''Match the indices of gx with those of ref_gx
           Both gx and ref_gx are of size G x N
           G = genes (gnames), N = donors
        '''
        gidx = np.array([gnames.index(x) for x in ref_gnames if x in gnames])
        didx = np.array([donors.index(x) for x in ref_donors if x in donors])
        if (gidx.shape[0] != len(ref_gnames)) or (didx.shape[0] != len(ref_donors)):
            self.logger.error("Gene expression files have different donors and / or gene names. Please check. Program cancelled!")
            raise
        return gx[:, didx][gidx, :]

    
    def HWEcheck(self, x):
        gt = x.tolist()
        f = np.array([0] * 3)
        f[0] = gt.count(0)
        f[1] = gt.count(1)
        f[2] = gt.count(2)
        n = sum(f)
        X2 = n * ( (4 * f[0] * f[2] - f[1] ** 2) / ((2 * f[0] + f[1]) * (2 * f[2] + f[1])) )**2
        pval = 1 - ss.chi2.cdf(X2, 1)
        return pval


    def filter_snps(self, snpinfo, dosage, maf_limit = 0.01, use_hwe = False):
        # Predixcan style filtering of snps
        newsnps = list()
        newdosage = list()
        npoly = 0
        nambi = 0
        nunkn = 0
        nlowf = 0
        nlowf_actual = 0
        nhwep = 0
        nalle = 0
        for i, snp in enumerate(snpinfo):
            pos = snp.bp_pos
            refAllele = snp.ref_allele
            effectAllele = snp.alt_allele
            rsid = snp.varid
            maf = round(snp.maf, 3)
            # Actual MAF is lower / higher than population MAF because some samples have been removed
            maf_actual = sum(dosage[i]) / 2 / len(dosage[i])
            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                npoly += 1
                continue
            # Skip unknown alleles
            if refAllele not in SNP_COMPLEMENT or effectAllele not in SNP_COMPLEMENT:
                nalle += 1
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                nambi += 1
                continue
            # Skip unknown RSIDs
            if rsid == '.':
                nunkn += 1
                continue
            # Skip low MAF
            if not (maf >= maf_limit and maf <= (1 - maf_limit)):
                nlowf += 1
                continue
            # Skip low actual MAF
            if not (maf_actual >= maf_limit and maf_actual <= (1 - maf_limit)):
                nlowf_actual += 1
                continue
            # Check HWE
            if use_hwe:
                # Convert to integers 0, 1 or 2
                bins = [0.66, 1.33]
                intdosage = np.digitize(dosage[i], bins)
                # Remove SNPs out of HWE
                hwep = self.HWEcheck(intdosage)
                if(hwep < 0.000001):
                    nhwep += 1
                    # self.logger.debug("SNP {:s} has a HWE p-value of {:g}".format(rsid, hwep))
                    continue
            new_snp = snp._replace(maf = maf_actual)
            newsnps.append(new_snp)
            newdosage.append(dosage[i])
        self.logger.debug("Removed {:d} SNPs because of non-single letter polymorphisms".format(npoly))
        self.logger.debug("Removed {:d} SNPs because of unknown allele symbol".format(nalle))
        self.logger.debug("Removed {:d} SNPs because of ambiguous strands".format(nambi))
        self.logger.debug("Removed {:d} SNPs because of unknown RSIDs".format(nunkn))
        self.logger.debug("Removed {:d} SNPs because of low MAF < {:g}".format(nlowf, maf_limit))
        self.logger.debug("Removed {:d} SNPs because of low MAF (current)".format(nlowf_actual))
        if use_hwe: self.logger.debug("Removed {:d} SNPs because of deviation from HWE".format(nhwep))
        return newsnps, np.array(newdosage)


    def normalize_and_center_dosage(self, dosage):
        f = [snp.maf for snp in self._snpinfo]
        f = np.array(f).reshape(-1, 1)
        gtnorm = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
        gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)
        return gtnorm, gtcent


    def load(self):
        ## Read Oxford File
        if self.args.oxf_file:
            oxf = ReadOxford(self.args.oxf_file, self.args.fam_file, self.args.startsnp, self.args.endsnp, isdosage=self.args.isdosage)
            dosage = oxf.dosage
            gt_donor_ids = oxf.samplenames
            snpinfo = oxf.snpinfo

        # Read VCF file
        if self.args.vcf_file:
            vcf = ReadVCF(self.args.vcf_file, self.args.startsnp, self.args.endsnp, samplefile=self.args.fam_file)
            dosage = vcf.dosage
            gt_donor_ids = vcf.donor_ids
            snpinfo = vcf.snpinfo

        # Read Gene Expression
        self.logger.debug("Reading expression levels for trans-eQTL discovery")
        rpkm = ReadRPKM(self.args.gx_file, self.args.gx_datafmt)
        expression = rpkm.expression
        expr_donors = rpkm.donor_ids
        gene_names = rpkm.gene_names

        # Read confounder corrected gene expression
        if self.args.gxcorr_file is not None:
            self.logger.debug("Reading expression levels for target gene discovery")
            rpkm_corr = ReadRPKM(self.args.gxcorr_file, self.args.gx_datafmt)
            exprcorr = self.match_gx_indices(expression, expr_donors, gene_names, rpkm_corr.expression, rpkm_corr.donor_ids, rpkm_corr.gene_names)

        self.logger.debug("Found {:d} genes of {:d} samples".format(expression.shape[0], expression.shape[1]))
        self.logger.debug("Reading gencode file for gene information")

        gene_info = readgtf.gencode(self.args.gtf_file, trim=self.args.gxtrim, biotype=self.args.biotype)
        if len(gene_info) == 0:
            self.logger.error("No gene annotations found in GTF file. Check feature 'gene' is present, with 'gene_id' and 'gene_name' annotations")

        # reorder donors gt and expr
        self.logger.debug("Selecting common samples of genotype and gene expression")
        self.logger.debug("Before expression selection: {:d} genes from {:d} samples".format(expression.shape[0], expression.shape[1]))
        vcfmask, exprmask = self.select_donors(gt_donor_ids, expr_donors)
        genes, indices = self.select_genes(gene_info, gene_names)       
        expression_selected = rpkm._normalize_expr(expression[:, exprmask][indices, :])
        if self.args.gxcorr_file is not None:
            exprcorr_selected = rpkm_corr._normalize_expr(exprcorr[:, exprmask][indices, :])
        self._geneinfo = genes

        dosage_masked = dosage[:, vcfmask]
        snpinfo_filtered, dosage_filtered_selected = self.filter_snps(snpinfo, dosage_masked)
        self.logger.debug("{:d} SNPs after filtering".format(len(snpinfo_filtered)))
        self._snpinfo = snpinfo_filtered

        self.logger.debug("After expression selection: {:d} genes from {:d} samples".format(expression_selected.shape[0], expression_selected.shape[1]))
        self.logger.debug("Retained {:d} samples".format(vcfmask.shape[0]))

        ### Until here, all filters have been applied and geneinfo and snpinfo reflect current data ###

        self._tgene_gtnorm, self._tgene_gtcent = self.normalize_and_center_dosage(dosage_filtered_selected)
        if self.args.gxcorr_file is not None:
            self._tgene_expr = exprcorr_selected
        else:
            self._tgene_expr = expression_selected

        if self.args.cismasking:
            self.logger.debug("Generate cis-masks for GX matrix for each SNP")
            self._cismasklist = cismasking.get_cismasklist(self._snpinfo, self._geneinfo, self.args.chrom, window=self.args.window)
            self._cismaskcomp = cismasking.compress_cismasklist(self._cismasklist)
            if self.args.crossmapfile is not None:
                self._cismaskcomp = cismasking.extend_cismask(self._geneinfo, self._cismaskcomp, self.args.crossmapfile )

        if self.args.knncorr:
            self.logger.debug("Applying KNN correction on gene expression and genotype")
            gx_corr, gt_corr = knn.knn_correction(expression_selected.T, dosage_filtered_selected, self.args.knn_nbr)
            self._expr = rpkm._normalize_expr(gx_corr.T)
            self._gtnorm, self._gtcent = self.normalize_and_center_dosage(gt_corr)
        else:
            self.logger.debug("No KNN correction.")
            self._expr = expression_selected
            self._gtnorm = self._tgene_gtnorm.copy()
            self._gtcent = self._tgene_gtcent.copy()
            # self._gtnorm, self._gtcent = self.normalize_and_center_dosage(dosage_filtered_selected)

        if self.args.shuffle:
            usedmask = [gt_donor_ids[i] for i in vcfmask]
            if self.args.shuffle_file is not None and os.path.isfile(self.args.shuffle_file):
                self.logger.warn("Shuffling genotype using supplied donor IDs")
                rand_donor_ids = [line.strip() for line in open(self.args.shuffle_file)]
            else:
                self.logger.warn("Shuffling genotype randomly")
                rand_donor_ids = usedmask.copy()
                random.shuffle(rand_donor_ids)
            rand_index = np.array([usedmask.index(x) for x in rand_donor_ids if x in usedmask])
            self._gtnorm = self._gtnorm[:, rand_index]
            self._gtcent = self._gtcent[:, rand_index]
