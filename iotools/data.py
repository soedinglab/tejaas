import numpy as np
import scipy.stats as ss

from iotools import simulate
from iotools import readgtf
from iotools.readOxford import ReadOxford
from iotools.readRPKM import ReadRPKM
from utils.containers import GeneInfo
import scipy.stats as ss
from collections import defaultdict
from utils.logs import MyLogger

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

from functools import wraps
import time

def timeit(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        print('{:s} took: {:.6f} seconds'.format(f.__name__, te-ts))
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
    def snp_cismasks(self):
        return self._snps_masks_lists

    @property
    def cismasks(self):
        return self._compressed_masks

    # @property
    # def geneinfo_dict(self):
    #     return self._geneinfo_dict

    @property
    def expression(self):
        return self._expr

    @timeit
    def get_cismasks(self, snpinfo, gene_info, chrom):
        # gene_info MUST be sorted as in the gene expression matrix
        # chrom   = snpinfo[0].chrom
        window  = 1e6

        chr_genes = list()
        chr_ix_genes = list()
        for i,g in enumerate(gene_info):
            if g.chrom == chrom:
                chr_genes.append(g)
                chr_ix_genes.append(i)
        chr_ix_genes = np.array(chr_ix_genes)
        
        cis_masks = list()
        prev_start = 0

        for s in snpinfo:
            #here starts for each snp
            pos     = s.bp_pos
            w_start = pos - window
            w_end   = pos + window
            cis_genes_ix = list()
            for i, g in enumerate(chr_genes[prev_start:]):
                if g.start >= w_start and g.start <= w_end:
                    # cis_genes_ix.append([prev_start + i, g])
                    cis_genes_ix.append(prev_start + i)
                elif g.end >= w_start and g.end <= w_end:
                    # cis_genes_ix.append([prev_start + i, g])
                    cis_genes_ix.append(prev_start + i)
                if g.start > w_end:
                    break
            cis_genes_ix = np.array(cis_genes_ix)
            if len(cis_genes_ix) > 0:
                gx_matrix_ix = chr_ix_genes[cis_genes_ix]
                if len(gx_matrix_ix) != len(cis_genes_ix):
                    raise
                if prev_start != cis_genes_ix[0]:
                    prev_start = cis_genes_ix[0]
            else:
                # print(pos, "empty cis")
                gx_matrix_ix = np.array([])
            cis_masks.append(gx_matrix_ix)
        return cis_masks

    @timeit
    def compress_cis_masks(self, cis_masks):
        snps_masks_lists = [] # contains lists of snp indices for each mask.
                              # i.e snps_masks_lists[0] is a list with all the snps 
                              # that have the same mask: compressed_masks[0]
        compressed_snps  = [] # temp list that will be used to build snps_masks_lists     
        compressed_masks = [] # contains all unique masks

        prev_mask = cis_masks[0]
        compressed_snps.append(0)
        for i in range(1, len(cis_masks)):
            # if masks are the same, add the snp index
            if (len(cis_masks[i]) == len(prev_mask)) and all(cis_masks[i] == prev_mask):
                compressed_snps.append(i)
            else:
                compressed_masks.append(cis_masks[i])
                prev_mask = cis_masks[i]
                snps_masks_lists.append(compressed_snps)
                compressed_snps = []
                compressed_snps.append(i)
        # add the last set of snps for the final mask
        snps_masks_lists.append(compressed_snps)
        compressed_masks.append(cis_masks[-1])
        return snps_masks_lists, compressed_masks


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

    
    def HWEcheck(self, x):
        gt = x.tolist()
        f = np.array([0] * 3)
        f[0] = gt.count(0)
        f[1] = gt.count(1)
        f[2] = gt.count(2)
        n = sum(f)
        #p_A = (2 * f[0] + f[1]) / (2 * n)
        #p_a = (2 * f[2] + f[1]) / (2 * n)
        X2 = n * ( (4 * f[0] * f[2] - f[1] ** 2) / ((2 * f[0] + f[1]) * (2 * f[2] + f[1])) )**2
        pval = 1 - ss.chi2.cdf(X2, 1)
        return pval


    def filter_snps(self, snpinfo, dosage):
        # Predixcan style filtering of snps
        newsnps = list()
        newdosage = list()
        npoly = 0
        nambi = 0
        nunkn = 0
        nlowf = 0
        nhwep = 0
        for i, snp in enumerate(snpinfo):
            pos = snp.bp_pos
            refAllele = snp.ref_allele
            effectAllele = snp.alt_allele
            rsid = snp.varid
            maf = snp.maf
            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                npoly += 1
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
            if not (maf >= 0.10 and maf <=0.90):
                nlowf += 1
                continue
            # Convert to integers 0, 1 or 2
            bins = [0.66, 1.33]
            intdosage = np.digitize(dosage[i], bins)
            # Remove SNPs out of HWE
            if(self.HWEcheck(intdosage) < 0.000001):
                nhwep += 1
                continue
            newsnps.append(snp)
            newdosage.append(intdosage)
        self.logger.debug("Removed {:d} SNPs because of non-single letter polymorphisms".format(npoly))
        self.logger.debug("Removed {:d} SNPs because of ambiguous strands".format(nambi))
        self.logger.debug("Removed {:d} SNPs because of unknown RSIDs".format(nunkn))
        self.logger.debug("Removed {:d} SNPs because of low MAF < 0.10".format(nlowf))
        self.logger.debug("Removed {:d} SNPs because of deviation from HWE".format(nhwep))
        return newsnps, np.array(newdosage)


    def normalize_and_center_dosage(self, dosage):
        f = [snp.maf for snp in self._snpinfo]
        f = np.array(f).reshape(-1, 1)
        self._gtnorm = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
        self._gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)


    def load(self):
        # Read Oxford File
        if self.args.oxf_file:
            oxf = ReadOxford(self.args.oxf_file, self.args.fam_file, self.args.startsnp, self.args.endsnp, isdosage=self.args.isdosage)
            dosage = oxf.dosage
            gt_donor_ids = oxf.samplenames
            snpinfo = oxf.snpinfo

        # Read VCF file
        if self.args.vcf_file:
            vcf = ReadVCF(self.args.vcf_file, self.args.startsnp, self.args.endsnp)
            dosage = vcf.dosage
            gt_donor_ids = vcf.donor_ids
            snpinfo = oxf.snpinfo

        snpinfo_filtered, dosage_filtered = self.filter_snps(snpinfo, dosage)
        self.logger.debug("{:d} SNPs after filtering".format(len(snpinfo_filtered)))
        self._snpinfo = snpinfo_filtered

        # Gene Expression
        rpkm = ReadRPKM(self.args.gx_file, "gtex")
        expression = rpkm.expression
        expr_donors = rpkm.donor_ids
        gene_names = rpkm.gene_names

        self.logger.debug("Completed reading expression levels of {:d} genes of {:d} samples".format(expression.shape[0], expression.shape[1]))

        ### for GTEx ###
        if(self.args.isdosage):
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=False)
        ### for Cardiogenics ###
        else:
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=True)

        # reorder donors gt and expr
        self.logger.debug("Selecting common samples of genotype and gene expression")
        vcfmask, exprmask = self.select_donors(gt_donor_ids, expr_donors)
        genes, indices = self.select_genes(gene_info, gene_names)

        self._expr = expression[:, exprmask][indices, :]
        self._geneinfo = genes
        dosage_filtered_selected = dosage_filtered[:, vcfmask]

        ### Until here, all filters have been applied and geneinfo and snpinfo reflect current data ###

        # This holds!! geneinfo is ordered by chrom and start position
        # just for testing
        # for i, g in enumerate(genes):
        #     if g.start > genes[i+1].start and g.chrom == genes[i+1].chrom:
        #         print(g, genes[i+1])

        if self.args.cismasking:
            self.logger.debug("Generate cis-masks for GX matrix for each SNP")
            cis_masks = self.get_cismasks(self._snpinfo, self._geneinfo, self.args.chrom)
            self._snps_masks_lists, self._compressed_masks = self.compress_cis_masks(cis_masks)

        if self.args.forcetrans:
            self.logger.debug("Forcing trans detection: removing genes from Chr {:d}".format(self.args.chrom))
            ix2keep = list()
            for i, g in enumerate(self._geneinfo):
                if g.chrom != self.args.chrom:
                    ix2keep.append(i)
            ix2keep = np.array(ix2keep)
            self._expr = self._expr[ix2keep, :]
            self._geneinfo = [self._geneinfo[i] for i in ix2keep]
        elif self.args.forcecis:
            self.logger.debug("Forcing cis detection: removing genes NOT from Chr {:d}".format(self.args.chrom))
            ix2keep = list()
            for i, g in enumerate(self._geneinfo):
                if g.chrom == self.args.chrom:
                    ix2keep.append(i)
            ix2keep = np.array(ix2keep)
            self._expr = self._expr[ix2keep, :]
            self._geneinfo = [self._geneinfo[i] for i in ix2keep]

        self.normalize_and_center_dosage(dosage_filtered_selected)

        #self._gtnorm = self._gtnorm[:, vcfmask]
        #self._gtcent = self._gtcent[:, vcfmask]

        self.logger.info("Completed data reading")

    def simulate(self):

        # Gene Expression
        if not self.args.gxsim:
            rpkm = ReadRPKM(self.args.gx_file, "gtex")
            expr = rpkm.expression
            gene_names = rpkm.gene_names
            geneinfo = list()
            for gene in gene_names:
                this_gene = GeneInfo(name = "x-gene", ensembl_id = gene, chrom = 1, start = 1, end   = 2)
                geneinfo.append(this_gene)
        else:
            expr = np.random.normal(0, 1, 2189 * 338).reshape(2189, 338)
            gene_names = ['ENS{:07d}'.format(i+1) for i in range(2189)]
            geneinfo = list()
            for gene in gene_names:
                this_gene = GeneInfo(name = "x-gene", ensembl_id = gene, chrom = 1, start = 1, end   = 2)
                geneinfo.append(this_gene)

        # Genotype
        fmin = float(self.args.simparams[0])
        fmax = float(self.args.simparams[1])
        nsnp = int(self.args.simparams[2])
        ntrans = int(self.args.simparams[3])
        cfrac = float(self.args.simparams[4])
        maketest = self.args.maketest
        snpinfo, gtnorm, gtcent = simulate.permuted_dosage(expr, nsnp = nsnp, fmin = fmin, fmax = fmax, maketest = maketest)
       
        # Trans-eQTL
        if ntrans > 0:
            newgx, hg2, nc = simulate.expression(gtnorm[-ntrans:, :], expr, cfrac = cfrac)
            expr = newgx
        
        self._gtnorm = gtnorm
        self._gtcent = gtcent
        self._snpinfo = snpinfo
        self._expr = expr
        self._geneinfo = geneinfo
