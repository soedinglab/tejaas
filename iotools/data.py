import numpy as np
import scipy.stats as ss
import random
import os

from iotools import simulate
from iotools import readgtf
from iotools.readOxford import ReadOxford
from iotools.readvcf import ReadVCF
from iotools.readRPKM import ReadRPKM
from utils.containers import GeneInfo, CisMask
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


    @timeit
    def get_cismasklist_old(self, snpinfo, geneinfo, chrom, window=1e6):
        chr_genes_ix = np.array([i for i, g in enumerate(geneinfo) if g.chrom == chrom])
        chr_genes = [geneinfo[ix] for ix in chr_genes_ix]
        genemasks = list()
        iprev = 0
        for snp in snpinfo:
            pos = snp.bp_pos
            left = pos - window
            right = pos + window
            thismask = list()
            for i, g in enumerate(chr_genes[iprev:]):
                gstart = g.start
                gend = g.end
                if gstart >= left and gstart < right:
                    thismask.append(iprev + i)
                elif gend >= left and gend < right:
                    thismask.append(iprev + i)
                if gstart > right:
                    break
            if len(thismask) > 0:
                genemasks.append(chr_genes_ix[np.array(thismask)])
                iprev = thismask[0]
            else:
                genemasks.append(np.array([]))
        return genemasks


    @timeit
    def get_cismasklist(self, snpinfo, geneinfo, chrom, window=1e6):
        chr_genes_ix = [[] for ichrm in range(22)] 
        chr_genes = [[] for ichrm in range(22)]
        if chrom is not None:
            chr_genes_ix[chrom - 1] = np.array([i for i, g in enumerate(geneinfo) if g.chrom == chrom])
            chr_genes[chrom - 1] = [geneinfo[ix] for ix in chr_genes_ix[chrom - 1]]
        else:
            for ichrm in range(22):
                chr_genes_ix[ichrm] = np.array([i for i, g in enumerate(geneinfo) if g.chrom == ichrm + 1])
                chr_genes[ichrm] = [geneinfo[ix] for ix in chr_genes_ix[ichrm]]
        genemasks = list()
        iprev = 0
        ichrmprev = 0
        for snp in snpinfo:
            pos = snp.bp_pos
            left = pos - window
            right = pos + window
            ichrm = chrom - 1 if chrom is not None else snp.chrom - 1
            iprev_started = False
            if ichrm != ichrmprev:
                iprev = 0
                ichrmprev = ichrm
            thismask = list()
            for i, g in enumerate(chr_genes[ichrm][iprev:]):
                gstart = g.start
                gend = g.end
                if gstart >= left and gstart <= right:
                    # thismask.append(iprev + i)
                    thismask.append(chr_genes_ix[ichrm][iprev + i])
                    if not iprev_started:
                        new_start_iloc = iprev
                        iprev_started = True
                elif gend >= left and gend <= right:
                    # thismask.append(iprev + i)
                    thismask.append(chr_genes_ix[ichrm][iprev + i])
                    if not iprev_started:
                        new_start_iloc = iprev
                        iprev_started = True
                if gstart > right:
                    break
            if len(thismask) > 0:
                #genemasks.append(chr_genes_ix[np.array(thismask)])
                #iprev = thismask[0]
                genemasks.append(np.array(thismask))
                iprev = new_start_iloc
            else:
                genemasks.append(np.array([]))
        return genemasks


    @timeit
    def compress_cismasklist_old(self, genemasks):
        cismasks = list()
        appendmask = False
        setprev = False
        snplist = list()
        for i, mask in enumerate(genemasks):
            if not setprev: 
                prev_mask = mask
                setprev = True
            if np.all(np.array_equal(mask, prev_mask)):
                snplist.append(i)
            else:
                appendmask = True

            if i == len(genemasks) - 1: appendmask = True # no more masks to process

            if appendmask:
                thismask = CisMask(rmv_id = prev_mask, apply2 = snplist)
                cismasks.append(thismask)
                snplist = list([i])
                prev_mask = mask
                appendmask = False
        return cismasks


    @timeit
    def compress_cismasklist(self, genemasks):
        cismasks = list()
        appendmask = False
        endmask = False
        setprev = False
        snplist = list()
        for i, mask in enumerate(genemasks):
            if not setprev:
                prev_mask = mask
                setprev = True
            if np.all(np.array_equal(mask, prev_mask)):
                snplist.append(i)
            else:
                appendmask = True

            if i == len(genemasks) - 1: endmask = True # no more masks to process

            if appendmask:
                thismask = CisMask(rmv_id = prev_mask, apply2 = snplist)
                cismasks.append(thismask)
                snplist = list([i])
                prev_mask = mask
                if not endmask:
                    appendmask = False

            if endmask:
                # if not appendmask:
                #     snplist.append(i)
                thismask = CisMask(rmv_id = mask, apply2 = snplist)
                cismasks.append(thismask)

        return cismasks


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
            maf = round(snp.maf, 3)
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
            # hwep = self.HWEcheck(intdosage)
            #if(hwep < 0.000001):
            #    nhwep += 1
            #    # self.logger.debug("SNP {:s} has a HWE p-value of {:g}".format(rsid, hwep))
            #    continue
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
            snpinfo = vcf.snpinfo

        snpinfo_filtered, dosage_filtered = self.filter_snps(snpinfo, dosage)
        self.logger.debug("{:d} SNPs after filtering".format(len(snpinfo_filtered)))
        self._snpinfo = snpinfo_filtered

        # Gene Expression
        self.logger.debug("Reading expression levels")
        rpkm = ReadRPKM(self.args.gx_file, "gtex")
        expression = rpkm.expression
        expr_donors = rpkm.donor_ids
        gene_names = rpkm.gene_names

        if self.args.selected_donors:
            self.logger.debug("Selecting samples from user supplied list")
            with open(self.args.selected_donors) as instream:
                donors_user_list = [l.strip() for l in instream.readlines()]
            donors_ix   = np.array([expr_donors.index(i.strip()) for i in donors_user_list if i in expr_donors])
            expr_donors = [expr_donors[ix] for ix in donors_ix]
            expression  = rpkm._normalize_expr(expression[:, donors_ix])

        self.logger.debug("Found {:d} genes of {:d} samples".format(expression.shape[0], expression.shape[1]))
        self.logger.debug("Reading gencode file for gene information")

        if(self.args.gxtrim):
            ### for Cardiogenics ###
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=True)
        else:
            ### for GTEx ###
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=False)

        # Shuffle genotype?
        if self.args.shuffle:
            if os.path.isfile(self.args.shuffle_file):
                self.logger.warn("Shuffling genotype using supplied donor IDs")
                gt_donor_ids = [line.strip() for line in open(self.args.shuffle_file)]
            else:
                logger.warn("Shuffling genotype randomly")
                random.shuffle(gt_donor_ids)

        # reorder donors gt and expr
        self.logger.debug("Selecting common samples of genotype and gene expression")
        vcfmask, exprmask = self.select_donors(gt_donor_ids, expr_donors)
        genes, indices = self.select_genes(gene_info, gene_names)

        self._expr = expression[:, exprmask][indices, :]
        self._geneinfo = genes
        dosage_filtered_selected = dosage_filtered[:, vcfmask]

        self.logger.debug("Retained {:d} samples".format(vcfmask.shape[0]))

        ### Until here, all filters have been applied and geneinfo and snpinfo reflect current data ###
        # This holds!! geneinfo is ordered by chrom and start position
        # just for testing
        # for i, g in enumerate(genes):
        #     if g.start > genes[i+1].start and g.chrom == genes[i+1].chrom:
        #         print(g, genes[i+1])

        if self.args.forcetrans:
            self.logger.debug("Forcing trans detection: removing genes from Chr {:d}".format(self.args.chrom))
            ix2keep = np.array([i for i, g in enumerate(self._geneinfo) if g.chrom != self.args.chrom])
            self._expr = self._expr[ix2keep, :]
            self._geneinfo = [self._geneinfo[i] for i in ix2keep]
        elif self.args.forcecis:
            self.logger.debug("Forcing cis detection: removing genes NOT from Chr {:d}".format(self.args.chrom))
            ix2keep = np.array([i for i, g in enumerate(self._geneinfo) if g.chrom == self.args.chrom])
            self._expr = self._expr[ix2keep, :]
            self._geneinfo = [self._geneinfo[i] for i in ix2keep]

        if self.args.cismasking:
            self.logger.debug("Generate cis-masks for GX matrix for each SNP")
            #cis_masks = self.get_cismaskcomp(self._snpinfo, self._geneinfo, self.args.chrom)
            #snps_masks_lists, compressed_masks = self.compress_cis_masks(cis_masks)
            self._cismasklist = self.get_cismasklist(self._snpinfo, self._geneinfo, self.args.chrom, window=self.args.window)
            self._cismaskcomp = self.compress_cismasklist(self._cismasklist)
            
        self.normalize_and_center_dosage(dosage_filtered_selected)

        #self._gtnorm = self._gtnorm[:, vcfmask]
        #self._gtcent = self._gtcent[:, vcfmask]

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
