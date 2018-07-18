import numpy as np
from iotools import simulate
from iotools import readgtf
from iotools.readOxford import ReadOxford
from iotools.readRPKM import ReadRPKM
from utils.containers import GeneInfo
import scipy.stats as ss

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

class Data():

    def __init__(self, args):
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
    def expression(self):
        return self._expr


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
    
    def HWEcheck(self, genotype):
        bins = [0.66, 1.33]
        digitised_geno = np.digitize(genotype, bins).tolist()
        observed_freqs = np.array([0]*3)
        observed_freqs[0] = digitised_geno.count(0)
        observed_freqs[1] = digitised_geno.count(1)
        observed_freqs[2] = digitised_geno.count(2)
        n = sum(observed_freqs)
        p_A = (2*observed_freqs[0] + observed_freqs[1])/(2*n)
        p_a = (2*observed_freqs[2] + observed_freqs[1])/(2*n)
        X2 = n * ((4*observed_freqs[0]*observed_freqs[2] - observed_freqs[1]**2)/((2*observed_freqs[0] + observed_freqs[1])*(2*observed_freqs[2] + observed_freqs[1])))**2
        pval = 1 - ss.chi2.cdf(X2, 1)
        return pval 
    
    def filter_snps(self, snpinfo, dosage):
        # Predixcan style filtering of snps
        newsnps = list()
        newdosage = list()
        for i, snp in enumerate(snpinfo):
            pos = snp.bp_pos
            refAllele = snp.ref_allele
            effectAllele = snp.alt_allele
            rsid = snp.varid
            maf = snp.maf
            # Skip non-single letter polymorphisms
            if len(refAllele) > 1 or len(effectAllele) > 1:
                continue
            # Skip ambiguous strands
            if SNP_COMPLEMENT[refAllele] == effectAllele:
                continue
            # Skip unknown RSIDs
            if rsid == '.':
                continue
            # Skip low MAF
            if not (maf >= 0.10 and maf <=0.90):
                continue
            if(self.HWEcheck(dosage[i]) < 0.000001):
                continue
            newsnps.append(snp)
            bins = [0.66, 1.33]
            newdosage.append(np.digitize(dosage[i], bins))
        return newsnps, np.array(newdosage)

    def normalize_and_center_dosage(self, dosage):
        f = [snp.maf for snp in self._snpinfo]
        f = np.array(f).reshape(-1, 1)
        self._gtnorm = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
        self._gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)

    def load(self):
        # Read Oxford File
        if self.args.oxf_file:
            oxf = ReadOxford(self.args.oxf_file, self.args.fam_file, self.args.startsnp, self.args.endsnp, isdosage=self.args.isdosage, data_columns=self.args.oxf_columns) #should be self.args.isdosage and self.args.oxf_columns
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
        self._snpinfo = snpinfo_filtered

        # Gene Expression
        rpkm = ReadRPKM(self.args.gx_file, "gtex")
        expression = rpkm.expression
        expr_donors = rpkm.donor_ids
        gene_names = rpkm.gene_names
        ### for GTEx ###
        if(self.args.isdosage):
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=False)
        ### for Cardiogenics ###
        else:
            gene_info = readgtf.gencode_v12(self.args.gtf_file, trim=True)
        self._geneinfo = gene_info

        # reorder donors gt and expr
        vcfmask, exprmask = self.select_donors(gt_donor_ids, expr_donors)
        genes, indices = self.select_genes(gene_info, gene_names)

        self._expr = expression[:, exprmask][indices, :]
        dosage_filtered_selected = dosage_filtered[:, vcfmask]

        self.normalize_and_center_dosage(dosage_filtered_selected)

        #self._gtnorm = self._gtnorm[:, vcfmask]
        #self._gtcent = self._gtcent[:, vcfmask]

        print ("Completed data reading")

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
