import numpy as np
from iotools import readrpkm
from iotools import simulate
from iotools.readOxford import ReadOxford
from iotools.readRPKM import ReadRPKM
from utils.containers import GeneInfo

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

    def normalize_and_center_dosage(self, dosage):
        f = [snp.maf for snp in self._snpinfo]
        f = np.array(f).reshape(-1, 1)
        self._gtnorm = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
        self._gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)

    def load(self):
        # Read Oxford File
        if self.args.oxf_file:
            oxf = ReadOxford(self.args.oxf_file, self.args.fam_file, self.args.startsnp, self.args.endsnp)
            dosage = oxf.dosage
            gt_donor_ids = oxf.samplenames
            self._snpinfo = oxf.snpinfo

        # Read VCF file
        if self.args.vcf_file:
            vcf = ReadVCF(self.args.vcf_file, self.args.startsnp, self.args.endsnp)
            dosage = vcf.dosage
            gt_donor_ids = vcf.donor_ids
            self._snpinfo = oxf.snpinfo

        self.normalize_and_center_dosage(dosage)

        # Gene Expression
        rpkm = ReadRPKM(self.args.gx_file, "gtex")
        expression = rpkm.expression
        expr_donors = rpkm.donor_ids
        gene_names = rpkm.gene_names

        # reorder donors gt and expr
        vcfmask, exprmask = self.select_donors(gt_donor_ids, expr_donors)

        self._expr = expression[:, exprmask]
        self._gtnorm = self._gtnorm[:, vcfmask]
        self._gtcent = self._gtcent[:, vcfmask]

        print ("Completed data reading")

    def simulate(self):

        maketest = self.args.maketest
        gx_file = self.args.gx_file
        sampleids, expr, gene_names = readrpkm.read_expression(gx_file)

        if self.args.simparams is not None:
            fmin = float(self.args.simparams[0])
            fmax = float(self.args.simparams[1])
            nsnp = int(self.args.simparams[2])
            snpinfo, gtnorm, gtcent = simulate.single_snp_permute(nsnp = nsnp, nsample = expr.shape[1], fmin = fmin, maketest = maketest)
        else:
            snpinfo, gtnorm, gtcent = simulate.single_snp_permute(maketest = maketest)

        self._expr = expr
        self._geneinfo = list()
        for i in range(expr.shape[0]):
            this_gene = GeneInfo(name = "x-gene",
                                 ensembl_id = str(i),
                                 chrom = 1,
                                 start = 1,
                                 end   = 2)
            self._geneinfo.append(this_gene)
        self._gtnorm = gtnorm
        self._gtcent = gtcent
        self._snpinfo = snpinfo
