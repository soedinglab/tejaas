import numpy as np
from iotools import readrpkm
from iotools import simulate
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


    def load(self):
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
