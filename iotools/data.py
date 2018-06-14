import numpy as np
from iotools import readrpkm
from iotools import simulate


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
            fmin = self.args.simparams[0]
            fmax = self.args.simparams[1]
            nsnp = self.args.simparams[2]
            snpinfo, gtnorm, gtcent = simulate.single_snp_permute(nsnp = nsnp, nsample = expr.shape[1], fmin = fmin, maketest = maketest)
        else:
            snpinfo, gtnorm, gtcent = simulate.single_snp_permute(maketest = maketest)

        self._expr = expr
        self._geneinfo = 0
        self._gtnorm = gtnorm
        self._gtcent = gtcent
        self._snpinfo = snpinfo
