import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper

class JPA:

    def __init__(self, x, y, start, end, mpi = False, rank = 0, comm = None, ncore = 1):
        self.gt = x
        self.gx = y
        self.start = start
        self.end = end
        self.mpi = mpi
        self._pvals = None
        self._qscores = None
        self.rank = rank
        self.comm = comm
        self.ncore = ncore


    @property
    def pvals(self):
        return self._pvals


    @property
    def scores(self):
        return self._qscores


    def jpascore(self, pvals):
        p = np.sort(pvals)
        n = p.shape[0]
        kmax = min(100, n)
        krange = [i + 1 for i in range(kmax)]
        digamma_n1 = special.digamma(n + 1)
        z = - ( np.log(p[:kmax]) - (special.digamma(krange) - digamma_n1) )
        zsum = np.cumsum(z)
        res = np.max(zsum)
        return res


    def ccomp(self, geno, expr):
        _path = os.path.dirname(__file__)
        clib = np.ctypeslib.load_library('../lib/linear_regression.so', _path)
        cfstat = clib.fit
        cfstat.restype = ctypes.c_int
        cfstat.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                           np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                           ctypes.c_int,
                           ctypes.c_int,
                           ctypes.c_int,
                           np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
                          ]
    
        x = geno.reshape(-1,)
        y = expr.reshape(-1,)
        xsize = x.shape[0]
        nsnps = geno.shape[0]
        nsample = geno.shape[1]
        ngene = expr.shape[0]
        fstat = np.zeros(nsnps * ngene)
        success = cfstat(x, y, nsnps, ngene, nsample, fstat)
        pvals = 1 - stats.f.cdf(fstat, 1, nsample-2)
        pvals = pvals.reshape(nsnps, ngene)

        qscores = np.array([self.jpascore(pvals[i,:]) for i in range(nsnps)])

        return pvals, qscores


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of genotypes for sending to your slaves
            geno = mpihelper.split_genotype(self.gt, self.ncore)
            expr = self.gx
        else:
            geno = None
            expr = None
        
        slave_geno = self.comm.scatter(geno, root = 0)
        expr = self.comm.bcast(expr, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        pvals, qscores = self.ccomp(slave_geno, expr)

        pvals   = self.comm.gather(pvals,   root = 0)
        qscores = self.comm.gather(qscores, root = 0)
        
        if self.rank == 0:
            self._pvals = np.vstack(pvals)
            self._qscores = np.concatenate(qscores)
        else:
            assert qscores is None
            assert pvals   is None
        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            pvals, qscores = self.ccomp(self.gt, self.gx)
            self._pvals = pvals
            self._qscores = qscores
        return
