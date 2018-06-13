import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
import logging
from utils import mpihelper
from qstats import rrstat
from qstats import crrstat
from utils.logs import MyLogger

class RevReg:


    def __init__(self, x, y, sigbeta2, comm, rank, ncore, null = 'perm', maf = None):
        self.gt = x
        self.gx = y
        self.sigbeta2 = sigbeta2
        self.comm = comm
        self.rank = rank
        self.ncore = ncore
        self.null = null
        self.maf = maf
        self.mpi = False
        if self.ncore > 1:
            self.mpi = True
        self._pvals = None
        self._qscores = None
        self._mu = None
        self._sigma = None

        if self.null == 'perm':
            self.sigx2 = np.var(self.gt, axis = 1)
        elif self.null == 'maf':
            self.sigx2 = np.ones(self.gt.shape[0])

        self.logger = MyLogger(__name__)

    @property
    def pvals(self):
        return self._pvals


    @property
    def scores(self):
        return self._qscores


    @property
    def null_mu(self):
        return self._mu


    @property
    def null_sigma(self):
        return self._sigma


#    def crevreg(self, geno, expr):
#        _path = os.path.dirname(__file__)
#        clib = np.ctypeslib.load_library('../lib/reverse_regression.so', _path)
#        cqscore = clib.qscore
#        cqscore.restype = ctypes.c_int
#        cqscore.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
#                            np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
#                            ctypes.c_int,
#                            ctypes.c_int,
#                            ctypes.c_int,
#                            np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
#                           ]
#    
#        x = geno.reshape(-1,)
#        y = expr.reshape(-1,)
#        xsize = x.shape[0]
#        nsnps = geno.shape[0]
#        nsample = geno.shape[1]
#        ngene = expr.shape[0]
#        fstat = np.zeros(nsnps * ngene)
#        success = cqscore(x, y, nsnps, ngene, nsample, fstat)
#        return pvals, qscores, mu, sigma


    def slavejob(self, gt, gx, sb2, sx2, maf, start, end):
        slv_gt  = gt[start:end, :]
        slv_gx  = gx
        slv_sb2 = sb2[start:end]
        slv_sx2 = sx2[start:end]
        if self.null == 'perm':
            p, q, mu, sig = rrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        elif self.null == 'maf':
            slv_maf = maf[start:end]
            p, q, mu, sig = rrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        return p, q, mu, sig


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of index for sending to your slaves
            nmax = int(self.gt.shape[0] / self.ncore)
            offset = 0
            for i in range(1, self.ncore):
                start = offset
                end   = offset + nmax
                self.comm.send(start, dest = i, tag = 10 + 3 * i - 2)
                self.comm.send(end,   dest = i, tag = 10 + 3 * i - 1)
                offset += nmax
            start = offset
            end = self.gt.shape[0]
        else:
            start = self.comm.recv(source = 0, tag = 10 + self.rank * 3 - 2)
            end   = self.comm.recv(source = 0, tag = 10 + self.rank * 3 - 1)

        if self.rank == 0:
            geno = self.gt
            expr = self.gx
            sb2  = self.sigbeta2
            sx2  = self.sigx2
            maf  = self.maf
        else:
            geno = None
            expr = None
            sb2  = None
            sx2  = None
            maf  = None
        
        sb2  = self.comm.bcast(sb2,  root = 0)
        sx2  = self.comm.bcast(sx2,  root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        expr = self.comm.bcast(expr, root = 0)
        geno = self.comm.bcast(geno, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        self.logger.debug("Reporting from node {:d}. Start: {:d} and End: {:d}".format(self.rank, start, end))
        pvals, qscores, mu, sigma = self.slavejob(geno, expr, sb2, sx2, maf, start, end)

        pvals   = self.comm.gather(pvals,   root = 0)
        qscores = self.comm.gather(qscores, root = 0)
        mu      = self.comm.gather(mu,      root = 0)
        sigma   = self.comm.gather(sigma,   root = 0)

        if self.rank == 0:
            self._pvals   = np.concatenate(pvals)
            self._qscores = np.concatenate(qscores)
            self._mu      = np.concatenate(mu)
            self._sigma   = np.concatenate(sigma)
        else:
            assert qscores is None
            assert pvals   is None
            assert mu      is None
            assert sigma   is None
        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            start = 0
            end = self.gt.shape[0]
            pvals, qscores, mu, sigma = self.slavejob(self.gt, self.gx, self.sigbeta2, self.sigx2, self.maf, start, end)
            self._pvals = pvals
            self._qscores = qscores
            self._mu = mu
            self._sigma = sigma
        return
