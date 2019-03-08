import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper
from utils.logs import MyLogger

class JPAFSTATS:

    def __init__(self, x, y, comm, rank, ncore):
        self.gt = x
        self.gx = y
        self._fstats = None
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.mpi = False
        if self.ncore > 1:
            self.mpi = True
        self.logger = MyLogger(__name__)


    @property
    def fstats(self):
        return self._fstats


    def clinreg(self, geno, expr, nrow):
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
        #res = 1 - stats.f.cdf(fstat, 1, nsample-2)
        return fstat


    def slavejob(self, geno, expr, nmax, offset):
        self.logger.debug('Rank {:d} calculating SNPs {:d} to {:d}'.format(self.rank, offset+1, nmax + offset))
        nsnps = geno.shape[0]
        ngene = expr.shape[0]
        fstat = self.clinreg(geno, expr, nmax)
        return fstat


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of genotypes for sending to your slaves
            geno, offset = mpihelper.split_genotype(self.gt, self.ncore)
            expr = self.gx
            nsnp = [x.shape[0] for x in geno]
        else:
            geno = None
            expr = None
            nsnp = None
            offset = None
            slave_geno = None
            slave_expr = None
            slave_nsnp = None
            slave_offs = None
        
        slave_geno = self.comm.scatter(geno, root = 0)
        slave_expr = self.comm.bcast(expr, root = 0)
        slave_nsnp = self.comm.scatter(nsnp, root = 0)
        slave_offs = self.comm.scatter(offset, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Do the calculations
        # ==================================
        fstat = self.slavejob(slave_geno, slave_expr, slave_nsnp, slave_offs)

        # ==================================
        # Collect the results
        # ==================================
        recvbuf = None
        if self.rank == 0:
            self.logger.debug("Number of SNPs sent to each slave: " + ", ".join(["{:d}".format(x) for x in nsnp]))
            ngene = self.gx.shape[0]
            flat_sizes = np.array([n * ngene for n in nsnp])
            recvbuf = np.zeros(sum(flat_sizes), dtype=np.float64)
        else:
            flat_sizes = None
        self.comm.Gatherv(sendbuf=fstat, recvbuf=(recvbuf, flat_sizes), root = 0)

        if self.rank == 0:
            self._fstats = recvbuf.reshape(sum(nsnp), ngene)
        else:
            assert qscores is None
            assert recvbuf is None

        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            fstats = self.slavejob(self.gt, self.gx, self.gt.shape[0], 0)
            self._fstats = fstats.reshape(self.gt.shape[0], self.gx.shape[0])
        return
