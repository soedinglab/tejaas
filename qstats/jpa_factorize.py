import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper
from utils.logs import MyLogger

class CPMA:

    def __init__(self, x, y, comm, rank, ncore, qcalc, masks):
        self.gt = x
        self.gx = y
        self._pvals = None
        self._qscores = None
        self._fstat = None
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.qcalc = qcalc
        self.mpi = False
        self.masks = masks
        self.usemask = True if masks is not None else False
        if self.ncore > 1:
            self.mpi = True
        self.logger = MyLogger(__name__)


    @property
    def pvals(self):
        return self._pvals

    @property
    def fstats(self):
        return self._fstat

    @property
    def scores(self):
        return self._qscores

    # @scores.setter
    # def scores(self, scores):
    #     self._qscores = scores

    def masked_jpascore(self, pvals, mask):
        usedgenes = np.ones_like(pvals, dtype=np.bool)
        if mask.shape[0] != 0: usedgenes[mask] = False
        res = self.jpascore(pvals[usedgenes])
        return res

    ##### Set of functions for calculating p-values from empirical null distrib
    ##### by fitting exp to empirical CDF
    def p_ecdf(self, score, random_scores):
        ecdf = ECDF(random_scores)
        res = 1 - ecdf(score)
        return res

    def p_qnullcdf(self, score, nullscores):
        ntop = min(100, int(len(nullscores)/10))
        Q_neg = nullscores[np.argsort(nullscores)[-ntop:]]
        cumsum = 0
        for i in range(ntop):
            cumsum += Q_neg[i] - Q_neg[0]
        lam = (1/ntop) * cumsum
        pval = (ntop / len(nullscores)) * np.exp(- (score - Q_neg[0])/lam)
        return pval

    def p_qscore(self, score, nullscores):
        if score >= np.sort(nullscores)[-100]:
            return self.p_qnullcdf(score, nullscores)
        else:
            return self.p_ecdf(score, nullscores)
    ##### end of pval functions

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
        return fstat

    # Given our true f-stat, generate the null Qscore distribution
    def null_CPMA(self, W, Q, permutations = 100):
        ### Generate Z* here (would be our fstat*)
        # Z = np.mean(self._fstat) + etc ... 
        #
        ### Get pvals
        # pvals = 1 - stats.f.cdf(Z, 1, nsample-2)
        ### Calc qscores from pvals
        null_qscoress = self.slavejob(Z, Z.shape[0], Z.shape[1], self.masks, self.usemask)
        return null_qscores

    # Calculate f-stat for the original data
    def slavejob_fstat(self, geno, expr, nmax, offset):
        self.logger.debug('Rank {:d} calculating SNPs {:d} to {:d}'.format(self.rank, offset+1, nmax + offset))
        nsnps = geno.shape[0]
        ngene = expr.shape[0]
        fstat = self.clinreg(geno, expr, nmax)
        return fstat

    # Calculate Q-scores for any f-stat data (real or simulated like Z*)
    def slavejob(self, fstat, nsnps, ngene, masks, usemask = False):
        pvals = 1 - stats.f.cdf(fstat, 1, nsample-2)
        if usemask:
            qscores = np.array([self.masked_jpascore(pvals[i*ngene : (i+1)*ngene], masks[i]) for i in range(nsnps)])
        else:
            qscores = np.array([self.jpascore(pvals[i*ngene : (i+1)*ngene]) for i in range(nsnps)])
        return qscores

    def mpicompute_pval(self):
        if self.rank == 0:
            nsnp = [x.shape[0] for x in geno]
            gmasks = mpihelper.split_genemasks(self.masks, nsnp, offset)
        else:
            gmasks = None
            slave_nsnp = None
            slave_offs = None
            slave_gmasks = None

        # if self.usemask:
        #     slave_gmasks = self.comm.scatter(gmasks, root = 0)
        # else:
        #     slave_gmasks = self.comm.bcast(gmasks, root = 0)

        if rank == 0:
            self._qscores = self.slavejob(self, self._fstat, self._fstat.shape[0], self._fstat.shape[1], self.masks, self.usemask)

        if rank == 0:
            C = np.cov(self._fstat.T)
            W, Q = np.linalg.eig(C)
        else:
            W = None
            Q = None

        slave_W = self.comm.bcast(W, root = 0)
        slave_Q = self.comm.bcast(Q, root = 0)
        self.comm.barrier()

        null_qscores = self.null_CPMA(slave_W, _slave_Q, permutations = 100)
        self.comm.barrier()

        # you will recieve NCORE x permutations qscores
        null_qscores = self.comm.gather(qscores, root = 0)

        if rank == 0:
            self._pvals = [self.p_qscore(q, null_qscores) for q in self._qscores]

        return

    def mpicompute_fstat(self):
        if self.rank == 0:
            # this is the master
            # create a list of genotypes for sending to your slaves
            geno, offset = mpihelper.split_genotype(self.gt, self.ncore)
            expr = self.gx
            nsnp = [x.shape[0] for x in geno]
            gmasks = mpihelper.split_genemasks(self.masks, nsnp, offset)
        else:
            geno = None
            expr = None
            nsnp = None
            offset = None
            gmasks = None
            slave_geno = None
            slave_expr = None
            slave_nsnp = None
            slave_offs = None
            slave_gmasks = None
        
        slave_geno = self.comm.scatter(geno, root = 0)
        slave_expr = self.comm.bcast(expr, root = 0)
        slave_nsnp = self.comm.scatter(nsnp, root = 0)
        slave_offs = self.comm.scatter(offset, root = 0)
        if self.usemask:
            slave_gmasks = self.comm.scatter(gmasks, root = 0)
        else:
            slave_gmasks = self.comm.bcast(gmasks, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Do the calculations
        # ==================================
        fstat = self.slavejob_fstat(slave_geno, slave_expr, slave_nsnp, slave_offs, slave_gmasks, jpa = self.qcalc, usemask = self.usemask)

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
            self._fstat = recvbuf.reshape(sum(nsnp), ngene)
        return
            
    def compute_fstat(self):
        if self.mpi:
            self.mpicompute_fstat()
        else:
            fstat = self.slavejob_fstat(self.gt, self.gx, self.gt.shape[0], 0, self.masks, jpa = self.qcalc, usemask = self.usemask)
            self._fstat = fstat.reshape(self.gt.shape[0], self.gx.shape[0])
        return

    def compute_pval(self):
        if self.mpi:
            self.mpicompute_pval()
        else:
            ### broken, needs to be finished
            # C = np.cov(self._fstat.T)
            # W, Q = np.linalg.eig(C)
        return
