import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from statsmodels.distributions.empirical_distribution import ECDF

from utils import mpihelper
from utils.logs import MyLogger

class JPA:

    def __init__(self, x, y, comm, rank, ncore, qcalc, masks, get_pvals = False, qnull_file = None, statmodel = 'zstat'):
        self.gt = x
        self.gx = y
        self._pvals = None
        self._qscores = None
        self._jpa_pvals = None
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.mpi = False
        if self.ncore > 1:
            self.mpi = True

        self.qcalc = qcalc
        self.masks = masks
        self.usemask = True if masks is not None else False
        self.get_empirical_pvals = False
        if get_pvals:
            self.get_empirical_pvals = True
        self.qnull_file = qnull_file
        self.statmodel = statmodel
        self.logger = MyLogger(__name__)


    @property
    def jpa_pvals(self):
        return self._jpa_pvals

    @property
    def pvals(self):
        return self._pvals

    @property
    def scores(self):
        return self._qscores


    def masked_jpascore(self, pvals, mask):
        usedgenes = np.ones_like(pvals, dtype=np.bool)
        if mask.shape[0] != 0: usedgenes[mask] = False
        res = self.jpascore(pvals[usedgenes])
        return res


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


    def clinreg_fstat(self, geno, expr, nrow):
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
        res = 1 - stats.f.cdf(fstat, 1, nsample-2)
        return res


    def clinreg_zstat(self, geno, expr, nrow):
        _path = os.path.dirname(__file__)
        clib = np.ctypeslib.load_library('../lib/linear_regression_zstat.so', _path)
        czstat = clib.fit
        czstat.restype = ctypes.c_int
        czstat.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
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
        zstat = np.zeros(nsnps * ngene)
        success = czstat(x, y, nsnps, ngene, nsample, zstat)
        res = 2.0 * (1 - stats.norm.cdf(np.abs(zstat)))
        return res 


    def get_qecdf_fit(self, q, ntop):
        qecdf = ECDF(q)
        qsort = q[np.argsort(q)]
        qneg = qsort[-ntop:]
        qcut = qneg[0]
        cumsum = 0
        for qnull in qneg:
            cumsum += qnull - qcut
        lam = (1 / ntop) * cumsum
        prefact = ntop / q.shape[0]
        return qecdf, qcut, lam, prefact


    def p_qscore(self, q, qecdf, qcut, lam, prefact):
        if q < qcut:
            res = 1 - qecdf(q)
        else:
            res = prefact * np.exp(- (q - qcut) / lam)
        return res


    def get_qnull(self):
        ''' 
        This function reads the qnull file, if provided. Otherwise, generates uniform random number.
        Must be called from master and broadcast to the slaves.
        '''
        qmod = np.array([])
        if self.get_empirical_pvals:
            if self.qnull_file is not None and os.path.isfile(self.qnull_file):
                #qnull = self.read_qnull(self.qnull_file)
                self.logger.debug("Read null Q-scores from {:s}".format(self.qnull_file))
                qnull = list()
                with open(self.qnull_file, 'r') as instream:
                    for line in instream:
                        lsplit = line.strip().split()
                        q = float(lsplit[0].strip())
                        qnull.append(q)
                qnull = np.array(qnull)
            else:
                self.logger.debug("Creating null Q-scores from uniform p-value distribution")
                ngene = 10000
                nsnps = 50000
                qnull = np.array([self.jpascore(np.random.uniform(0, 1, size = ngene)) for i in range(nsnps)])
            qmod = qnull[np.isfinite(qnull)]
            self.logger.debug("Obtained {:d} null Q-scores".format(qmod.shape[0]))
        return qmod


    def slavejob(self, geno, expr, nmax, offset, masks, qnull):
        '''
        Outputs.
            pvals: p-values from every SNP-gene pair linear regression. Dimension I x G.
            qscores: JPA-score from the above p-values. Dimension I.
            p_jpa: p-values for the significance of JPA, calculated empirically. Dimension I.
        '''
        self.logger.debug('Rank {:d} calculating SNPs {:d} to {:d}'.format(self.rank, offset+1, nmax + offset))
        nsnps = geno.shape[0]
        ngene = expr.shape[0]

        # Simple linear regression calculating either f-statistic or Z-statistic for every SNP-gene pair
        if self.statmodel == 'fstat':
            pvals = self.clinreg_fstat(geno, expr, nmax)
        elif self.statmodel == 'zstat':
            pvals = self.clinreg_zstat(geno, expr, nmax)

        # Calculate JPA-score
        if self.qcalc:
            # calculate JPA for each SNP (using ngene pvals)
            if self.usemask:
                qscores = np.array([self.masked_jpascore(pvals[i*ngene : (i+1)*ngene], masks[i]) for i in range(nsnps)])
            else:
                qscores = np.array([self.jpascore(pvals[i*ngene : (i+1)*ngene]) for i in range(nsnps)])
        else:
            qscores = np.zeros(nsnps)

        # Calculate empirical p-values for the JPA-scores
        if self.get_empirical_pvals:
            ntop = min(500, int(qnull.shape[0] / 10))
            qecdf, qcut, lam, prefact = self.get_qecdf_fit(qnull, ntop)
            p_jpa = np.array([self.p_qscore(q, qecdf, qcut, lam, prefact) for q in qscores])
        else:
            p_jpa = np.array([1 for q in qscores])

        return pvals, qscores, p_jpa


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of genotypes for sending to your slaves
            geno, offset = mpihelper.split_genotype(self.gt, self.ncore)
            expr = self.gx
            nsnp = [x.shape[0] for x in geno]
            gmasks = mpihelper.split_genemasks(self.masks, nsnp, offset)
            qnull = self.get_qnull()
        else:
            geno = None
            expr = None
            nsnp = None
            offset = None
            gmasks = None
            qnull  = None
            slave_geno = None
            slave_expr = None
            slave_nsnp = None
            slave_offs = None
            slave_gmasks = None
            slave_qnull  = None
        
        slave_geno = self.comm.scatter(geno, root = 0)
        slave_expr = self.comm.bcast(expr, root = 0)
        slave_nsnp = self.comm.scatter(nsnp, root = 0)
        slave_offs = self.comm.scatter(offset, root = 0)
        if self.usemask:
            slave_gmasks = self.comm.scatter(gmasks, root = 0)
        else:
            slave_gmasks = self.comm.bcast(gmasks, root = 0)
        slave_qnull = self.comm.bcast(qnull, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Do the calculations
        # ==================================
        pvals, qscores, p_jpa = self.slavejob(slave_geno, slave_expr, slave_nsnp, slave_offs, slave_gmasks, slave_qnull) 

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
        self.comm.Gatherv(sendbuf=pvals, recvbuf=(recvbuf, flat_sizes), root = 0)

        if self.rank == 0:
            self._pvals = recvbuf.reshape(sum(nsnp), ngene)

        qscores = self.comm.gather(qscores, root = 0)
        p_jpa = self.comm.gather(p_jpa, root = 0)

        if self.rank == 0:
            self._qscores = np.concatenate(qscores)
            if self.get_empirical_pvals:
                self._jpa_pvals = np.concatenate(p_jpa)
        else:
            assert qscores is None
            assert recvbuf is None
            assert p_jpa is None

        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            qnull = self.get_qnull()
            pvals, qscores, p_jpa = self.slavejob(self.gt, self.gx, self.gt.shape[0], 0, self.masks, qnull) 
            self._pvals = pvals.reshape(self.gt.shape[0], self.gx.shape[0])
            self._qscores = qscores
            if self.get_empirical_pvals:
                self._jpa_pvals = p_jpa
        return
