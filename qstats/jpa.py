import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper
from utils.logs import MyLogger

class JPA:

    def __init__(self, x, y, comm, rank, ncore, qcalc, masks):
        self.gt = x
        self.gx = y
        self._pvals = None
        self._qscores = None
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

    # def apply_mask(self, cismasks, snpblocks):
    #     # cismasks:  each list element is a mask for cis genes
    #     # snpblocks: each list element is a list of snps corresponding to the mask
    #     newqscores = list()
    #     ngene = self.gx.shape[0]
    #     for j, mask in enumerate(cismasks):
    #         for i in snpblocks[j]:
    #             pvals = self.pvals[i*ngene:(i+1)*ngene]
    #             allgenes = np.ones(ngene, dtype=bool)
    #             if mask.shape[0] != 0: allgenes[mask] = False
    #             masked_pvals = pvals[allgenes]
    #             newqscore = self.jpascore(masked_pvals)
    #             # pvals_crop = np.delete(self.pvals[i,:], mask, axis = 0)
    #             # if len(mask) > 0:
    #             #     self._pvals[i,:][np.array(mask)] = 1
    #             # print("Iterating on SNP {:d} with {:d}/{:d} pvals".format(i, len(pvals_crop), len(self.pvals[i,:])))
    #             # newqscore = self.jpascore(pvals_crop)
    #             newqscores.append(newqscore)
    #     self._qscores = np.array(newqscores)

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
        res = 1 - stats.f.cdf(fstat, 1, nsample-2)
        return res


    def slavejob(self, geno, expr, nmax, offset, masks, jpa = True, usemask = False):
        self.logger.debug('Rank {:d} calculating SNPs {:d} to {:d}'.format(self.rank, offset+1, nmax + offset))
        nsnps = geno.shape[0]
        ngene = expr.shape[0]
        pvals = self.clinreg(geno, expr, nmax)
        if jpa:
            # calculate JPA for each SNP (using ngene pvals)
            if usemask:
                qscores = np.array([self.masked_jpascore(pvals[i*ngene : (i+1)*ngene], masks[i]) for i in range(nsnps)])
            else:
                qscores = np.array([self.jpascore(pvals[i*ngene : (i+1)*ngene]) for i in range(nsnps)])
        else:
            qscores = np.zeros(nsnps)
        return pvals, qscores


    def mpicompute(self):
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
        pvals, qscores = self.slavejob(slave_geno, slave_expr, slave_nsnp, slave_offs, slave_gmasks, jpa = self.qcalc, usemask = self.usemask)

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

        if self.rank == 0:
            self._qscores = np.concatenate(qscores)
        else:
            assert qscores is None
            assert recvbuf is None

        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            pvals, qscores = self.slavejob(self.gt, self.gx, self.gt.shape[0], 0, self.masks, jpa = self.qcalc, usemask = self.usemask)
            self._pvals = pvals.reshape(self.gt.shape[0], self.gx.shape[0])
            self._qscores = qscores
        return
