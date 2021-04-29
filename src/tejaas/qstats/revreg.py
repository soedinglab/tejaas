import numpy as np
import os
import time
import ctypes

from scipy import special
from scipy import stats
from scipy.optimize import minimize
#from statsmodels.distributions.empirical_distribution import ECDF

from tejaas.utils import mpihelper
from tejaas.qstats import rrstat
from tejaas.qstats import crrstat
from tejaas.utils.logs import MyLogger

class RevReg:

    def __init__(self, x, y, sigbeta2, comm, rank, ncore, null = 'perm', maf = None, masks = None):
        self.gt = x
        self.gx = y
        self.sigbeta2 = sigbeta2
        self.comm = comm
        self.rank = rank
        self.ncore = ncore
        self.null = null
        self.maf = maf
        self.mpi = False
        self.masks = masks
        self.usemask = False
        if self.masks is not None:
            self.usemask = True
        if self.ncore > 1:
            self.mpi = True
        self._pvals = None
        self._qscores = None
        self._mu = None
        self._sigma = None
        self._betas = None

        if self.null == 'perm':
            self.sigx2 = np.var(self.gt, axis = 1)
        elif self.null == 'maf':
            self.sigx2 = np.ones(self.gt.shape[0])

        self.logger = MyLogger(__name__)

    @property
    def sb2(self):
        return self.sigbeta2

    @sb2.setter
    def sb2(self, value):
        self.sigbeta2 = value

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

    @property
    def betas(self):
        return self._betas

    def slavejob(self, gt, gx, sb2, sx2, maf, masks, start, end, usemask, get_betas = False):
        if usemask:
            if len(masks) == 0: return [], [], [], [], np.array([])
            startsnp = min([min(x.apply2) for x in masks])
            endsnp = max([max(x.apply2) for x in masks])
            totsnp = sum(x.nsnp for x in masks)
            self.logger.debug("Rank {:d} using {:d} masks on {:d} SNPs [{:d} to {:d}]".format(self.rank, len(masks), totsnp, startsnp, endsnp))
            stime = time.time()
            p, q, mu, sig, b = self.maskjob(gt, gx, sb2, sx2, maf, masks, get_betas)
            self.logger.debug("Rank {:d} took {:g} seconds".format(self.rank, time.time() - stime))
        else:
            self.logger.debug("Rank {:d}. Using {:d} SNPs [{:d} to {:d}]".format(self.rank, end - start, start, end-1))
            applyon = np.arange(start,end)
            p, q, mu, sig, b = self.basejob(gt, gx, sb2, sx2, maf, applyon, get_betas)
        return p, q, mu, sig, b

    def maskjob(self, gt, gx, sb2, sx2, maf, masks, get_betas):
        p = np.array([])
        q = np.array([])
        mu = np.array([])
        sig = np.array([])
        b = np.array([])
        for mask in masks:
            usegenes = np.ones(gx.shape[0], dtype=bool)
            if mask.rmv_id.shape[0] > 0: usegenes[mask.rmv_id] = False
            masked_gx = np.ascontiguousarray(gx[usegenes])
            _p, _q, _mu, _sig, _b = self.basejob(gt, masked_gx, sb2, sx2, maf, np.array(mask.apply2), get_betas)
            p = np.append(p, _p)
            q = np.append(q, _q)
            mu = np.append(mu, _mu)
            sig = np.append(sig, _sig)

            if get_betas:
                # set beta value for masked genes to zero
                betas = self.reshape_masked_betas(_b, mask, gx.shape[0])
                b = np.append(b, betas)

        return p, q, mu, sig, b

    def basejob(self, gt, gx, sb2, sx2, maf, applyon, get_betas):
        slv_gt  = np.ascontiguousarray(gt[applyon, :])
        slv_gx  = gx
        slv_sb2 = sb2[applyon]
        slv_sx2 = sx2[applyon]
        b = []
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        elif self.null == 'maf':
            slv_maf = maf[applyon]
            p, q, mu, sig = crrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        if get_betas:
            b = crrstat.crrbetas(slv_gt, slv_gx, slv_sb2)
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig, b


    def reshape_masked_betas(self, b, mask, ngenes):
        self.logger.debug("Rank {:d}: reshaping {:d} betas into ({:d},{:d}) with {:d} masked genes out of {:d}".format(self.rank, len(b), len(mask.apply2), (ngenes - len(mask.rmv_id)), len(mask.rmv_id), ngenes)) 
        _b = b.reshape(len(mask.apply2), ngenes-len(mask.rmv_id))
        paddedBeta = np.zeros( (len(mask.apply2), ngenes) )
        inv_ind = np.delete(np.arange(ngenes), mask.rmv_id)
        paddedBeta[:, inv_ind] = _b
        return paddedBeta.reshape(-1)


    def mpicompute(self, get_betas = False):
        gmasks = None
        pstart = None
        pend = None
        geno = None
        expr = None
        sb2  = None
        sx2  = None
        maf  = None
        if self.rank == 0:
            # this is the master
            # create a list of index for sending to your slaves
            if self.usemask:
                self.logger.debug("Masks on: " + ", ".join(['{:d}'.format(x.nsnp) for x in self.masks]))
                gmasks = mpihelper.split_maskcomp(self.masks, self.ncore)
            else:
                pstart, pend = mpihelper.split_n(self.gt.shape[0], self.ncore)

            self.logger.debug("Get betas set to "+str(get_betas))
            geno = self.gt
            expr = self.gx
            sb2  = self.sigbeta2
            sx2  = self.sigx2
            maf  = self.maf

        sb2  = self.comm.bcast(sb2,  root = 0)
        sx2  = self.comm.bcast(sx2,  root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        expr = self.comm.bcast(expr, root = 0)
        geno = self.comm.bcast(geno, root = 0)
        if self.usemask:
            gmasks = self.comm.scatter(gmasks, root = 0)
        else:
            pstart = self.comm.scatter(pstart, root = 0)
            pend = self.comm.scatter(pend, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        pvals, qscores, mu, sigma, betas = self.slavejob(geno, expr, sb2, sx2, maf, gmasks, pstart, pend, self.usemask, get_betas = get_betas)

        pvals   = self.comm.gather(pvals,   root = 0)
        qscores = self.comm.gather(qscores, root = 0)
        mu      = self.comm.gather(mu,      root = 0)
        sigma   = self.comm.gather(sigma,   root = 0)

        if get_betas:
            recvbuf = None
            betalength = len(betas)
            self.comm.barrier()   # is it necessary?
            received_counts = self.comm.gather(betalength)
            if self.rank == 0:
                self.logger.debug("Number of coefficients from each node: {:s}".format(", ".join(['{:d}'.format(x) for x in received_counts])))
                recvbuf = np.zeros(np.sum(received_counts), dtype=np.float64)
            self.comm.Gatherv(sendbuf=betas, recvbuf=(recvbuf, received_counts), root = 0)

        if self.rank == 0:
            self._pvals   = np.concatenate(pvals)
            self._qscores = np.concatenate(qscores)
            self._mu      = np.concatenate(mu)
            self._sigma   = np.concatenate(sigma)
            if get_betas:
                self._betas   = recvbuf.reshape(self.gt.shape[0], self.gx.shape[0])
                self.logger.debug("All nodes computed a total of {:d} pvalues and {:s} betas".format(len(self._pvals), str(self._betas.shape)))
        else:
            assert qscores is None
            assert pvals   is None
            assert mu      is None
            assert sigma   is None
        return

    def compute(self, get_betas = False):
        if self.mpi:
            self.mpicompute(get_betas)
        else:
            start = 0
            end = self.gt.shape[0]
            self._pvals, self._qscores, self._mu, self._sigma, self._betas = self.slavejob(self.gt, self.gx,
                                                                              self.sigbeta2, self.sigx2, self.maf,
                                                                              self.masks, start, end, self.usemask, get_betas)
            if get_betas:
                self._betas = self._betas.reshape(self.gt.shape[0], self.gx.shape[0])
        return
