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
import time

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


    def slavejob(self, gt, gx, sb2, sx2, maf, masks, start, end, usemask):
        if usemask:
            startsnp = masks[0].apply2[0]
            endsnp = masks[-1].apply2[-1]
            totsnp = endsnp - startsnp + 1
            self.logger.debug("Rank {:d} using {:d} masks on {:d} SNPs [{:d} to {:d}]".format(self.rank, len(masks), totsnp, startsnp, endsnp))
            stime = time.time()
            p, q, mu, sig = self.maskjob(gt, gx, sb2, sx2, maf, masks)
            self.logger.debug("Rank {:d} took {:g} seconds".format(self.rank, time.time() - stime))
        else:
            self.logger.debug("Reporting from node {:d}. Using {:d} SNPs [{:d} to {:d}]".format(self.rank, end - start + 1, start, end))
            p, q, mu, sig = self.basejob(gt, gx, sb2, sx2, maf, start, end)
        return p, q, mu, sig


    def maskjob(self, gt, gx, sb2, sx2, maf, masks):
        p = np.array([])
        q = np.array([])
        mu = np.array([])
        sig = np.array([])
        for mask in masks:
            start = min(mask.apply2)
            end = max(mask.apply2) + 1
            usegenes = np.ones(gx.shape[0], dtype=bool)
            if mask.rmv_id.shape[0] > 0: usegenes[mask.rmv_id] = False
            masked_gx = np.ascontiguousarray(gx[usegenes])
            _p, _q, _mu, _sig = self.basejob(gt, masked_gx, sb2, sx2, maf, start, end)
            p = np.append(p, _p)
            q = np.append(q, _q)
            mu = np.append(mu, _mu)
            sig = np.append(sig, _sig)
        return p, q, mu, sig


    def basejob(self, gt, gx, sb2, sx2, maf, start, end):
        slv_gt  = np.ascontiguousarray(gt[start:end, :])
        slv_gx  = gx
        slv_sb2 = sb2[start:end]
        slv_sx2 = sx2[start:end]
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        elif self.null == 'maf':
            slv_maf = maf[start:end]
            p, q, mu, sig = crrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig


    def mpicompute(self):
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
                self.logger.debug("Masks on: " + ", ".join(['{:d}'.format(len(x.apply2)) for x in self.masks]))
                gmasks = mpihelper.split_maskcomp(self.masks, self.ncore)
            else:
                pstart, pend = mpihelper.split_n(self.gt.shape[0], self.ncore)
            geno = self.gt
            expr = self.gx
            sb2  = self.sigbeta2
            sx2  = self.sigx2
            maf  = self.maf

        ##     offset = nmax # this is a hack to keep things in order. Node 0 must get 0:nmax
        ##     for i in range(1, self.ncore):
        ##         start = offset
        ##         end   = offset + nmax if i < (self.ncore - 1) else self.gt.shape[0]
        ##         self.comm.send(start, dest = i, tag = 10 + 3 * i - 2)
        ##         self.comm.send(end,   dest = i, tag = 10 + 3 * i - 1)
        ##         offset += nmax
        ##     start = 0
        ##     end = nmax
        ## else:
        ##     start = self.comm.recv(source = 0, tag = 10 + 3 * self.rank - 2)
        ##     end   = self.comm.recv(source = 0, tag = 10 + 3 * self.rank - 1)

        ## if self.rank == 0:
        ## else:
        
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
        pvals, qscores, mu, sigma = self.slavejob(geno, expr, sb2, sx2, maf, gmasks, pstart, pend, self.usemask)

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
            self._pvals, self._qscores, self._mu, self._sigma = self.slavejob(self.gt, self.gx, 
                                                                              self.sigbeta2, self.sigx2, self.maf, 
                                                                              self.masks, start, end, self.usemask)
        return
