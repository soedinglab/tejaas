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
# from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import expon


from scipy.optimize import minimize

class SBoptimizer:

    def __init__(self, GT, GX, sx2):

        self._GT  = np.ascontiguousarray(GT)
        self._GX  = np.ascontiguousarray(GX)
        self._sx2 = np.ascontiguousarray(sx2)
        self._nsnps = GT.shape[0]
        self._nsample = GX.shape[1]
        
        U, S, VT = np.linalg.svd(GX.T)
        self._S = S
        self._U = U
        self._S2 = np.square(S)
        self._opt_sb2 = np.zeros(self._nsnps)
    
    @property
    def sb2(self):
        return self._opt_sb2

    def get_ML(self, _sb2, i):
        # sb2 = sb * sb
        sb2 = np.exp(_sb2)
        S2mod = self._S2 + (self._sx2[i] / sb2)
        Rscore = np.sum(np.square(np.dot(self._U.T, self._GT[i,:])) * (self._S2 / S2mod)) / self._sx2[i]
        MLL = -0.5*np.sum(np.log( self._S2 * (sb2 / self._sx2[i]) + 1 )) + 0.5*Rscore

        denom = (self._S2 * sb2 + self._sx2[i])
        der = 0.5* np.sum( ( self._S2 / denom ) * ( (np.square(np.dot(self._U.T, self._GT[i,:])) / denom ) - 1 ) )
        return -MLL, sb2*np.array([-der])

    def fit(self):
        st = time.time()
        
        sb_init = np.exp(0.01)
        for i in range(self._nsnps):
            res = minimize(   self.get_ML,
                              sb_init, 
                              args = i,
                              method='L-BFGS-B',
                              jac = True,
                              #bounds = [[0,1]],
                              options={'maxiter': 200000,
                                       'maxfun': 2000000,
                                       #'ftol': 1e-9,
                                       #'gtol': 1e-9,
                                       'disp': False})

            # print(res)
            self._opt_sb2[i] = np.exp(res.x[0])
        et = time.time()
        print("MML sb2 optimization took: ",et-st)

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
        self._selected_snps = np.array([])
        self._selected_genes = None
        self._null_qscores = None

        if self.null == 'perm' or self.null == "no_null":
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

    @property
    def selected_snps(self):
        return self._selected_snps

    @property
    def selected_genes(self):
        return self._selected_genes

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

            ## Optimize sb2 (if its None, then optimize it)
            """ if sb2[0] is None:
                optimizer = SBoptimizer(gt[mask.apply2], masked_gx, sx2[mask.apply2])
                optimizer.fit()
                sb2_opt = np.zeros(gt.shape[0])
                sb2_opt[mask.apply2] = optimizer.sb2 * 100 # increase 10-fold sigma_beta heuristic

                _p, _q, _mu, _sig, _b = self.basejob(gt, masked_gx, sb2_opt, sx2, maf, np.array(mask.apply2), get_betas)
            else: """
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
        if self.null == 'no_null':
            p, q, mu, sig = crrstat.no_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        elif self.null == 'maf':
            slv_maf = maf[applyon]
            p, q, mu, sig = crrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        if get_betas:
            b = crrstat.crrbetas(slv_gt, slv_gx, slv_sb2)
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig, b

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

    # def fit_exp(self, nulldistrib):
    #     halfnull = nulldistrib[nulldistrib >= np.mean(nulldistrib)]
    #     loc, scale = expon.fit(halfnull)
    #     lambd = 1 / scale
    #     return lambd

    # def p_qscore_old(self, score, random_scores, p0 = 0.001):
    #     n1 = int(p0 * len(random_scores))
    #     maxscore = random_scores[np.argsort(-random_scores)][n1]
    #     if score > maxscore:
    #         lambd = self.fit_exp(random_scores)
    #         res = np.exp(- lambd * score)
    #     else:
    #         ecdf = ECDF(random_scores)
    #         res = 1.0 - ecdf(score)
    #     return res

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
            print("received_counts ", received_counts)
            if self.rank == 0:
                recvbuf = np.zeros(np.sum(received_counts), dtype=np.float64)
            self.comm.Gatherv(sendbuf=betas, recvbuf=(recvbuf, received_counts), root = 0)

        if self.rank == 0:
            self._pvals   = np.concatenate(pvals)
            self._qscores = np.concatenate(qscores)
            self._mu      = np.concatenate(mu)
            self._sigma   = np.concatenate(sigma)
            if get_betas:
                self._betas   = recvbuf.reshape(self.gt.shape[0], self.gx.shape[0])
                self.logger.debug("Rank {:d}: Sparse: all nodes computed a total of {:d} pvalues and {:s} betas".format(self.rank, len(self._pvals), str(self._betas.shape)))
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

    def prune_masks(self, masks, snpi):
        selected_masks = [x for x in masks if len(list(set(x.apply2) & set(snpi)))]
        newmasks = []
        for m in selected_masks:
            pruned_apply2 = sorted(list(set(snpi) & set(m.apply2)))
            newmasks.append(m._replace(apply2=[snpi.index(x) for x in pruned_apply2]))
        return newmasks

    def select_best_genes(self, betas, n=1000, selected_genes=None):
        # This function selectes the best n beta values for each snp
        # it returns a matrix snp x n with the indices of the genes with highest effect size
        # selected_genes holds the matrix of a previous iteration of selected_best_genes
        gene_indices = np.array([], dtype=int)
        for j in range(betas.shape[0]):  # iterate over all snps
            beta_j        = np.abs(betas[j])
            bestbetas_ind = sorted(np.argpartition(beta_j, -n)[-n:])
            if selected_genes is not None:
                ix = selected_genes[j,:][bestbetas_ind]
                gene_indices = np.append(gene_indices, ix)
            else:
                gene_indices = np.append(gene_indices, bestbetas_ind)
        gene_indices = gene_indices.reshape(betas.shape[0], n)
        self.logger.debug("Rank {:d}: Selected best {:d} beta values for the best {:d} SNPs".format(self.rank, gene_indices.shape[1], gene_indices.shape[0]))
        return gene_indices

    def select_best_SNPs(self, n=1000, pval_thres = 1e-2, use_pvals = True):
        if use_pvals:
            bestsnps_ind = np.argwhere(self._pvals < pval_thres).reshape(-1,)
            # worstsnps_ind = np.argwhere(self._pvals > 0.8).reshape(-1,)[:2] # used to select worst pvals
            # bestsnps_ind = np.append(bestsnps_ind, worstsnps_ind)
            self.logger.debug("Selected {:d} pvalues below {:f} threshold".format(len(bestsnps_ind), pval_thres))
            # self.logger.debug("Selecting best {:d} pvalues out of {:d}".format(n, len(self.pvals)))
            # bestsnps_ind = np.argpartition(self.pvals, n)[:n]
        else:
            self.logger.debug("Selecting best {:d} Qscores out of {:d}".format(n, len(self.qscores)))
            bestsnps_ind = sorted(np.argpartition(self.qscores, -n)[-n:])
        return bestsnps_ind

