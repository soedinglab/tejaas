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
from statsmodels.distributions.empirical_distribution import ECDF

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

    def p_qscore(self, score, random_scores):
        ecdf = ECDF(random_scores)
        res = 1 - ecdf(score)
        return res

    def basejob_sparse(self, gt, gx, sb2, sx2, maf, get_betas):
        slv_gt  = np.ascontiguousarray(gt)
        slv_gx  = gx
        slv_sb2 = sb2
        slv_sx2 = sx2
        b = []
        stime = time.time()
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
            qarray = np.array([])
            self.logger.debug("Rank {:d}: Sparse: {:f}, {:f}, {:f}, {:f}".format(self.rank, p[0], q[0], mu[0], sig[0]))
            for i in range(1000):
                np.random.shuffle(slv_gx.T)
                _p, _q, _mu, _sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
                qarray = np.append(qarray, _q)
            p   = self.p_qscore(q, qarray)
            mu  = np.mean(qarray)
            sig = np.std(qarray)
            self.logger.debug("Rank {:d}: Sparse: newpval = {:f}".format(self.rank, p[0]))
            self.logger.debug("Rank {:d}: Sparse: took {:g} seconds for shuffling {:s} snps".format(self.rank, time.time() - stime, str(slv_gt.shape)))
        elif self.null == 'maf':
            slv_maf = maf
            p, q, mu, sig = crrstat.maf_null(slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        if get_betas:
            b = crrstat.crrbetas(slv_gt, slv_gx, slv_sb2)
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig, b

    # def slavejob_sparse(self, slave_geno, expr, sb2, slave_sx2, maf, slave_betas, nbetas, selected_genes = None):
    def slavejob_sparse(self, gt, gx, sb2, sx2, maf, gene_indices, start, end, get_betas):
        self.logger.debug("Rank {:d}: Using {:d} SNPs [{:d} to {:d}]".format(self.rank, end - start, start, end))
        p  = np.array([])
        q  = np.array([])
        mu = np.array([])
        s  = np.array([])
        b  = np.array([])
        stime = time.time()
        for snp_i in range(start, end): #range(gene_indices.shape[0]):  # iterate over all snps
            local_expr = gx[gene_indices[snp_i],:]
            snpgt      = gt[snp_i,:][np.newaxis].copy()
            _p, _q, _mu, _s, _b = self.basejob_sparse(snpgt, local_expr, sb2[snp_i].reshape(-1,), sx2[snp_i].reshape(-1,), maf[snp_i].reshape(-1), get_betas)
            p  = np.append(p, _p)
            q  = np.append(q, _q)
            mu = np.append(mu, _mu)
            s  = np.append(s, _s)
            b  = np.append(b, _b)
        self.logger.debug("Rank {:d} Sparse: Computed {:d} pvals and {:s} betas in {:g} seconds".format(self.rank, len(p), str(b.shape), time.time() - stime))
        return p, q, mu, s, b

    # def write_rr_out(self, outdir, suffix, snpinfo, geneinfo):
    #     mysnpinfo = snpinfo
    #     if len(self._selected_snps):
    #         mysnpinfo = [snpinfo[int(i)] for i in self._selected_snps]
    #     fname = outdir + "_rr_" + suffix + ".txt"
    #     with open(fname, "w") as f:
    #         f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
    #         print("write rr out: ", len(mysnpinfo), len(self.pvals))
    #         for i, x in enumerate(mysnpinfo):
    #             f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, self.scores[i], self.null_mu[i], self.null_sigma[i], self.pvals[i]))
    #     # np.savetxt(outdir + "_betas_" + suffix + ".txt", self.betas, fmt='%1.4e')
    #     if self._selected_genes is not None:
    #         np.savetxt(outdir + "_selected_genes_" + suffix + ".txt", self._selected_genes, fmt='%i')

    def mpicompute_sparse(self, get_betas = False):
        gt    = None
        expr  = None
        sx2   = None
        sb2   = None
        maf   = None
        gene_indices = None
        pstart = None
        pend = None
        finalize = False
        if self.rank == 0:
            # get best SNPs
            best_snp_indices = self.select_best_SNPs()
            if len(best_snp_indices) > 0:
                self._selected_snps = best_snp_indices
                # get best beta values for EACH snp
                gene_indices = self.select_best_genes(self.betas[best_snp_indices,:], n=1000)
                self._selected_genes = gene_indices
                gt   = self.gt[best_snp_indices,:]
                expr = self.gx
                sx2  = self.sigx2[best_snp_indices]
                sb2  = self.sigbeta2[best_snp_indices]
                maf  = self.maf
                pstart, pend = mpihelper.split_n(gene_indices.shape[0], self.ncore)
            else:
                self.logger.info("No significant SNPs available. No calculations to make.")
                finalize = True

        finalize = self.comm.bcast(finalize, root = 0)
        if finalize:
            return

        gt   = self.comm.bcast(gt, root = 0)
        expr = self.comm.bcast(expr, root = 0)
        sb2  = self.comm.bcast(sb2, root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        sx2  = self.comm.bcast(sx2, root = 0)
        gene_indices = self.comm.bcast(gene_indices, root = 0)
        pstart = self.comm.scatter(pstart, root = 0)
        pend   = self.comm.scatter(pend, root = 0)
        self.comm.barrier()

        ### Start sparse iteration of RR for each SNP ###
        self.logger.debug("Rank {:d}: Sparse: Calculating RR iteration over best {:d} genes (biggest beta values)".format(self.rank, gene_indices.shape[1]))
        self.logger.debug("Rank {:d}: Sparse: I have {:s} gene_indices, {:d} snps, geno is {:s} and {:s} sx2.".format(self.rank, str(gene_indices.shape), gt.shape[0], str(gt.shape), str(sx2.shape)))

        pvals, qscores, mu, sigma, betas = self.slavejob_sparse(gt, expr, sb2, sx2, maf, gene_indices, pstart, pend, get_betas)  

        pvals   = self.comm.gather(pvals,   root = 0)
        qscores = self.comm.gather(qscores, root = 0)
        mu      = self.comm.gather(mu,      root = 0)
        sigma   = self.comm.gather(sigma,   root = 0)

        if get_betas:
            recvbuf = None
            betalength = len(betas)
            self.comm.barrier()   # is it necessary?
            received_counts = self.comm.gather(betalength)
            self.logger.debug("Rank {:d}: Sparse: I have {:d} betas".format(self.rank, len(betas)))
            if self.rank == 0:
                recvbuf = np.zeros(np.sum(received_counts), dtype=np.float64)
            self.comm.Gatherv(sendbuf=betas, recvbuf=(recvbuf, received_counts), root = 0)

        if self.rank == 0:
            self._pvals   = np.concatenate(pvals)
            self.logger.debug("Rank {:d}: Sparse: Overwriting {:d} p-values".format(self.rank, len(self._pvals)))
            self._qscores = np.concatenate(qscores)
            self._mu      = np.concatenate(mu)
            self._sigma   = np.concatenate(sigma)
            if get_betas:
                self._betas   = recvbuf.reshape(gene_indices.shape)
                self.logger.debug("Rank {:d}: Sparse: all nodes computed a total of {:d} pvalues, {:s} betas, and {:s} selected genes".format(self.rank, len(self._pvals), str(self._betas.shape), str(gene_indices.shape)))
        else:
            assert qscores is None
            assert pvals   is None
            assert mu      is None
            assert sigma   is None
        return

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

    def compute_sparse(self, get_betas = False):
        if self.mpi:
            self.mpicompute_sparse(get_betas)
        else:
            start = 0
            end = self.gt.shape[0]
            best_snp_indices = self.select_best_SNPs()
            if len(best_snp_indices) == 0:
                self.logger.info("No significant SNPs available. No calculations to make.")
            else:
                self._selected_snps = best_snp_indices
                # get best beta values for EACH snp
                gene_indices = self.select_best_genes(self.betas[best_snp_indices,:], n=5000)
                self._selected_genes = gene_indices

                gt   = self.gt[best_snp_indices,:]
                expr = self.gx
                sx2  = self.sigx2[best_snp_indices]
                sb2  = self.sigbeta2[best_snp_indices]
                maf  = self.maf

                ### Start sparse iteration of RR for each SNP ###
                self.logger.debug("Sparse: Calculating RR iteration over best {:d} genes (biggest beta values)".format(gene_indices.shape[1]))
                self.logger.debug("Sparse: I have {:s} gene_indices, {:d} snps, geno is {:s} and {:s} sx2.".format(str(gene_indices.shape), gt.shape[0], str(gt.shape), str(sx2.shape)))

                self._pvals, self._qscores, self._mu, self._sigma, betas = self.slavejob_sparse(gt, self.gx, sb2, sx2, maf, gene_indices, start, end, get_betas)  
                if get_betas:
                    self._betas = betas.reshape(gene_indices.shape)
                self.logger.debug("Rank {:d}: Sparse: all nodes computed a total of {:d} pvalues, {:s} betas, and {:s} selected genes".format(self.rank, len(self._pvals), str(self._betas.shape), str(gene_indices.shape)))
        return

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
        self.logger.debug("Selected best {:d} beta values for the best {:d} SNPs".format(gene_indices.shape[1], gene_indices.shape[0]))
        return gene_indices

    def select_best_SNPs(self, n=1000, pval_thres = 1e-2, use_pvals = True):
        if use_pvals:
            bestsnps_ind = np.argwhere(self._pvals < pval_thres).reshape(-1,)
            # bestsnps_ind = np.argwhere(self._pvals > 0.8).reshape(-1,)[:3] # used to select worst pvals
            self.logger.debug("Selected {:d} pvalues below {:f} threshold".format(len(bestsnps_ind), pval_thres))
            # self.logger.debug("Selecting best {:d} pvalues out of {:d}".format(n, len(self.pvals)))
            # bestsnps_ind = np.argpartition(self.pvals, n)[:n]
        else:
            self.logger.debug("Selecting best {:d} Qscores out of {:d}".format(n, len(self.qscores)))
            bestsnps_ind = np.argpartition(self.qscores, -n)[-n:]
        # return bestsnps_ind
        return sorted(bestsnps_ind)

