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
import copy
import pdb

class RevReg:


    def __init__(self, x, y, sigbeta2, comm, rank, ncore, null = 'perm', maf = None, cismasks = [], snps_cismasks = [], outdir = "out"):
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
        self._betas = None
        self._cismasks = cismasks
        self._snps_cismasks = snps_cismasks
        self._selected_snps = []
        self._selected_genes = None
        self.outdir = outdir

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

    @property
    def betas(self):
        return self._betas

    @property
    def selected_snps(self):
        return self._selected_snps

    @property
    def selected_genes(self):
        return self._selected_genes

    def write_rr_out(self, suffix, snpinfo, geneinfo):
        mysnpinfo = snpinfo
        if len(self.selected_snps):
            mysnpinfo = [snpinfo[int(i)] for i in self.selected_snps]
        fname = self.outdir + "_rr_" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
            print("write rr out: ", len(mysnpinfo), len(self.pvals))
            for i, x in enumerate(mysnpinfo):
                f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, self.scores[i], self.null_mu[i], self.null_sigma[i], self.pvals[i]))
        np.savetxt(self.outdir + "_betas_" + suffix + ".txt", self.betas, fmt='%1.4e')
        if self.selected_genes is not None:
            np.savetxt(self.outdir + "_selected_genes_" + suffix + ".txt", self.selected_genes, fmt='%i')

    def basejob(self, cis_gt, cis_gx, sb2, sx2, maf, getb=False):
        slv_gt  = cis_gt
        slv_gx  = cis_gx
        slv_sb2 = sb2
        slv_sx2 = sx2
        nsnps   = slv_gt.shape[0]
        ngenes  = slv_gx.shape[0]
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
            if getb:
                B = crrstat.crrbetas(slv_gt, slv_gx, slv_sb2)
            else:
                B = np.zeros(nsnps * ngenes)    
            B = B.reshape(nsnps, ngenes)
        elif self.null == 'maf':
            slv_maf = maf
            p, q, mu, sig = crrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
            B = np.zeros(nsnps * ngenes)
            B = B.reshape(nsnps, ngenes)
        return p, q, mu, sig, B


    def slavejob_cismasks_wrapper(self, gt, gx, sb2, sx2, maf, snp_masks, gene_masks, getb=False):
        local_p   = np.array([])
        local_q   = np.array([])
        local_mu  = np.array([])
        local_s   = np.array([])
        local_b   = np.array([])

        self.logger.debug("Rank {:d}: wrapper: I have {:d} snps, geno is {:s} and gx is {:s} .".format(self.rank, gt.shape[0], str(gt.shape), str(gx.shape)))
        for i in range(len(gene_masks)): #range(3):
            current_mask = gene_masks[i]
            current_snp_mask = snp_masks[i]
            gt_crop  = gt[current_snp_mask, :]
            sb2_crop = sb2[current_snp_mask]
            sx2_crop = sx2[current_snp_mask]
            if len(current_mask) > 0:
                gx_crop = np.delete(gx, current_mask, axis=0)
            else:
                gx_crop = gx
            p, q, mu, s, b = self.basejob(gt_crop, gx_crop, sb2_crop, sx2, maf, getb=True)
            self.logger.debug("Rank {:d}: Computed {:d} pvals and {:s} betas".format(self.rank, len(p), str(b.shape)))

            # make array of betas compatible with previous expr data matrix (NxG), fill in with zeros the masked genes
            padBeta = np.zeros( (len(current_snp_mask), gx.shape[0]) )
            inv_ind = np.delete(np.arange(gx.shape[0]), current_mask)
            padBeta[:, inv_ind] = b

            local_p  = np.append(local_p, p)
            local_q  = np.append(local_q, q)
            local_mu = np.append(local_mu, mu)
            local_s  = np.append(local_s, s)
            if local_b.shape[0] == 0:
                local_b = padBeta
            else:
                local_b = np.vstack((local_b, padBeta))
        return local_p, local_q, local_mu, local_s, local_b

    def mpicompute(self):
        if self.rank == 0:
            geno = self.gt
            expr = self.gx
            sb2  = self.sigbeta2
            sx2  = self.sigx2
            maf  = self.maf

            # split masks in 4 sublists
            masks_list = mpihelper.split_1darray(self._cismasks, self.ncore)
            snp_list = mpihelper.split_1darray(self._snps_cismasks, self.ncore)
        else:
            geno = None
            expr = None
            sb2  = None
            sx2  = None
            maf  = None
            masks_list = None
            snp_list = None

        slave_cismasks      = self.comm.scatter(masks_list, root = 0)
        slave_snps_cismasks = self.comm.scatter(snp_list, root = 0)
        geno = self.comm.bcast(geno, root = 0)
        expr = self.comm.bcast(expr, root = 0)
        sb2  = self.comm.bcast(sb2,  root = 0)
        sx2  = self.comm.bcast(sx2,  root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        self.comm.barrier()

        # ==================================
        # Data sent. Now do the calculations
        # ==================================

        self.logger.debug("Rank {:d}: calculating over {:d} gene masks and {:d} slave_snps_cismasks".format(self.rank, len(slave_cismasks), len(slave_snps_cismasks)))

        local_p, local_q, local_mu, local_s, local_b = self.slavejob_cismasks_wrapper(geno, expr, sb2, sx2, maf, slave_snps_cismasks, slave_cismasks, getb=True)

        self.logger.debug("Rank {:d}: Done. Obtained {:s} betas".format(self.rank, str(local_b.shape)))

        ## Synchronize all slaves and collect results
        self.comm.barrier()
        apvals   = self.comm.gather(local_p,   root = 0)
        aqscores = self.comm.gather(local_q, root = 0)
        amu      = self.comm.gather(local_mu,      root = 0)
        asigma   = self.comm.gather(local_s,   root = 0)


        ## Collect beta values this way to prevent overflow (arrays are to big)
        # received counts collects the number of betas calculated by each core (snps * ngenes)
        flatbeta = local_b.reshape(-1,).copy()
        flatlength = len(flatbeta)
        self.comm.barrier()
        received_counts = self.comm.gather(flatlength)
        
        if self.rank == 0:
            recvbuf = np.zeros(np.sum(received_counts), dtype=np.float64)
        else:
            recvbuf = None  
        
        self.comm.Gatherv(sendbuf=flatbeta, recvbuf=(recvbuf, received_counts), root = 0)
        
        #########

        if self.rank == 0:
            self._pvals   = np.concatenate(apvals)
            self._qscores = np.concatenate(aqscores)
            self._mu      = np.concatenate(amu)
            self._sigma   = np.concatenate(asigma)
            self._betas   = recvbuf.reshape(self.gt.shape[0], self.gx.shape[0])
            self.logger.debug("Rank {:d}: all nodes computed a total of {:d} pvalues and {:s} betas".format(self.rank, len(self._pvals), str(self._betas.shape)))
        else:
            assert aqscores is None
            assert apvals   is None
            assert amu      is None
            assert asigma   is None
            newarr = None
        return True

    def slavejob_sparse(self, slave_geno, expr, sb2, slave_sx2, maf, slave_betas, nbetas, selected_genes = None):
        local_p  = np.array([])
        local_q  = np.array([])
        local_mu = np.array([])
        local_s  = np.array([])
        local_b  = np.array([])
        local_indices = np.array([], dtype=int)
        print("Rank {:d}: expr: {:s} ".format(self.rank, str(expr.shape)))
        print("Rank {:d}: slave_geno: {:s}".format(self.rank, str(slave_geno.shape)))
        print("Rank {:d}: slave_betas: {:s}".format(self.rank, str(slave_betas.shape)))
        for j in range(slave_betas.shape[0]):  # iterate over all snps
            beta_j        = np.abs(slave_betas[j])
            bestbetas_ind = np.argpartition(beta_j, -nbetas)[-nbetas:]
            print("bestbetas_ind: ", bestbetas_ind.shape)
            if selected_genes is not None:
                print("Rank {:d}: selected_genes size: ".format(self.rank), selected_genes.shape)
                ix = selected_genes[j,:][bestbetas_ind]
                print("ix:", ix.shape, ix)
                ix = ix.reshape(-1,)
                print("ix reshape:", ix )
                local_indices = np.append(local_indices, ix)
                local_expr    = expr[ix,:]
            else:
                print("Rank {:d}: selected_genes: ".format(self.rank), selected_genes)
                local_indices = np.append(local_indices, bestbetas_ind)
                local_expr    = expr[bestbetas_ind,:]                
            # snpgt         = np.ascontiguousarray(slave_geno[j,:][np.newaxis]) # because is only one row, we need to make it properly 1 x nsample -> ascontiguous didn't work for me, more testint needed
            snpgt         = slave_geno[j,:][np.newaxis].copy()
            p, q, mu, s, b = self.basejob(snpgt, local_expr, sb2[j].reshape(-1,), slave_sx2[j].reshape(-1,), maf, getb=True)
            local_p   = np.append(local_p, p)
            local_q = np.append(local_q, q)
            local_mu      = np.append(local_mu, mu)
            local_s   = np.append(local_s, s)
            local_b   = np.append(local_b, b)
        self.logger.debug("Rank {:d}: Sparse: Computed {:d} pvals and {:s} betas".format(self.rank, len(local_p), str(local_b.shape)))  
        return local_p, local_q, local_mu, local_s, local_b, local_indices

    def compute_sparse(self, ncutoff = None, nbetas = 1000):
        if self.rank == 0:
            if self.pvals is None:
                self.logger.error("No pvalues found, first run compute()")
                return False
            
            if ncutoff is not None and ncutoff > 0:
                ncutoff = min(ncutoff, self.gt.shape[0])
                self.logger.debug("Selecting best {:d}/{:d} Q scores".format(ncutoff, len(self.scores)))
                bestsnps_ind = np.argpartition(self._qscores, -ncutoff)[-ncutoff:]
            else:
                self.logger.error("Specify N cutoff")
                return False
            # # Implementation for SNP p-value cutoff
            # elif pval_cutoff is not None and pval_cutoff > 0:
            #     bestsnps_ind = np.argwhere(self._pvals < pval_cutoff).reshape(-1,)
            #     if len(bestsnps_ind) == 0:
            #         self.logger.error("No SNP meets the pval cutoff criteria")
            #         return False
            #     if len(bestsnps_ind) == len(self._pvals):
            #         self.logger.error("Number of selected SNPs is equal to current number. No point in running sparse computation")
            #         return False
            #     self.logger.debug("Selected {:d} SNPs with pval < {:d} ".format(len(bestsnps_ind), pval_cutoff))

            if len(self._selected_snps):
                self._selected_snps = self._selected_snps[bestsnps_ind]
            else:
                self._selected_snps = bestsnps_ind

            if nbetas >= self._betas.shape[1]:
                self.logger.error("Number of best betas to select is to large (max {:d}). No point in running sparsity".format(self._betas.shape[1]))
                return False

            expr     = self.gx
            sb2      = self.sigbeta2
            maf      = self.maf
            split_geno  = mpihelper.split_genotype(self.gt[bestsnps_ind,:], self.ncore)
            split_betas = mpihelper.split_genotype(self._betas[bestsnps_ind,:], self.ncore)
            split_sx2   = mpihelper.split_1darray(self.sigx2[bestsnps_ind], self.ncore)
            snp_per_node = [x.shape[0] for x in split_geno]
            if self._selected_genes is not None:
                print("I have {:s} selected genes".format(str(self._selected_genes.shape)))
                split_selected_genes = mpihelper.split_genotype(self._selected_genes[bestsnps_ind,:], self.ncore)
            else:
                split_selected_genes = [None] * self.ncore
            self.logger.debug("Send SNPs to each node: {:s}".format(str(snp_per_node)))
        else:
            expr         = None
            sb2          = None
            maf          = None
            split_geno   = None
            split_betas  = None
            split_sx2    = None
            nmax         = None
            snp_per_node = None
            split_selected_genes = None            

        expr        = self.comm.bcast(expr, root = 0)  
        sb2         = self.comm.bcast(sb2, root = 0)
        maf         = self.comm.bcast(maf,  root = 0)
        slave_geno  = self.comm.scatter(split_geno, root = 0)
        slave_betas = self.comm.scatter(split_betas, root = 0)
        slave_sx2   = self.comm.scatter(split_sx2, root = 0)
        nmax        = self.comm.scatter(snp_per_node, root = 0)
        slave_selected_genes = self.comm.scatter(split_selected_genes, root = 0)
        self.comm.barrier()

        ### Start sparse iteration of RR for each SNP ###
        self.logger.debug("Rank {:d}: Sparse: Calculating RR iteration over best {:d} genes (biggest beta values)".format(self.rank, nbetas))
        self.logger.debug("Rank {:d}: Sparse: I have {:s} betas, {:d} snps, geno is {:s} and {:s} sx2.".format(self.rank, str(slave_betas.shape), nmax, str(slave_geno.shape), str(slave_sx2.shape)))
 
        local_p, local_q, local_mu, local_s, local_b, local_indices = self.slavejob_sparse(slave_geno, expr, sb2, slave_sx2, maf, slave_betas, nbetas, slave_selected_genes)  

        self.comm.barrier()
        apvals   = self.comm.gather(local_p, root = 0)
        aqscores = self.comm.gather(local_q, root = 0)
        amu      = self.comm.gather(local_mu, root = 0)
        asigma   = self.comm.gather(local_s, root = 0)
        
        # gather betas and indices
        flatbeta = local_b.reshape(-1,).copy()
        flatindices = local_indices.reshape(-1,).copy()
        flatlengths = len(flatbeta)
        self.logger.debug("Rank {:d}: Sparse: I also have {:s} -> {:d} indices".format(self.rank, str(local_indices.shape), len(flatindices)))

        self.comm.barrier()
        received_counts = np.array(self.comm.gather(flatlengths))
        print("received counts: ", received_counts)

        if self.rank == 0:
            recvbuf_betas = np.zeros(np.sum(received_counts), dtype=np.float64)
            recvbuf_indices = np.empty(np.sum(received_counts), dtype=int)
        else:
            recvbuf_betas = None
            recvbuf_indices = None
        self.comm.Gatherv(sendbuf=flatbeta, recvbuf=(recvbuf_betas, received_counts), root = 0)
        self.comm.Gatherv(sendbuf=flatindices, recvbuf=(recvbuf_indices, received_counts), root = 0)
        
        if self.rank == 0:
            self.logger.debug("Rank {:d}: Sparse: Overwriting {:d} p-values with new {:d} pvals.".format(self.rank, len(self._pvals), len(np.concatenate(apvals))))
            self._pvals   = np.concatenate(apvals)
            self._qscores = np.concatenate(aqscores)
            self._mu      = np.concatenate(amu)
            self._sigma   = np.concatenate(asigma)
            self._betas   = recvbuf_betas.reshape(len(self._selected_snps), nbetas)
            self._selected_genes = recvbuf_indices.reshape(len(self._selected_snps), nbetas)

            self.logger.debug("Rank {:d}: Sparse: all nodes computed a total of {:d} pvalues, {:s} betas, and {:s} selected genes".format(self.rank, len(self._pvals), str(self._betas.shape), str(self._selected_genes.shape)))
        else:
            assert aqscores is None
            assert apvals   is None
            assert amu      is None
            assert asigma   is None
        return

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            for i in range(len(self._cismasks)):
                # print("Computing for cismask {:d}/{:d}".format(i, len(self._cismasks)))
                # print("Current mask has {:d} cis-genes removed".format(len(self._cismasks[i])))
                gtcent_crop   = self.gt[self._snps_cismasks[i], :]
                sigbeta2_crop = self.sigbeta2[self._snps_cismasks[i]]
                if len(self._cismasks[i]) > 0:
                    expr_crop = np.delete(self.gx, self._cismasks[i], axis=0)
                else:
                    expr_crop = expr
                pvals, qscores, mu, sigma, betas = self.basejob(gtcent_crop, expr_crop, sigbeta2_crop, self.sigx2, self.maf)
                self._pvals = pvals
                self._qscores = qscores
                self._mu = mu
                self._sigma = sigma
        return
