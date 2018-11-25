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

class RevReg:


    def __init__(self, x, y, sigbeta2, comm, rank, ncore, null = 'perm', maf = None, cismasks = [], snps_cismasks = [], outdir = "out", snpinfo = None):
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
        self.outdir = outdir
        self.snpinfo = snpinfo

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
    def betas(self):
        return self._betas

    @property
    def selected_snps(self):
        return self._selected_snps

    @property
    def null_sigma(self):
        return self._sigma

    def write_rr_out(self, suffix):
        mysnpinfo = self.snpinfo
        if len(self.selected_snps):
            mysnpinfo = [self.snpinfo[int(i)] for i in self.selected_snps]
        fname = self.outdir + "_rr_" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
            print("write rr out: ", len(mysnpinfo), len(self.pvals))
            for i, x in enumerate(mysnpinfo):
                print(x.varid, x.bp_pos, self.scores[i], self.null_mu[i], self.null_sigma[i], self.pvals[i])
                f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, self.scores[i], self.null_mu[i], self.null_sigma[i], self.pvals[i]))

    def slavejob_cismasks(self, cis_gt, cis_gx, sb2, sx2, maf, getb=False):
        slv_gt  = cis_gt
        slv_gx  = cis_gx
        slv_sb2 = sb2
        slv_sx2 = sx2
        nsnps   = slv_gt.shape[0]
        ngenes  = slv_gx.shape[0]
        # print("Rank {:d}: slavejob_cismasks: Reporting {:d} SNPs, {:d} samples, {:d} genes".format(self.rank, slv_gt.shape[0], slv_gx.shape[1], slv_gx.shape[0]))
        # print("Rank {:d}: slavejob_cismasks: slv_gt is {:s} and slv_gx is {:s}".format(self.rank, str(slv_gt.shape), str(slv_gx.shape)))
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
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig, B


    def slavejob_cismasks_wrapper(self, gt, gx, sb2, sx2, maf, snp_masks, gene_masks, getb=False):
        local_pvals   = np.array([])
        local_qscores = np.array([])
        local_mu      = np.array([])
        local_sigma   = np.array([])
        local_betas   = np.array([])

        print("Rank {:d}: wrapper: I have {:d} snps, geno is {:s} and gx is {:s} .".format(self.rank, gt.shape[0], str(gt.shape), str(gx.shape)))
        print("Rank {:d}: wrapper: calculating over {:d} gene masks".format(self.rank, len(gene_masks)))
        # self.comm.barrier()
        for i in range(len(gene_masks)): #range(3):
            current_mask = gene_masks[i]
            current_snp_mask = snp_masks[i]
            # print("Rank {:d}: Computing for cismask {:d}/{:d}".format(self.rank, i, len(gene_masks)-1))
            # print("Rank {:d}: Current mask has {:d} cis-genes removed affecting {:d} SNPs".format(self.rank, len(current_mask), len(current_snp_mask)))
            # self.comm.barrier()
            gt_crop  = gt[current_snp_mask, :]
            sb2_crop = sb2[current_snp_mask]
            sx2_crop = sx2[current_snp_mask]
            if len(current_mask) > 0:
                gx_crop = np.delete(gx, current_mask, axis=0)
            else:
                gx_crop = gx
            print("Rank {:d}: wrapper: gt_crop is {:s} and gx_crop is {:s}".format(self.rank, str(gt_crop.shape), str(gx_crop.shape)))
            pvals, qscores, mu, sigma, betas = self.slavejob_cismasks(gt_crop, gx_crop, sb2_crop, sx2, maf, getb=True)
            # self.comm.barrier()
            print("Rank {:d}: Computed {:d} pvals".format(self.rank, len(pvals)))
            print("Rank {:d}: Computed {:s} betas".format(self.rank, str(betas.shape)))

            # make beta array size compatible with previous expr data, fill in with zeros
            padBeta = np.zeros( (len(current_snp_mask), gx.shape[0]) )
            inv_ind = np.delete(np.arange(gx.shape[0]), current_mask)
            padBeta[:, inv_ind] = betas

            local_pvals = np.append(local_pvals, pvals)
            local_qscores = np.append(local_qscores, qscores)
            local_mu = np.append(local_mu, mu)
            local_sigma = np.append(local_sigma, sigma)
            if local_betas.shape[0] == 0:
                local_betas = padBeta
            else:
                local_betas = np.vstack((local_betas, padBeta))
        return local_pvals, local_qscores, local_mu, local_sigma, local_betas

    def mpicompute_cismasks(self):
        if self.rank == 0:
            geno = self.gt
            expr = self.gx
            sb2  = self.sigbeta2
            sx2  = self.sigx2
            maf  = self.maf

            # split masks in 4 sublists
            nmasks = len(self._cismasks)
            maxmasks = int(nmasks / self.ncore)
            offset = 0
            masks_list = [None for i in range(self.ncore)]
            snp_list = [None for i in range(self.ncore)]
            # tags = [None for i in range(self.ncore)]
            for dest in range(self.ncore-1):
                start = offset
                end = offset + maxmasks
                masks_list[dest] = self._cismasks[start:end]
                snp_list[dest]   = self._snps_cismasks[start:end]
                # tags[dest] = np.arange(1,100) * (dest+1)
                offset += maxmasks
            masks_list[self.ncore-1] = self._cismasks[offset:]
            snp_list[self.ncore-1]   = self._snps_cismasks[offset:]
            # tags[self.ncore-1] = np.arange(1,100) * (self.ncore)
        else:
            geno = None
            expr = None
            sb2  = None
            sx2  = None
            maf  = None
            masks_list = None
            snp_list = None
            # tags = None

        slave_cismasks      = self.comm.scatter(masks_list, root = 0)
        slave_snps_cismasks = self.comm.scatter(snp_list, root = 0)
        # slave_tags          = self.comm.scatter(tags, root = 0)
        geno = self.comm.bcast(geno, root = 0)
        expr = self.comm.bcast(expr, root = 0)
        sb2  = self.comm.bcast(sb2,  root = 0)
        sx2  = self.comm.bcast(sx2,  root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        self.comm.barrier()

        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        # local_pvals   = np.array([])
        # local_qscores = np.array([])
        # local_mu      = np.array([])
        # local_sigma   = np.array([])
        # local_betas   = np.array([])

        print("Rank {:d}: calculating over {:d} gene masks".format(self.rank, len(slave_cismasks)))
        print("Rank {:d}: slave_snps_cismasks {:d}".format(self.rank, len(slave_snps_cismasks)))
        flatened_list = [item for sublist in slave_snps_cismasks for item in sublist]
        print("Rank {:d}: flatened slave_snps_cismasks {:d}".format(self.rank, len(flatened_list)))
        self.comm.barrier()
        local_pvals, local_qscores, local_mu, local_sigma, local_betas = self.slavejob_cismasks_wrapper(geno, expr, sb2, sx2, maf, slave_snps_cismasks, slave_cismasks, getb=True)
        
        # for i in range(len(slave_cismasks)): #range(3):
        #     current_mask = slave_cismasks[i]
        #     current_snp_mask = slave_snps_cismasks[i]
        #     # print("Rank {:d}: Computing for cismask {:d}/{:d}".format(self.rank, i, len(slave_cismasks)-1))
        #     # print("Rank {:d}: Current mask has {:d} cis-genes removed affecting {:d} SNPs".format(self.rank, len(current_mask), len(current_snp_mask)))
        #     gtcent_crop   = geno[current_snp_mask, :]
        #     sb2_crop = sb2[current_snp_mask]
        #     sx2_crop = sx2[current_snp_mask]
        #     if len(current_mask) > 0:
        #         expr_crop = np.delete(expr, current_mask, axis=0)
        #     else:
        #         expr_crop = expr
        #     pvals, qscores, mu, sigma, betas = self.slavejob_cismasks(gtcent_crop, expr_crop, sb2_crop, sx2, maf, getb=True)
        #     print("Rank {:d}: Computed {:d} pvals".format(self.rank, len(pvals)))
        #     print("Rank {:d}: Computed {:s} betas".format(self.rank, str(betas.shape)))

        #     # make beta array size compatible with previous expr data, fill in with zeros
        #     padBeta = np.zeros( (len(current_snp_mask), expr.shape[0]) )
        #     inv_ind = np.delete(np.arange(expr.shape[0]), current_mask)
        #     padBeta[:, inv_ind] = betas

        #     local_pvals = np.append(local_pvals, pvals)
        #     local_qscores = np.append(local_qscores, qscores)
        #     local_mu = np.append(local_mu, mu)
        #     local_sigma = np.append(local_sigma, sigma)
        #     if local_betas.shape[0] == 0:
        #         local_betas = padBeta
        #     else:
        #         local_betas = np.vstack((local_betas, padBeta))
        # raise

        print("Rank {:d}: Finally {:s} betas obtained".format(self.rank, str(local_betas.shape)))
        print("Rank {:d}: Computed so far {:d} pvals".format(self.rank, len(local_pvals)))
        print("Rank {:d}: Computed so far {:d} qscores".format(self.rank, len(local_qscores)))

        flatbeta = local_betas.reshape(-1,).copy()

        ## Synchronize all slaves and collect results
        self.comm.barrier()
        apvals   = self.comm.gather(local_pvals,   root = 0)
        aqscores = self.comm.gather(local_qscores, root = 0)
        amu      = self.comm.gather(local_mu,      root = 0)
        asigma   = self.comm.gather(local_sigma,   root = 0)


        ## FIXED: Sometime it fails to gather flatbeta due to overflow (arrays are to big)
        recvbuf = None
        if self.rank == 0:
            buffer_sizes = []
            for i in range(self.ncore):
                flatlist = []
                for corelist in snp_list[i]:
                    for sublist in corelist:
                        flatlist.append(sublist)
                buffer_sizes.append(len(flatlist))
            print("Buffer sizes: ", buffer_sizes)
            print("Number of snp-betas to receive from each slave: ", buffer_sizes)
            ngenes = self.gx.shape[0]
            flat_sizes = np.array([rows * ngenes for rows in buffer_sizes]) # rows = nsnps per node (accounting for all snp_masks)
            total_rows = sum(np.array(buffer_sizes))
            recvbuf = np.zeros(sum(flat_sizes), dtype=np.float64)
        else:
            flat_sizes = None        
        
        # recvbuf would be equivalent to abetas
        self.comm.Gatherv(sendbuf=flatbeta, recvbuf=(recvbuf, flat_sizes), root = 0)
        
        #########

        # abetas   = self.comm.gather(flatbeta,   root = 0)
        # self.comm.barrier()
        
        if self.rank == 0:

            self._pvals   = np.concatenate(apvals)
            self._qscores = np.concatenate(aqscores)
            self._mu      = np.concatenate(amu)
            self._sigma   = np.concatenate(asigma)
            # tmp_betas     = np.concatenate(abetas)
            # self._betas   = tmp_betas.reshape(self.gt.shape[0], self.gx.shape[0])
            newarr        = recvbuf.reshape(total_rows, ngenes)
            self._betas   = newarr

            # write data for iteration 0
            self.write_rr_out("it0")


            print("Rank {:d}: all nodes computed a total of {:d} pvalues and {:s} betas".format(self.rank, len(self._pvals), str(self._betas.shape)))
            signif_snps_ind = np.argwhere(self._pvals < 10e-3).reshape(-1,)
            self._selected_snps = signif_snps_ind
            print("Rank {:d}: selected {:d} signifcative (<10e-3) pvalues ".format(self.rank, len(signif_snps_ind)))

            newgeno  = mpihelper.split_genotype(self.gt[signif_snps_ind,:], self.ncore)
            newbetas = mpihelper.split_genotype(self._betas[signif_snps_ind,:], self.ncore)
            newsx2   = mpihelper.split_1darray(self.sigx2[signif_snps_ind], self.ncore)
            snp_per_node = [x.shape[0] for x in newgeno]
        else:
            assert aqscores is None
            assert apvals   is None
            assert amu      is None
            assert asigma   is None
            # assert abetas   is None
            newgeno = None
            newbetas = None
            newsx2 = None
            snp_per_node = None
            nmax = None

        slave_geno  = self.comm.scatter(newgeno, root = 0)
        slave_betas = self.comm.scatter(newbetas, root = 0)
        slave_sx2   = self.comm.scatter(newsx2, root = 0)
        # already broadcasted before # expr = self.comm.bcast(expr, root = 0)
        nmax = self.comm.scatter(snp_per_node, root = 0)
        self.comm.barrier()

        ### Second iteration of RR ###
        Gs = [1000]
        for i, G in enumerate(Gs):
            print("Rank {:d}: Calculating RR iteration {:d}/{:d}, over best {:d} genes".format(self.rank, i+1, len(Gs), G))
            print("Rank {:d}: I have {:s} betas, {:d} snps, geno is {:s} and {:s} sx2.".format(self.rank, str(slave_betas.shape), nmax, str(slave_geno.shape), str(slave_sx2.shape)))
            local_Gpvals   = np.array([])
            local_Gqscores = np.array([])
            local_Gmu      = np.array([])
            local_Gsigma   = np.array([])
            for j in range(slave_betas.shape[0]):  # iterate over all snps basically
                beta_j = slave_betas[j]
                bestG_indices = np.argpartition(np.abs(beta_j), -G)[-G:]
                bestG_expr    = expr[bestG_indices,:]
                # snpgt         = np.ascontiguousarray(slave_geno[j,:][np.newaxis]) # because is only one row, we need to make it properly 1 x nsample
                snpgt         = slave_geno[j,:][np.newaxis].copy()
                # single_sb2    = sb2_crop[j] # this doesn't mind because all snps have same sb2
                Gpvals, Gqscores, Gmu, Gsigma, _tmp = self.slavejob_cismasks(snpgt, bestG_expr, sb2[j].reshape(-1,), slave_sx2[j].reshape(-1,), maf, getb=False)
                local_Gpvals   = np.append(local_Gpvals, Gpvals)
                local_Gqscores = np.append(local_Gqscores , Gqscores)
                local_Gmu      = np.append(local_Gmu , Gmu)
                local_Gsigma   = np.append(local_Gsigma , Gsigma)
            print("Rank {:d}: Computed {:d} pvals".format(self.rank, len(local_Gpvals)))
            local_pvals   = local_Gpvals
            local_qscores = local_Gqscores
            local_mu      = local_Gmu
            local_sigma   = local_Gsigma
        

        apvals   = self.comm.gather(local_Gpvals,   root = 0)
        aqscores = self.comm.gather(local_Gqscores, root = 0)
        amu      = self.comm.gather(local_Gmu,      root = 0)
        asigma   = self.comm.gather(local_Gsigma,   root = 0)
        # abetas   = self.comm.gather(flatbeta,   root = 0)
        # atag     = self.comm.gather(slave_tags,  root = 0)

        # print(atag)
        if self.rank == 0:
            print("Rank {:d}: Overwriting {:d} p-values with new {:d} pvals.".format(self.rank, len(self._pvals), len(np.concatenate(apvals))))
            self._pvals   = np.concatenate(apvals)
            print(self._pvals)
            self._qscores = np.concatenate(aqscores)
            self._mu      = np.concatenate(amu)
            self._sigma   = np.concatenate(asigma)
            # tmp_betas     = np.concatenate(abetas)
            # self._betas   = tmp_betas.reshape(self.gt.shape[0], self.gx.shape[0])
            print("Rank {:d}: all nodes computed a total of {:d} pvalues".format(self.rank, len(self._pvals)))

            self.write_rr_out("it1")
        else:
            assert aqscores is None
            assert apvals   is None
            assert amu      is None
            assert asigma   is None
            # assert abetas   is None

        return

    def compute(self):
        if self.mpi:
            self.mpicompute_cismasks()
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
                pvals, qscores, mu, sigma, betas = self.slavejob_cismasks(gtcent_crop, expr_crop, sigbeta2_crop, self.sigx2, self.maf)
                self._pvals = pvals
                self._qscores = qscores
                self._mu = mu
                self._sigma = sigma
        return
