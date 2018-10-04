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

class RevReg:


    def __init__(self, x, y, sigbeta2, comm, rank, ncore, null = 'perm', maf = None, cismasks = [], snps_cismasks = []):
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
        self._cismasks = cismasks
        self._snps_cismasks = snps_cismasks

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

    def slavejob_cismasks(self, cis_gt, cis_gx, sb2, sx2, maf):
        slv_gt  = cis_gt
        slv_gx  = cis_gx
        slv_sb2 = sb2
        slv_sx2 = sx2
        # print("Rank {:d}: Reporting {:d} SNPs, {:d} samples, {:d} genes".format(self.rank, slv_gt.shape[0], slv_gx.shape[1], slv_gx.shape[0]))
        if self.null == 'perm':
            p, q, mu, sig = crrstat.perm_null(slv_gt, slv_gx, slv_sb2, slv_sx2)
        elif self.null == 'maf':
            slv_maf = maf
            p, q, mu, sig = crrstat.maf_null (slv_gt, slv_gx, slv_sb2, slv_sx2, slv_maf)
        #self.logger.debug("Reporting from node {:d}. Sigma = ".format(self.rank) + np.array2string(sig) + "\n" )
        return p, q, mu, sig

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
        geno = self.comm.bcast(self.gt)
        expr = self.comm.bcast(expr, root = 0)
        sb2  = self.comm.bcast(sb2,  root = 0)
        sx2  = self.comm.bcast(sx2,  root = 0)
        maf  = self.comm.bcast(maf,  root = 0)
        self.comm.barrier()

        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        all_pvals   = np.array([])
        all_qscores = np.array([])
        all_mu      = np.array([])
        all_sigma   = np.array([])
        for i in range(len(slave_cismasks)):
            current_mask = slave_cismasks[i]
            current_snp_mask = slave_snps_cismasks[i]
            # print("Rank {:d}: Computing for cismask {:d}/{:d}".format(self.rank, i, len(slave_cismasks)-1))
            # print("Rank {:d}: Current mask has {:d} cis-genes removed".format(self.rank, len(current_mask)))
            gtcent_crop   = geno[current_snp_mask, :]
            sb2_crop = sb2[current_snp_mask]
            sx2_crop = sx2[current_snp_mask]
            if len(current_mask) > 0:
                expr_crop = np.delete(expr, current_mask, axis=0)
            else:
                expr_crop = expr
            pvals, qscores, mu, sigma = self.slavejob_cismasks(gtcent_crop, expr_crop, sb2_crop, sx2, maf)
            # print("Rank {:d}: Computed {:d} pvals".format(self.rank, len(pvals)))
            all_pvals = np.append(all_pvals, pvals)
            all_qscores = np.append(all_qscores, qscores)
            all_mu = np.append(all_mu, mu)
            all_sigma = np.append(all_sigma, sigma)
            # print("Rank {:d}: Computed so far {:d} pvals".format(self.rank, len(all_pvals)))
            # print("Rank {:d}: Computed so far {:d} qscores".format(self.rank, len(all_qscores)))
        
        apvals   = self.comm.gather(all_pvals,   root = 0)
        aqscores = self.comm.gather(all_qscores, root = 0)
        amu      = self.comm.gather(all_mu,      root = 0)
        asigma   = self.comm.gather(all_sigma,   root = 0)
        # atag     = self.comm.gather(slave_tags,  root = 0)

        # print(atag)
        if self.rank == 0:
            self._pvals   = np.concatenate(apvals)
            self._qscores = np.concatenate(aqscores)
            self._mu      = np.concatenate(amu)
            self._sigma   = np.concatenate(asigma)
            # print(np.concatenate(atag))
            print("Rank {:d}: all nodes computed a total of {:d} pvalues".format(self.rank, len(self._pvals)))
        else:
            assert aqscores is None
            assert apvals   is None
            assert amu      is None
            assert asigma   is None

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
                pvals, qscores, mu, sigma = self.slavejob_cismasks(gtcent_crop, expr_crop, sigbeta2_crop, self.sigx2, self.maf)
                self._pvals = pvals
                self._qscores = qscores
                self._mu = mu
                self._sigma = sigma
        return
