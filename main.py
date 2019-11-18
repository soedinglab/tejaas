#/usr/bin/env python

import numpy as np
import logging
import time
import itertools
import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from qstats.jpa import JPA
from qstats.revreg import RevReg
from qstats.jpa_pvals import JPAPVALS

from iotools.data import Data
from iotools.outhandler import Outhandler
from iotools import readmaf
from iotools import readqnull

# ==================================================
# Start MPI calculation
# =================================================


MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

if rank == 0: start_time = time.time()

args = Args(comm, rank)
logger = MyLogger(__name__)

gtcent = None
gtnorm = None
expr   = None
maf = None
snpinfo = None
masklist = None
maskcomp = None
if args.onlyjpa: qnull = None

if rank == 0:
    logger.debug("Using {:d} cores".format(ncore))
    data = Data(args)
    if args.simulate:
        data.simulate()
    else:
        data.load()
    gtcent = data.geno_centered
    gtnorm = data.geno_normed
    snpinfo = data.snpinfo
    expr = data.expression
    geneinfo = data.geneinfo
    masklist = data.cismasks_list
    maskcomp = data.cismasks_comp
    logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
    
    maf = readmaf.load(snpinfo, args.nullmodel, args.maf_file)
    if args.onlyjpa: qnull = readqnull.load(args.qnullfile)

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
snpinfo = comm.bcast(snpinfo, root = 0)
maf  = comm.bcast(maf, root = 0)
masklist = comm.bcast(masklist, root = 0)
maskcomp = comm.bcast(maskcomp, root = 0)
if args.onlyjpa: qnull = comm.bcast(qnull, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()

if rank == 0: logger.debug("Computing JPA")

if args.onlyjpa:
    jpa = JPAPVALS(gtnorm, expr, qnull, comm, rank, ncore, args.onlyjpa, masklist)
    jpa.compute()
else:
    jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist)
    jpa.compute()

##### NULL for JPA
# if rank == 0: jpanull_time = time.time()
# rand_gtnorm = gtnorm.copy()
# jpanull = None
# for i in range(333):
#     np.random.shuffle(rand_gtnorm.T)
#     jpa_null = JPA(rand_gtnorm, expr, comm, rank, ncore, args.jpa, masklist)
#     jpa_null.compute()
#     if rank == 0:
#         if jpanull is not None:
#             jpanull = np.vstack((jpanull, jpa_null.scores))
#         else:
#             jpanull = jpa_null.scores
# if rank == 0: jpanull_endtime = time.time()
# if rank == 0: logger.info("JPA Null calculation time: {:g} seconds".format (jpanull_endtime - jpanull_time))
# if rank == 0: print(jpanull.shape) # times x qscores

# if rank == 0: np.savetxt("/cbscratch/franco/tejaas_output/tests/jpanull_scores.txt", jpanull) # for testing purposes
# 2019-03-05 21:08:24,372 | __main__ | INFO | JPA Null calculation time: 8886.86 seconds
##### End of NULL for JPA

if rank == 0: jpa_time = time.time()

# # Select the SNPs with JPA score above threshold for RevReg
if args.jpa and args.rr:
    select = None
    if rank == 0:
        select = np.where(jpa.scores > args.jpacut)[0]
        logger.debug('JPA threshold: {:g}, and number of SNPs retained: {:d}'.format(args.jpacut, select.shape[0]))
    qselect = comm.bcast(select, root = 0)
    comm.barrier()
    #geno = geno[qselect, :]
    #broadcast this new genotype

if args.rr:

    if args.dynamic:
        logger.debug("-- DINAMICALLY ADJUSTING SIGMA BETA -- ")
        # adjust sigma_beta2 for each SNP (aprox)
        Yt = expr.T
        U, S, Vt  = np.linalg.svd(Yt, full_matrices=False)
        sigmax2   = np.var(gtnorm, axis = 1)
        S2_median = np.median(np.square(S))
        sigbeta2  = sigmax2 / S2_median
    else:
        sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])

    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)

    if args.pms:
        rr.compute(get_betas = True)
        if rank == 0:
            ohandle = Outhandler(args, snpinfo, geneinfo)
            ohandle.write_rr_out(jpa, rr, prefix = "_it0")
    else:
        rr.compute(get_betas = False)

if rank == 0: rr_time = time.time()

if args.pms:
    ###### New implementation of sparsity and null model calculation
    best_snp_indices = None
    gene_indices     = None
    if rank == 0:
        best_snp_indices = rr.select_best_SNPs(pval_thres = 1e-6, use_pvals = True)
        gene_indices     = rr.select_best_genes(rr.betas[best_snp_indices,:], n=100)
    gene_indices     = comm.bcast(gene_indices, root = 0)
    best_snp_indices = comm.bcast(best_snp_indices, root = 0)
    comm.barrier()

    if len(best_snp_indices):
        logger.debug("Prunning {:d} masks".format(len(maskcomp)))
        newmasks = rr.prune_masks(maskcomp, list(best_snp_indices))
        logger.debug("Got {:d} pruned masks".format(len(newmasks)))
        
        stime = time.time()
        qnull = None
        if args.nullmodel == 'perm':
            gt   = gtcent[best_snp_indices,:]

        sb2        = sigbeta2[best_snp_indices]
        sb2_sparse = np.repeat(0.05 ** 2, len(best_snp_indices))
        for i in range(1000):
            # THIS IS WRONG: each core shuffles whatever
            np.random.shuffle(gt.T)
            if rank == 0: logger.debug(str(rank) + ":---> computing rr for null")
            rr_null = RevReg(gt, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = newmasks)
            rr_null.compute(get_betas = True)

            null_gene_indices = None
            if rank == 0:
                null_gene_indices = rr_null.select_best_genes(rr_null.betas, n=1000)
            null_gene_indices = comm.bcast(null_gene_indices, root = 0)
            comm.barrier()

            rr_null.sb2 = sb2_sparse
            if rank == 0: logger.debug("---> computing rr for sparse null")
            rr_null.compute_sparse(null_gene_indices, get_betas = False)

            if rank == 0:
                if qnull is not None:
                    qnull = np.vstack((qnull, rr_null.scores))
                else:
                    qnull = rr_null.scores
            qnull = comm.bcast(qnull, root = 0)

        #np.savetxt("/cbscratch/franco/tejaas_output/tests/qnull_scores.txt", qnull) # for testing purposes
        logger.debug("TIME-->MPIcompute_null took {:g} seconds".format(time.time() - stime))
        
        logger.debug("---> computing rr for sparse")
        rr_sparse = RevReg(gtcent[best_snp_indices,:], expr, sb2_sparse, comm, rank, ncore, null = args.nullmodel, maf = maf[best_snp_indices], masks = newmasks)
        rr_sparse.compute_sparse(gene_indices, qnull, get_betas = False)
        ##### end

if rank == 0: pms_time = time.time()
 
# Output handling only from master node // move it to module
if rank == 0: 
    ohandle = Outhandler(args, snpinfo, geneinfo)
    if args.onlyjpa:
        ohandle.write_jpa_pvals(jpa)
    if args.rr:
        ohandle.write_rr_out(jpa, rr, write_betas = True)
    if args.pms:
        ohandle.write_rr_out(jpa, rr_sparse, selected_snps = best_snp_indices, selected_genes = gene_indices, prefix = "_it1")

if rank == 0: write_time = time.time()

if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("JPA calculation time: {:g} seconds".format (jpa_time - read_time))
    logger.info("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    logger.info("PMS calculation time: {:g} seconds".format(pms_time - rr_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - rr_time))
    logger.info("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
