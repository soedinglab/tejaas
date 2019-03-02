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
from iotools.data import Data
from iotools.outhandler import Outhandler
from iotools import readmaf

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

#if args.cismasking:
#    from qstats.revreg_cis import RevReg
#else:
#    from qstats.revreg import RevReg

gtcent = None
gtnorm = None
expr   = None
maf = None
snpinfo = None
masklist = None
maskcomp = None

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

    ### Franco shuffle on genes expr
    # if args.randomize_file is not None:
    #     rindices = np.loadtxt(args.randomize_file, dtype=int)
    #     if len(rindices) == data.expression.shape[1]:
    #         expr = data.expression[:, rindices]
    #     else:
    #         print("Random indices number do not match with current number of gx donors")
    #         raise

    ### Saikat shuffle on GTs
    #if args.shuffle:
    #    logger.warn("Shuffling genotype.") 
    #    gtcent_shuf = np.zeros_like(gtcent)
    #    gtnorm_shuf = np.zeros_like(gtnorm)
    #    for i in range(gtcent.shape[0]):
    #        idx = np.random.permutation(np.arange(0,gtcent.shape[1]))
    #        np.random.shuffle(idx)
    #        gtcent_shuf[i,:] = gtcent[i,idx]
    #        gtnorm_shuf[i,:] = gtnorm[i,idx]
    #    gtcent = gtcent_shuf
    #    gtnorm = gtnorm_shuf

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
snpinfo = comm.bcast(snpinfo, root = 0)
maf  = comm.bcast(maf, root = 0)
masklist = comm.bcast(masklist, root = 0)
maskcomp = comm.bcast(maskcomp, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()

if rank == 0: logger.debug("Computing JPA")
jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist)
jpa.compute()

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
    sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])

    if rank == 0: 
        if args.outprefix is None:
            args.outprefix = "out"    
        ohandle = Outhandler(args, snpinfo, geneinfo)
    
    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
    rr.compute(get_betas = True)

    if rank == 0:
        ohandle.write_rr_out(jpa, rr, prefix = "_it0")


    ###### New implementation of sparsity and null model calculation
    best_snp_indices = None
    gene_indices     = None
    if rank == 0:
        best_snp_indices = rr.select_best_SNPs(pval_thres = 1e-2, use_pvals = True)
        gene_indices     = rr.select_best_genes(rr.betas[best_snp_indices,:], n=5000)
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
                null_gene_indices = rr_null.select_best_genes(rr_null.betas, n=5000)
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

        np.savetxt("/cbscratch/franco/tejaas_output/tests/qnull_scores.txt", qnull) # for testing purposes
        logger.debug("TIME-->MPIcompute_null took {:g} seconds".format(time.time() - stime))
        
        logger.debug("---> computing rr for sparse")
        rr_sparse = RevReg(gtcent[best_snp_indices,:], expr, sb2_sparse, comm, rank, ncore, null = args.nullmodel, maf = maf[best_snp_indices], masks = newmasks)
        rr_sparse.compute_sparse(gene_indices, qnull, get_betas = False)
        ##### end

        if rank == 0:
            ohandle.write_rr_out(jpa, rr_sparse, selected_snps = best_snp_indices, selected_genes = gene_indices, prefix = "_it1")

if rank == 0: rr_time = time.time()
    
# Output handling only from master node // move it to module
# if rank == 0: 
#     if args.outprefix is None:
#         args.outprefix = "out"    
#     ohandle = Outhandler(args, snpinfo, geneinfo)
#     if args.jpa:
#         ohandle.write_jpa_out(jpa)
#     if args.rr:
#         ohandle.write_rr_out(jpa, rr)

if rank == 0: write_time = time.time()

if rank == 0:
    logger.info("File reading time: {:g} seconds".format(read_time - start_time))
    logger.info("JPA calculation time: {:g} seconds".format (jpa_time - read_time))
    logger.info("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    logger.info("Result writing time: {:g} seconds".format(write_time - rr_time))
    logger.info("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
