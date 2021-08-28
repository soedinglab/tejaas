#/usr/bin/env python
'''
Tejaas: Discover Trans-eQTLs!

This initiates the MPI calls for Tejaas.
It is a wrapper for all the options used in Tejaas.
Currently, the following arguments are processed:
1) args.jpa
    This is very similar to the CPMA-score of Brynedal et. al.
    Requires a file containing null JPA-scores.
    The user must provide a filepath using the command: --jpanull filepath
        - If file exists, then it uses the file.
        - If file does not exist, then it creates a new null file at the specified filepath.
    Performs JPA-score analysis and calculate p-values using the null JPA-scores.
3) args.rr
    Performs RR-score analysis and calculate p-values from an analytical null model.
    Two analytical null models are implemented:
        a) --null maf
            MAF null (old deprecated version, not discussed in the manuscript)
            Requires 1000G / population MAF file, see documentation.
        b) --null perm
            PERM null (null model obtained by permuting SNPs)
'''

import numpy as np
import time
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from tejaas.utils.args import Args
from tejaas.utils.logs import MyLogger
from tejaas.utils import project
from tejaas.qstats.jpa import JPA
from tejaas.qstats.jpa_null import JPANULL
from tejaas.qstats.revreg import RevReg

from tejaas.iotools.data import Data
from tejaas.iotools.outhandler import Outhandler
from tejaas.iotools import readmaf
from tejaas.iotools import readqnull

from tejaas.tester import UnittestTester
from tejaas.tests.test_output_files import TestOutputFiles

def run_unittests(comm, args):
    tester = UnittestTester(TestOutputFiles)
    tester.execute()
    del tester
    return

def run_tejaas(comm, args):

    # ==================================
    # Start MPI calculation
    # ==================================
    rank = comm.Get_rank()
    ncore = comm.Get_size()
    if rank == 0: start_time = time.time()


    # ==================================
    # Input Processing
    # ==================================
    logger = MyLogger(__name__)

    # List of variables that are broadcast over all slave nodes
    gtcent       = None   # Centered genotype (x - mu). Note: Not scaled (divided) by standard deviation. Dimension I x N.
    gtnorm       = None   # Centered and scaled genotype (x - mu) / sigma. Dimension I x N.
    expr         = None   # Expression matrix. Dimension G x N.
    tgene_gtnorm = None   # Centered and scaled genotype without KNN correction. Dimension I x N.
    tgene_expr   = None   # Expression matrix for finding target genes after RR-score. Dimension G x N.
    masklist     = None   # List of gene indices masked for each SNP. List of length I. Each element contains a list of indices for the cis-genes of a candidate SNP.
    maskcomp     = None   # List of CisMasks. Each element is a collection of (a) gene indices to be masked, and (b) all SNP indices which require this mask.
    maf          = None   # List of MAF of each SNP as observed in the sample (or read separately from the population if file is provided). Length I.
    if rank == 0:
        logger.debug("Using {:d} cores".format(ncore))
        data = Data(args)
        data.load()
        gtcent       = data.geno_centered
        gtnorm       = data.geno_normed
        expr         = data.expression
        tgene_gtnorm = data.tgene_geno_normed 
        tgene_expr   = data.tgene_expression
        masklist     = data.cismasks_list
        maskcomp     = data.cismasks_comp
        snpinfo      = data.snpinfo  # slaves don't need this. Only needed for maf + outhandler in master node.
        geneinfo     = data.geneinfo # slaves don't need this. Only needed for outhandler in master node.
        maf          = readmaf.load(snpinfo, maf_file = args.maf_file)
        logger.debug("After prefilter: {:d} SNPs and {:d} genes in {:d} samples".format(gtcent.shape[0], expr.shape[0], gtcent.shape[1]))
    gtcent       = comm.bcast(gtcent, root = 0)
    gtnorm       = comm.bcast(gtnorm, root = 0)
    expr         = comm.bcast(expr,  root = 0)
    tgene_gtnorm = comm.bcast(tgene_gtnorm, root = 0)
    tgene_expr   = comm.bcast(tgene_expr, root = 0)
    masklist     = comm.bcast(masklist, root = 0)
    maskcomp     = comm.bcast(maskcomp, root = 0)
    maf          = comm.bcast(maf, root = 0)

    comm.barrier()
    if rank == 0: read_time = time.time()


    # ==================================
    # Calculate the null JPA-scores
    # ==================================
    if args.jpa_calc_null:
        if rank == 0: logger.debug("Calculating null values for JPA")
        jpa_null = JPANULL(gtnorm, expr, comm, rank, ncore, args.jpanull_file, niter=args.jpanull_iter, seed=args.seed)
        jpa_null.compute()
    if rank == 0: jpanull_time = time.time()


    # ==================================
    # Calculate JPA-score
    # with empirical p-values from the null file
    # ==================================
    if args.jpa:
        if rank == 0: logger.debug("Computing JPA")
        jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa, masklist, get_pvals = True, qnull_file = args.jpanull_file)
        jpa.compute()

        ### Calculate target genes
        teqtl_id = None
        if rank == 0: teqtl_id = np.where(jpa.jpa_pvals <= args.psnpcut)[0]
        teqtl_id = comm.bcast(teqtl_id, root = 0)

        tgjpa = JPA(gtnorm[teqtl_id, :], expr, comm, rank, ncore, False, None, get_pvals = False, statmodel = 'fstat')
        tgjpa.compute()

    if rank == 0: jpa_time = time.time()


    # ==================================
    # Calculate the RR-score
    #   - with perm-null (default)
    #   - with MAF-null (old and deprecated, kept for legacy)
    #
    # Find target genes 
    #   - Use KNN corrected gene expression and SNPs to find SNP-gene association p-values
    #   - Use separate confounder-corrected gene expression to find SNP-gene association p-values
    #   - Run JPA only on the trans-eQTLs (discovered with user-provided cutoff) with following options:
    #       a) qcalc = False // JPA-score is not calculated
    #       b) masklist = None // do not mask any gene. we also want cis-associations if any.
    #       c) get_pvals = False // p-values for significance of JPA-score is not calculated
    #       d) statmodel = 'fstat' // Use F-statistic as used in MatrixEQTL. Default is Z-statistic as used in CPMA.
    #
    # TO-DO:
    #  - The theory of MAF-null is not included in the manuscript, create a document
    # ==================================
    if args.rr:
        sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])
        if rank == 0: logger.debug("Sigma_beta fixed to {:g}".format(args.sigmabeta))
        if args.nullmodel == 'maf':
            rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
        elif args.nullmodel == 'perm':
            rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf, masks = maskcomp)
        rr.compute(get_betas = False)  # Set get_betas = True to obtain the coefficients of multiple regression (stored in rr.betas)

        teqtl_id = None
        if rank == 0: teqtl_id = np.where(rr.pvals <= args.psnpcut)[0]
        teqtl_id = comm.bcast(teqtl_id, root = 0)

        tgknn = JPA(gtnorm[teqtl_id, :], expr, comm, rank, ncore, False, None, get_pvals = False, statmodel = 'fstat')
        tgknn.compute()

        if args.usefdr:
            # add masks
            tgjpa = JPA(tgene_gtnorm[teqtl_id, :], tgene_expr, comm, rank, ncore, False, [masklist[i] for i in teqtl_id], get_pvals = False, statmodel = 'fstat', target_fdr=0.5)
            tgjpa.compute()
        else:
            tgjpa = JPA(tgene_gtnorm[teqtl_id, :], tgene_expr, comm, rank, ncore, False, None, get_pvals = False, statmodel = 'fstat')
            tgjpa.compute()

    if rank == 0: rr_time = time.time()


    # ==================================
    # Write output
    # ==================================
    if rank == 0: 
        ohandle = Outhandler(args, snpinfo, geneinfo)
        if args.jpa:
            ohandle.write_jpa_out(jpa, tgjpa, teqtl_id)
        if args.rr:
            ohandle.write_rr_out(rr, tgknn, tgjpa, teqtl_id, write_betas = False)
    if rank == 0: write_time = time.time()


    # ==================================
    # Log the run times
    # ==================================
    if rank == 0:
        logger.info("File reading time: {:g} seconds".format(read_time - start_time))
        logger.info("JPA null calculation time: {:g} seconds".format (jpanull_time - read_time))
        logger.info("JPA calculation time: {:g} seconds".format (jpa_time - jpanull_time))
        logger.info("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
        logger.info("Result writing time: {:g} seconds".format(write_time - rr_time))
        logger.info("Total run time: {:g} seconds".format(time.time() - start_time))

def main():

    # ==================================
    # Start MPI calculation
    # ==================================
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncore = comm.Get_size()

    # ==================================
    # Input Processing
    # ==================================
    args = Args(comm, rank)

    if args.maketest:
        run_unittests(comm, args)
    else:
        run_tejaas(comm, args)

    MPI.Finalize()


if __name__ == "__main__":
    main()
