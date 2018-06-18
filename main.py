#/usr/bin/env python

import numpy as np
import logging
import time
import itertools
import sqlite3

import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from utils import iotools
from qstats.jpa import JPA
from qstats.revreg import RevReg
from iotools.data import Data
import iotools.outhandler as outhandler
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

if rank == 0:
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
    maf = np.array([x.maf for x in snpinfo])

gtnorm = comm.bcast(gtnorm, root = 0)
gtcent = comm.bcast(gtcent, root = 0)
expr   = comm.bcast(expr,  root = 0)
maf  = comm.bcast(maf, root = 0)
comm.barrier()

if rank == 0: read_time = time.time()

jpa = JPA(gtnorm, expr, comm, rank, ncore, args.jpa)
jpa.compute()

if rank == 0: jpa_time = time.time()

# Select the SNPs with JPA score above threshold for RevReg
if args.jpa and args.rr:
    select = None
    if rank == 0:
        select = np.where(jpa.scores > args.jpacut)[0]
        logger.debug('JPA threshold: {:g}, and number of SNPs retained: {:d}'.format(args.jpacut, select.shape[0]))

    qselect = comm.bcast(select, root = 0)
    comm.barrier()
    geno = geno[qselect, :]


if args.rr:
    sigbeta2 = np.repeat(args.sigmabeta ** 2, gtnorm.shape[0])
    if args.nullmodel == 'maf':
        rr = RevReg(gtnorm, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf)
    elif args.nullmodel == 'perm':
        rr = RevReg(gtcent, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel)
    rr.compute()
    
    
    
    

if rank == 0: rr_time = time.time()

# Output handling only from master node // move it to module
if rank == 0:
    pvals = jpa.pvals
    logger.debug('Pval matrix size: {:d} x {:d}'.format(pvals.shape[0], pvals.shape[1]))
    if(args.outprefix is None):
        logger.debug('No outfile prefix given, using default names\n')
        args.outprefix = "out"
    #args.outprefix = args.outprefix + 
    #dbwriter = outhandler.DBwriter(snpinfo, geneinfo, jpa.pvals, args.outprefix + ".db" )
    #dbwriter.write()

    if args.jpa:
        jpascores = jpa.scores
        jpawriter = outhandler.jpaOutWriter(snpinfo, jpascores, args.outprefix + "_jpa.txt")
        jpawriter.write()
        logger.debug('Scores size: {:d}'.format(jpascores.shape[0]))
    if args.rr:
        rrscores = rr.scores
        indices  = rr.pvals < args.psnpcut
        indices  = list(itertools.compress(range(len(indices)), indices))
        mu = np.mean(rr.null_mu)
        sigma = np.mean(rr.null_sigma)
        logger.debug('Mean of RR scores: {:g}, Mean of RR null: {:g}\n'.format(np.mean(rrscores), mu))
        logger.debug('Variance of RR scores: {:g}, Variance of RR null: {:g}\n'.format(np.std(rrscores), sigma))
        
        selected_snps = [snpinfo[i] for i in indices]
        #selected_pvals = list(itertools.compress(indices, g))        
        

        dbwriter = outhandler.DBwriter(selected_snps, geneinfo, jpa.pvals[indices,:], args.outprefix + ".db" )
        dbwriter.write()
        rrwriter = outhandler.rrOutWriter(selected_snps, [rr.pvals[i] for i in indices], [rr.scores[i] for i in indices], [rr.null_mu[i] for i in indices], [rr.null_sigma[i] for i in indices], args.outprefix + "_rr.txt")
        rrwriter.write()
        db      = sqlite3.connect(args.outprefix + ".db")
        cursor  = db.cursor()
        cursor.execute('''SELECT geneid, snpid, pval FROM pvals WHERE pval < (?) ORDER BY geneid''',(args.pgenecut,))  #argument for the ? mark needs to be given in tuple
        f = open(args.outprefix + "_gene_snp_list.txt", "w")
        f.write("geneid\tsnpid\tpval\n")
        for row in cursor:
            f.write(row[0] + "\t" + row[1] + "\t" + str(row[2]) + "\n")
        f.close()

if rank == 0: write_time = time.time()

if rank == 0:
    print ("File reading time: {:g} seconds".format(read_time - start_time))
    print ("JPA calculation time: {:g} seconds".format (jpa_time - read_time))
    print ("RR calculation time: {:g} seconds".format (rr_time - jpa_time))
    print ("Result writing time: {:g} seconds".format(write_time - rr_time))
    print ("Total execution time: {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
