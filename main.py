#/usr/bin/env python

import numpy as np
import logging
import time

import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from utils.logs import MyLogger
from utils import project
from qstats.jpa import JPA
from qstats.revreg import RevReg


nsnps = 10000
ngene = 23000
nsample = 338
start = 1
end = nsnps
#geno = np.random.rand(nsnps * nsample).reshape(nsnps, nsample)
expr = np.random.rand(ngene * nsample).reshape(ngene, nsample)

# ==================================================
# Start MPI calculation
# =================================================

start_time = time.time()

MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

args = Args(comm, rank)
logger = MyLogger(__name__)

dosage = None
geno = None
maf = None
if rank == 0:
    fmin = 0.4
    mafratios = np.array([(1 - fmin)**2, 2 * fmin * (1 - fmin), fmin**2])
    
    nfreq  = np.random.multinomial(nsample, mafratios, size=1)
    f1 = np.repeat(0, nfreq[0][0])
    f2 = np.repeat(1, nfreq[0][1])
    f3 = np.repeat(2, nfreq[0][2])
    x  = np.concatenate((f1,f2,f3))
    
    dosage = np.zeros((nsnps, nsample))
    for i in range(nsnps):
        #dosage[i, :] = np.random.permutation(x)
        dosage[i, :] = np.random.binomial(2, fmin, nsample)
        
    maf = np.repeat(fmin, nsnps)
    geno = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)

dosage = comm.bcast (dosage, root = 0)
geno = comm.bcast(geno, root = 0)
maf  = comm.bcast(maf,  root = 0)
comm.barrier()

jpa = JPA(geno, expr, comm, rank, ncore, args.jpa)
jpa.compute()

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
    if args.nullmodel == 'maf':
        f = maf.reshape(-1, 1)
        geno = (dosage - (2 * f)) / np.sqrt(2 * f * (1 - f))
    sigbeta2 = np.repeat(0.001 ** 2, geno.shape[0])
    rr = RevReg(geno, expr, sigbeta2, comm, rank, ncore, null = args.nullmodel, maf = maf)
    rr.compute()

# Output handling only from master node // move it to module
if rank == 0:
    pvals = jpa.pvals
#    logger.debug('Pval matrix size: {:d} x {:d}'.format(pvals.shape[0], pvals.shape[1]))
#    outhandler = OutputHandler(outprefix, compute_jpa, compute_rr)
#    outhandler.writedb(snpinfo, geneinfo, pvals)
    if args.jpa:
        jpascores = jpa.scores
#        outhandler.writejpa(snpinfo, jpascores)
#        logger.debug('Scores size: {:d}'.format(jpascores.shape[0]))
    if args.rr:
        rrscores = rr.scores
        mu = np.mean(rr.null_mu)
        sigma = np.mean(rr.null_sigma)
        logger.debug('Mean of RR scores: {:g}, Mean of RR null: {:g}\n'.format(np.mean(rrscores), mu))
        logger.debug('Variance of RR scores: {:g}, Variance of RR null: {:g}\n'.format(np.std(rrscores), sigma))
#        outhandler.writerr(snpinfo, rrscores)

if rank == 0:
    print ("Job completed in {:g} seconds".format(time.time() - start_time))

MPI.Finalize()
