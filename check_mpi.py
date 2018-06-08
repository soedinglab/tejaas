#!/usr/bin/env python

import numpy as np
from mpi4py import MPI
import ctypes
from scipy import stats
from sklearn import linear_model
import os
import time
import qstat

def cpvalcomp(geno, expr, qcal):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('lib/linear_regression.so', _path)
    cfstat = clib.fit
    cfstat.restype = ctypes.c_int
    cfstat.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                       ctypes.c_int,
                       ctypes.c_int,
                       ctypes.c_int,
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
                      ]

    x = geno.reshape(-1,)
    y = expr.reshape(-1,)
    xsize = x.shape[0]
    nsnps = geno.shape[0]
    nsample = geno.shape[1]
    ngene = expr.shape[0]
    fstat = np.zeros(nsnps * ngene)
    success = cfstat(x, y, nsnps, ngene, nsample, fstat)
    pvals = 1 - stats.f.cdf(fstat, 1, nsample-2)
    pvals = pvals.reshape(nsnps, ngene)
    qscores = np.array([qstat.qscore(pvals[i,:]) for i in range(nsnps)])
    pqvals  = np.array([qstat.p_qscore(qscores[i], qcal) for i in range(nsnps)])
    gene_indices = list()
    for snp in range(nsnps):
        if pqvals[snp] < 0.05:
            gene_indices.append(np.where(pvals[snp, :] < 0.05)[0])
        else:
            gene_indices.append(np.array([], dtype=int))
    return pvals, qscores, pqvals, gene_indices

nsnps = 100
ngene = 23000
nsample = 338

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

if rank == 0:
    # this is the master
    # create a list of genotypes for sending to your slaves

    geno = np.random.rand(nsnps * nsample).reshape(nsnps, nsample)
    expr = np.random.rand(ngene * nsample).reshape(ngene, nsample)

    start_time = time.time()
    qemp = np.random.normal(0, 1, 50000)
    print(qemp)
    qcal = qstat.q_calibrate(qemp)
    #qcal = qstat.q_calibrate(G = 5000, S = 10000, p0 = 0.01)

    maxsnp = int(nsnps / ncore)
    offset = 0
    data = [None for i in range(ncore)]
    for dest in range(ncore-1):
        start = offset
        end = offset + maxsnp
        data[dest] = geno[start:end, :]
        offset += maxsnp

    data[ncore-1] = geno[offset:nsnps, :]

else:
    data = None
    expr = None
    qcal = None

slave_geno = comm.scatter(data, root = 0)
expr = comm.bcast(expr, root = 0)
qcal = comm.bcast(qcal, root = 0)
comm.barrier()

# Now perform the following tasks at every node, before we start gathering the results
pvals, qscores, pqvals, gene_indices = cpvalcomp(slave_geno, expr, qcal)

pvals   = comm.gather(pvals,   root = 0)
qscores = comm.gather(qscores, root = 0)
pqvals  = comm.gather(pqvals,  root = 0)
gene_indices = comm.gather(gene_indices,  root = 0)

if rank == 0:
    pvals = np.vstack(pvals)
    qscores = np.concatenate(qscores)
    pqvals  = np.concatenate(pqvals)
    genes = list()
    for x in gene_indices:
        genes += x
    #gene_indices = np.concatenate(gene_indices)
    #output.write(filepath, snpinfo, gene_names, pvals, qscore, p_qscore, gene_indices)
        print (x)
    print (pvals.shape)
    print (time.time() - start_time)
else:
    assert gene_indices is None
    assert qscores is None
    assert pvals   is None
    assert pqvals  is None

#print("Hello MPI from %d of %d" % (rank, ncore))
