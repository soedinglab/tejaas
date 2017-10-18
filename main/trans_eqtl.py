#!/usr/bin/env python

import numpy as np
from mpi4py import MPI
import ctypes
from scipy import stats
import os
import time
import argparse

from dosageparser import DosageParser
import qstat
import output

def parse_args():

    parser = argparse.ArgumentParser(description='Trans-eQTL joint analysis from all SNPs (TEJAS)')

    parser.add_argument('--genotype',
                        type=str,
                        dest='genotype_filename',
                        metavar='FILE',
                        help='input genotype file')

    parser.add_argument('--sample',
                        type=str,
                        dest='sample_filename',
                        metavar='FILE',
                        help='input fam file')

    parser.add_argument('--expression',
                        type=str,
                        dest='expression_filename',
                        metavar='FILE',
                        help='input expression file')

    parser.add_argument('--output',
                        type=str,
                        dest='output_fileprefix',
                        metavar='FILE',
                        help='output file prefix')

    parser.add_argument('--start',
                        type=int,
                        dest='startsnp',
                        help='starting SNP index')

    parser.add_argument('--end',
                        type=int,
                        dest='endsnp',
                        help='ending SNP index')

    opts = parser.parse_args()
    return opts


def read_expression(filename):
    gene_names = list()
    expression = list()
    with open(filename, 'r') as mfile:
        header = mfile.readline()
        donorids = header.strip().split()[1:]
        for line in mfile:
            linesplit = line.strip().split()
            expression.append(np.array(linesplit[1:], dtype=float))
            gene_names.append(linesplit[0])
    expression = np.array(expression)
    return donorids, expression, gene_names

def norm_binom(gt, freq):
    f = freq.reshape(-1, 1)
    gt = (gt - (2 * f)) / np.sqrt(2 * f * (1 - f))
    return gt

def cpvalcomp(geno, expr, qcal, shuffle=False):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('../lib/linear_regression.so', _path)
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
    if shuffle:
        np.random.shuffle(x)
        np.random.shuffle(y)
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


# ==================================================
# Start MPI calculation
# =================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()


if rank == 0:
    # this is the master
    # create a list of genotypes for sending to your slaves

    opts = parse_args()
    out_fileprefix = opts.output_fileprefix
    genotype_filename = opts.genotype_filename
    sample_filename = opts.sample_filename
    expression_filename = opts.expression_filename
    startsnp = opts.startsnp
    endsnp = opts.endsnp

    #genotype_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/dosages/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr21.gz"
    #sample_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/dosages/donor_ids.fam"
    #expression_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/gtex_wholeblood_normalized.expression.txt"

    start_time = time.time()

    ds = DosageParser(genotype_filename, sample_filename, startsnp, endsnp)
    dosage = ds.dosage
    snpinfo = ds.snpinfo
    donorids = ds.sample_id
    nsnps = ds.nsnps
    nsample = ds.nsample

    sampleids, expression, gene_names = read_expression(expression_filename)
    
    choose_ids = [x for x in sampleids if x in donorids]
    dosage_indices = [i for i, x in enumerate(donorids)  if x in choose_ids]
    exprsn_indices = [i for i, x in enumerate(sampleids) if x in choose_ids]
    
    geno = dosage[:, dosage_indices]
    freq = np.array([x.freq for x in snpinfo])

    print ("Completed data reading")

    # These are the inputs to cpvalcomp
    geno = norm_binom(geno, freq)
    expr = expression[:, exprsn_indices]
    #qcal = qstat.q_calibrate(G = 5000, S = 10000, p0 = 0.01)
    qcal = qstat.q_calibrate()

    print ("Completed Qstat calibration")

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

# ==================================
# Everything sent. Now do the calculations
# ==================================
pvals, qscores, pqvals, gene_indices = cpvalcomp(slave_geno, expr, qcal, shuffle=True)

pvals   = comm.gather(pvals,   root = 0)
qscores = comm.gather(qscores, root = 0)
pqvals  = comm.gather(pqvals,  root = 0)
gene_indices = comm.gather(gene_indices,  root = 0)

if rank == 0:
    pvals = np.vstack(pvals)
    np.save(out_fileprefix, pvals)
    qscores = np.concatenate(qscores)
    pqvals  = np.concatenate(pqvals)
    gene_indices_list = list()
    for x in gene_indices:
        gene_indices_list += x
    output.write(out_fileprefix, snpinfo, gene_names, pvals, qscores, pqvals, gene_indices_list)
    print ("Job completed in {:g} seconds".format(time.time() - start_time))
else:
    assert gene_indices is None
    assert qscores is None
    assert pvals   is None
    assert pqvals  is None

#print ("Done")
