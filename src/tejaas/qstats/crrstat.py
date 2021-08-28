import numpy as np
import ctypes
import os

from tejaas.utils import project

def crevreg(geno, expr, sb2, null, maf):
    clib = project.get_clib('reverse_regression')
    cqscore = clib.qscore
    cqscore.restype = ctypes.c_bool
    cqscore.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int,
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
                       ]
    x = geno.reshape(-1,)
    y = expr.reshape(-1,)
    xsize = x.shape[0]
    nsnps = geno.shape[0]
    nsamples = geno.shape[1]
    ngenes = expr.shape[0]
    R = np.zeros(nsnps)
    p = np.zeros(nsnps)
    mu = np.zeros(nsnps)
    sigma = np.zeros(nsnps)
    success = cqscore(x, y, sb2, ngenes, nsnps, nsamples, null, maf, R, p, mu, sigma)
    return p, R, mu, sigma

def crrbetas(geno, expr, sb2):
    clib = project.get_clib('reverse_regression')
    cbetas = clib.betas
    cbetas.restype = ctypes.c_bool
    cbetas.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int,
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
                       ]
    x = geno.reshape(-1,)
    y = expr.reshape(-1,)
    nsnps = geno.shape[0]
    nsamples = geno.shape[1]
    ngenes = expr.shape[0]
    B = np.zeros(ngenes * nsnps)
    success = cbetas(x, y, sb2, ngenes, nsnps, nsamples, B)
    return B

def maf_null(gt, gx, sb2, sx2, maf):
    p, R, mu, sigma = crevreg(gt, gx, sb2, 1, maf)
    return p, R, mu, sigma

def perm_null(gt, gx, sb2, sx2):
    maf = np.zeros(gt.shape[0])
    p, R, mu, sigma = crevreg(gt, gx, sb2, 2, maf)
    return p, R, mu, sigma

def no_null(gt, gx, sb2, sx2):
    maf = np.zeros(gt.shape[0])
    p, R, mu, sigma = crevreg(gt, gx, sb2, 0, maf)
    return p, R, mu, sigma
