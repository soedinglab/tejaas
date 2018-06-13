import numpy as np
import ctypes
import os

def crevreg(geno, expr, sb2, null, maf):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('../lib/reverse_regression.so', _path)
    cqscore = clib.driver
    cqscore.restype = ctypes.c_int
    cqscore.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        ctypes.c_int,
                        np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                        ctypes.c_int,
                        ctypes.c_int,
                        ctypes.c_int,
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
    sig2 = np.zeros(nsnps)

    success = cqscore(x, y, sb2, null, maf, ngenes, nsnps, nsamples, R, p, mu, sig2)
    sigma = np.sqrt(sig2)
    return p, R, mu, sigma


def maf_null(gt, gx, sb2, sx2, maf):
    p, R, mu, sigma = crevreg(gt, gx, sb2, 1, maf)
    return p, R, mu, sigma


def perm_null(gt, gx, sb2, sx2):
    maf = np.zeros(gt.shape[0])
    p, R, mu, sigma = crevreg(gt, gx, sb2, 2, maf)
    return p, R, mu, sigma
