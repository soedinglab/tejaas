import numpy as np
import ctypes
import os

def crevreg(geno, expr, sb2, null, maf):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('../lib/reverse_regression.so', _path)
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

def crevreg_corr(geno, expr, sb2, null, maf):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('../lib/reverse_regression.so', _path)
    cqscore = clib.qscore_corr
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
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library('../lib/reverse_regression.so', _path)
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

################################################
##### S values correction ######################
################################################
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

def fit_line(S):
    rank = len(S)
    start= int(rank*55/100)
    end  = int(rank*55/100+3)
    x    = np.arange(rank)[start:end].reshape(-1, 1)
    y    = S[start:end].reshape(-1, 1)
    reg = LinearRegression().fit(x, y)
    b  = reg.coef_[0]
    y0 = reg.intercept_
    return b, y0

def func(x, a, b):
    return a * x**2 + b

def correct_S(S):
    rank = len(S)
    x_vals = np.arange(rank/2, rank)
    b, y0 = fit_line(S)
    y_vals = y0 + b * x_vals
    err  = np.abs(y_vals - S[int(rank/2): rank])
    popt, pcov = curve_fit(func, x_vals, err/np.max(err), bounds=(0,[1,1]))
    corr = np.zeros(len(S))
    corr[:int(rank/2)] = S[:int(rank/2)]
    corr[int(rank/2):] = S[int(rank/2): rank] + err*func(x_vals, *popt)
    corr[-1] = S[-1]
    return corr, popt, err
##########################################


def perm_null(gt, gx, sb2, sx2):
    maf = np.zeros(gt.shape[0])
    ###### Original function call
    p, R, mu, sigma = crevreg(gt, gx, sb2, 2, maf)

    ######### S value correction - experimental ###########
    # p, R, mu, sigma = crevreg_corr(gt, gx, sb2, 2, maf)
    ########################################

    return p, R, mu, sigma

def no_null(gt, gx, sb2, sx2):
    maf = np.zeros(gt.shape[0])
    p, R, mu, sigma = crevreg(gt, gx, sb2, 0, maf)
    return p, R, mu, sigma
