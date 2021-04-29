import numpy as np
#from scipy.stats import chi2
from scipy.stats import norm
import logging
from tejaas.utils.logs import MyLogger

logger = MyLogger(__name__)

def perm_null(GT, GX, sigmabeta2, sigmax2):
    nsnps = GT.shape[0]
    nsamples = GT.shape[1]
    ngenes = GX.shape[0]
    Rscore = np.zeros(nsnps)
    pvals  = np.zeros(nsnps)
    muQ    = np.zeros(nsnps)
    sigmaQ = np.zeros(nsnps)

    logger.debug("Number of SNPs: {:d}".format(nsnps))
    logger.debug("Number of samples: {:d}".format(nsamples))
    logger.debug("Number of genes: {:d}".format(ngenes))

    Yt = GX.T # shape N x G
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    S2 = np.square(S)
    for i in range(nsnps):
        S2mod = S2 + sigmax2[i] / sigmabeta2[i]
        W = np.dot(U, np.dot(np.diag(S2 / S2mod), U.T)) #/ sigmax2[i]
        Rscore[i] = np.sum(np.square(np.dot(U.T, GT[i,:])) * S2 / S2mod)
        pvals[i], muQ[i], sigmaQ[i] = pvals_perm(GT[i, :].reshape(1, -1), Rscore[i], W)

    return pvals, Rscore, muQ, sigmaQ

def pvals_perm(GT, R, W):
    mu2, mu4 = moment_data(GT)
    N = GT.shape[1]
    q11 = np.sum(W)
    q2  = np.sum(np.diag(W))
    muQ = mu2 * (N * q2 - q11) / (N - 1)

    v31 = - mu4 / (N - 1)
    v22 = v31 + (N * mu2 * mu2 / (N - 1)) #(N*(mu2**2) - mu4)/(N-1)
    v211 = - (v31 + v22) / (N - 2)
    v1111 = - 3 * v211 / (N - 3)

    q31 = np.dot(np.diag(W),np.sum(W,axis = 1))
    q4 = np.sum(np.square(np.diag(W)))
    q22 = np.sum(np.square(W))
    q211 = np.sum(np.square(np.sum(W,axis = 1)))

    sigma2 = v1111*(q11**2 - 2*q2*q11 - 4*q211 + 8*q31 + 2*q22 + q2**2 - 6*q4) + 2*v211*(q2*q11 + 2*q211 - 6*q31 - 2*q22 - q2**2 + 6*q4) + v22*(q2**2 + 2*q22 - 3*q4) + 4*v31*(q31 - q4) + mu4*q4

    sigma2 = sigma2 - muQ**2
    #mscale = sigma2 / muQ / 2.0
    #df = muQ / mscale
    #Rscaled = R / mscale
    #p = 1 - chi2.cdf(Rscaled, df)
    sigmaQ = np.sqrt(sigma2)
    p = 1 - norm.cdf(R, loc=muQ, scale=sigmaQ)
    return p, muQ, sigmaQ

def moment_data(GT):   #GT ixN
    GT2 = np.square(GT)
    GT4 = np.square(GT2)
    mu2 = np.mean(GT2)
    mu4 = np.mean(GT4)
    return mu2, mu4


def maf_null(GT, GX, sigmabeta2, sigmax2, maf):
    nsnps = GT.shape[0]
    nsamples = GT.shape[1]
    ngenes = GX.shape[0]

    Yt = GX.T # shape N x G
    U, S, Vt = np.linalg.svd(Yt, full_matrices=False)
    S2 = np.square(S)
    S2mod = S2 + sigmax2[0] / sigmabeta2[0]
    W = np.dot(U, np.dot(np.diag(S2 / S2mod), U.T))
    Rscore = np.zeros(nsnps)
    for i in range(nsnps):
        Rscore[i] = np.sum(np.square(np.dot(U.T, GT[i,:])) * S2 / S2mod)
    pvals, muQ, sigmaQ = pvals_maf(S, Rscore, sigmabeta2, sigmax2, maf, W)
    return pvals, Rscore, muQ, sigmaQ


def pvals_maf(S, R, sigmabeta2, sigmax2, maf, W):
    nsnps = R.shape[0]
    S2 = np.square(S)
    S2mod = S2 + sigmax2[0] / sigmabeta2[0]
    mu = np.sum(S2 / S2mod)
    sigma2 = 2 * np.sum(np.square(S2 / S2mod))

    x4 = gt_fourth_moment(maf)
    corr_maf = (x4 - 3) * np.sum(np.diag(np.square(W)))
    sigma2 += corr_maf

    #mscale = sigma / mu / 2.0
    #df = mu / mscale
    #Rscaled = R / mscale
    #p = 1 - chi2.cdf(Rscaled, df)

    sigma = np.sqrt(sigma2)
    p = 1 - norm.cdf(R, loc=mu, scale=sigma)
    return p, np.repeat(mu, nsnps), sigma


def gt_fourth_moment(maf):
    f0 = np.square(1 - maf)
    f1 = 2.0 * maf * (1 - maf)
    f2 = np.square(maf)
    mu = 2 * maf
    sig = np.sqrt(f1)
    x4 = f0 *  (-mu/sig) ** 4 + f1 * ((1 - mu)/sig)**4 + f2 * ((2 - mu)/sig)**4
    return x4
