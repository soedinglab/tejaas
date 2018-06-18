import numpy as np
from utils.containers import SnpInfo
from utils.logs import MyLogger

logger = MyLogger(__name__)

def permuted_dosage(expr, nsnp = 5000, fmin = 0.1, fmax = 0.9, maketest = False):

    f = np.random.uniform(fmin, fmax, nsnp)
    if maketest:
        f = np.repeat(0.1, nsnp)
    nsample = expr.shape[1]

    dosage = np.zeros((nsnp, nsample))
    snpinfo = list()
    for i in range(nsnp):
        if maketest:
            nfreq = np.array([[279,  54,   5]])[0]
        else:
            mafratios = np.array([(1 - f[i])**2, 2 * f[i] * (1 - f[i]), f[i]**2])
            nfreq  = np.random.multinomial(nsample, mafratios, size=1)[0]
        f1 = np.repeat(0, nfreq[0])
        f2 = np.repeat(1, nfreq[1])
        f3 = np.repeat(2, nfreq[2])
        x  = np.concatenate((f1,f2,f3))
        dosage[i, :] = np.random.permutation(x)
        this_snp = SnpInfo(chrom      = 1,
                           bp_pos     = i,
                           varid      = 'rs{:d}'.format(i),
                           ref_allele = 'A',
                           alt_allele = 'T',
                           maf        = f[i])
        snpinfo.append(this_snp)

    maf2d = f.reshape(-1, 1)
    gtnorm = (dosage - (2 * maf2d)) / np.sqrt(2 * maf2d * (1 - maf2d))
    gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)

    return snpinfo, gtnorm, gtcent


def expression(gt, gx, cfrac = 0.001):

    ntrans  = gt.shape[0]
    ngene   = gx.shape[0]
    nsample = gx.shape[1]

    liabilities = np.zeros((ngene, nsample))
    cindex = np.zeros((ntrans, ngene))                                           # Index matrix of gene / trans-eQTLs pairs
    nc  = np.random.gamma(ngene * cfrac, scale = 1.0, size = ntrans).astype(int) # number of target genes for each trans-eQTL
    for i in range(ntrans):
        ncausal = min(ngene, nc[i])                                              # do something better, trans-eQTLs cannot target all genes
        choose = np.random.choice(ngene, ncausal, replace = False)
        cindex[i, choose] = 1                                                    # mark these genes as causal

    gtarget = list()
    H2  = np.square(np.random.normal(2, 1, ngene))                               # do something better, target variance of sum(x_i * beta_i)
    H2 /= np.max(H2)
    for i in range(ngene):
        csnps = np.where(cindex[:, i] == 1)[0]
        if csnps.shape[0] > 0: # then we got a trans-eQTL
            beta = np.random.normal(0, 1, size = csnps.shape[0])
            beta *= np.sqrt( H2[i] / np.sum(np.square(beta)) )
            liabilities[i, :] = np.dot(gt[csnps, :].T, beta)
            gtarget.append(i)

    logger.debug("Created {:d} target genes".format(len(gtarget)))
    newGX = gx + liabilities
    newGX = (newGX - np.mean(newGX, axis = 1).reshape(-1, 1)) / np.std(newGX, axis = 1).reshape(-1, 1)
    return newGX, H2[np.array(gtarget)], nc
