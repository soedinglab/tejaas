import numpy as np
from utils.containers import SnpInfo

def single_snp_permute(nsnp = 10, nsample = 338, fmin = 0.1, maketest = False):

    if maketest:
        nsample = 338
        fmin = 0.1
        nfreq = np.array([[279,  54,   5]])
    else:
        mafratios = np.array([(1 - fmin)**2, 2 * fmin * (1 - fmin), fmin**2])
        nfreq  = np.random.multinomial(nsample, mafratios, size=1)

    f1 = np.repeat(0, nfreq[0][0])
    f2 = np.repeat(1, nfreq[0][1])
    f3 = np.repeat(2, nfreq[0][2])
    x  = np.concatenate((f1,f2,f3))

    dosage = np.zeros((nsnp, nsample))
    snpinfo = list()
    for i in range(nsnp):
        dosage[i, :] = np.random.permutation(x)
        this_snp = SnpInfo(chrom      = 1,
                           bp_pos     = i,
                           varid      = 'rs{:d}'.format(i),
                           ref_allele = 'A',
                           alt_allele = 'T',
                           maf        = fmin)
        snpinfo.append(this_snp)
         
    maf2d = np.repeat(fmin, nsnp).reshape(-1, 1)
    gtnorm = (dosage - (2 * maf2d)) / np.sqrt(2 * maf2d * (1 - maf2d))
    gtcent = dosage - np.mean(dosage, axis = 1).reshape(-1, 1)

    return snpinfo, gtnorm, gtcent

def multiple_snp_permute(nsnp = 10, nsample = 338, fmin = 0.1, fmax = 0.9):

    f = np.random.uniform(fmin, fmax, nsnp)
    dosage = np.zeros((nsnp, nsample))
    snpinfo = list()
    for i in range(nsnp):
        dosage[i, :] = np.random.binomial(2, f[i], nsample)
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


#def permute(gx, nsnp = 10, ntrans = 10, fmin = 0.1, fmax = 0.9, maketest = False):
#    return snpinfo, gtnorm, gtcent
