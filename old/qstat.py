#!/bin/usr/env python

import numpy as np
import collections
from statsmodels.distributions.empirical_distribution import ECDF
from scipy import special
from scipy.stats import expon

QCAL_FIELDS = ['join', 'lam', 'coef', 'ecdf']
class Q_Cal(collections.namedtuple('_Q_Cal', QCAL_FIELDS)):
    __slots__ = ()

def qscore(pvals):
    p = np.sort(pvals)
    n = p.shape[0]
    kmax = min(100, n)
    krange = [i + 1 for i in range(kmax)]
    digamma_n1 = special.digamma(n + 1)
    z = - ( np.log(p[:kmax]) - (special.digamma(krange) - digamma_n1) )
    zsum = np.cumsum(z)
    res = np.max(zsum)
    return res

def q_calibrate(qemp, p0=0.001):
    #qneg = np.load(qemp)
    qneg = qemp
    qneg_sorted = qneg[np.argsort(qneg)]
    S = qneg_sorted.shape[0]
    ntop = int(p0 * S)
    #cdf = np.arange(S)/float(S - 1)
    match = S - ntop - 1
    print(qneg)
    qjoin = qneg_sorted[ match  ]
    lam = ntop / np.sum(qneg_sorted[match:] - qjoin) # maximum likelihood estimate of exponential distribution of Q
    qcal = Q_Cal(join = qjoin,
                 lam  = lam,
                 coef = ntop / (S - 1),
                 ecdf = ECDF(qneg_sorted[::-1]))
    return qcal

def p_qscore(qscore, qcal):
    if qscore > qcal.join:
        res = qcal.coef * np.exp( - (qscore - qcal.join) * qcal.lam)
    else:
        res = 1 - qcal.ecdf(qscore)
    return res

## ==== old Q-score calibration using uniform distribution of p-values 
##      (deprecated after Anubhav showed that gene expressions are correlated) 
## ====
## def q_calibrate(G = 20000, S = 100000, p0 = 0.001):
##     qneg = [qscore(np.random.rand(G)) for i in range(S)]
##     qneg = np.array(qneg)
##     qneg_sorted = qneg[np.argsort(-qneg)]
##     n1 = int(p0 * S)
##     #lam = S / np.sum(qneg) # maximum likelihood estimate of exponential distribution of Q
##     loc, scale = expon.fit(qneg, floc = 0)
##     lam = 1 / scale
##     qcal = Q_Cal(qmax = qneg_sorted[n1],
##                  lam  = lam,
##                  p0 = p0,
##                  ecdf = ECDF(qneg_sorted))
##     return qcal
