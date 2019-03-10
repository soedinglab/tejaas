import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper
from utils.logs import MyLogger

class JPANULL:

    def __init__(self, W, Q, N, comm, rank, ncore):
        self._pvals = None
        self._qscores = None
        self._W = W
        self._Q = Q
        self._N = N
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.mpi = False
        if self.ncore > 1:
            self.mpi = True
        self.logger = MyLogger(__name__)
        

    @property
    def pvals(self):
        return self._pvals


    @property
    def scores(self):
        return self._qscores


    def jpascore(self, pvals):
        p = np.sort(pvals)
        n = p.shape[0]
        kmax = min(100, n)
        krange = [i + 1 for i in range(kmax)]
        digamma_n1 = special.digamma(n + 1)
        z = - ( np.log(p[:kmax]) - (special.digamma(krange) - digamma_n1) )
        zsum = np.cumsum(z)
        res = np.max(zsum)
        return res


    def slavejob(self, W, Q, n):
        self.logger.debug('Rank {:d} calculating for {:d} SNPs'.format(self.rank, n))
        ngene = W.shape[0]
        pvals = np.zeros(n * ngene)
        for i in range(n):
            ngene = W.shape[0]
            znull = np.random.normal(0, 1, size=ngene)
            pvals[i*ngene : (i+1)*ngene] = 1 - stats.norm.cdf(znull)
            #znull = zmu + Q * np.sqrt(W) * N(0,1)
            #pvals[i*ngene : (i+1)*ngene] = 1 - stats.chi2.cdf(znull)
        qnull = np.array([self.jpascore(pvals[i*ngene : (i+1)*ngene]) for i in range(n)])
        return pvals, qnull


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of N for sending to your slaves
            thisW = self._W
            thisQ = self._Q
            start, end = mpihelpher.split_n(self._N, self.ncore)
            nlist = [x - y for x, y in zip(end, start)]
        else:
            nlist = None
            slave_W = None
            slave_Q = None
            slave_n = None
        
        slave_W = self.comm.bcast(W, root = 0)
        slave_Q = self.comm.bcast(Q, root = 0)
        slave_n = self.comm.scatter(nlist, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Do the calculations
        # ==================================
        pvals, qscores = self.slavejob(slave_W, slave_Q, slave_n)

        # ==================================
        # Collect the results
        # ==================================
        recvbuf = None
        if self.rank == 0:
            self.logger.debug("Number of SNPs sent to each slave: " + ", ".join(["{:d}".format(x) for x in nlist]))
            ngene = self._W.shape[0]
            flat_sizes = np.array([n * ngene for n in nlist])
            recvbuf = np.zeros(sum(flat_sizes), dtype=np.float64)
        else:
            flat_sizes = None
        self.comm.Gatherv(sendbuf=pvals, recvbuf=(recvbuf, flat_sizes), root = 0)

        if self.rank == 0:
            self._pvals = recvbuf.reshape(sum(nlist), ngene)

        qscores = self.comm.gather(qscores, root = 0)

        if self.rank == 0:
            self._qscores = np.concatenate(qscores)
        else:
            assert qscores is None
            assert recvbuf is None
        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            pvals, qscores = self.slavejob(self._W, self._Q, self._N)
            self._pvals = pvals.reshape(self._N, self._W.shape[0])
            self._qscores = qscores
        return
