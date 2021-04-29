import numpy as np
from scipy.linalg import eigh
from scipy import special
from scipy import stats

from tejaas.qstats.zstats import ZSTATS
from tejaas.utils.logs import MyLogger
from tejaas.utils import mpihelper

class JPANULL:


    def __init__(self, x, y, comm, rank, ncore, outfile, niter=100000, seed=None):
        self.gt = x
        self.gx = y
        self._niter = niter
        self.outfile = outfile
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.mpi = False
        if seed is not None:
            np.random.seed(seed)
        if self.ncore > 1:
            self.mpi = True
        self.logger = MyLogger(__name__)


    def write_qscores(self):
        qmax = np.max(self._qscores)
        valid_qscores = self._qscores[np.where(self._qscores < qmax)]
        with open(self.outfile, 'w') as fout:
            for qnull in valid_qscores:
                fout.write(f"{qnull}\n")


    def jpascore(self, pvals):
        min_nonzero = np.min(pvals[np.nonzero(pvals)])
        pvals[pvals == 0] = min_nonzero
        p = np.sort(pvals)
        n = p.shape[0]
        kmax = min(100, n)
        krange = [i + 1 for i in range(kmax)]
        digamma_n1 = special.digamma(n + 1)
        z = - ( np.log(p[:kmax]) - (special.digamma(krange) - digamma_n1) )
        zsum = np.cumsum(z)
        res = np.max(zsum)
        return res


    def slavejob(self, W, Q, Zmean, n):
        self.logger.debug('Rank {:d} calculating {:d} null JPA scores'.format(self.rank, n))
        ngene = W.shape[0]
        pvals = np.zeros(n * ngene)
        for i in range(n):
            ngene = W.shape[0]
            zrand = np.random.normal(0, 1, size = ngene)
            znull = Zmean + np.einsum('ij, j, j', Q, np.sqrt(W), zrand)
            pvals[i*ngene : (i+1)*ngene] = 2.0 * (1 - stats.norm.cdf(np.abs(znull)))
        qnull = np.array([self.jpascore(pvals[i*ngene : (i+1)*ngene]) for i in range(n)])
        return pvals, qnull


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of N for sending to your slaves
            thisW = self._W
            thisQ = self._Q
            thisZmean = self._Zmean
            start, end = mpihelper.split_n(self._niter, self.ncore)
            nlist = [x - y for x, y in zip(end, start)]
        else:
            thisW = None
            thisQ = None
            thisZmean = None
            nlist = None
            slave_W = None
            slave_Q = None
            slave_Zmean = None
            slave_n = None

        slave_W = self.comm.bcast(thisW, root = 0)
        slave_Q = self.comm.bcast(thisQ, root = 0)
        slave_Zmean = self.comm.bcast(thisZmean, root = 0)
        slave_n = self.comm.scatter(nlist, root = 0)
        self.comm.barrier()

        if self.rank == 0: self.logger.debug("Broadcast W, Q and Zmean to the slave nodes")

        # ==================================
        # Data sent. Do the calculations
        # ==================================
        pvals, qscores = self.slavejob(slave_W, slave_Q, slave_Zmean, slave_n)

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


    def WQ_mpiwrap(self):
        '''
        Populates self._W, self._Q and self._Zmean
        Learns the null Zstats from the ZSTATS class
        and performs eigendecomposition in the master node.
        '''

        self._W = None
        self._Q = None
        self._Zmean = None

        if self.rank == 0: self.logger.debug("Computing Z-stats")
        zstats = ZSTATS(self.gt, self.gx, self.comm, self.rank, self.ncore)
        zstats.compute()
        
        if self.rank == 0:
            self.logger.debug("Computing W and Q")
            zscores = zstats.scores
            C = np.cov(zscores.T)
            # Numpy gives imaginary eigenvalues, use eigh from scipy
            # for decomposition of real symmetric matrix
            W, Q = eigh(C)
            self.logger.debug("Eigendecomposition done")
            # still some eigenvalues are negative. force them to zero if they are negligible. (!!!!!!!!!!!)
            # check if everything is ok
            #Wsparse = W.copy()
            #Wsparse[np.where(W < 0)] = 0
            #W = Wsparse
            W[np.where(W < 0)] = 0
            self.logger.debug("Forced negative eigenvalues to zero")
            #if not np.allclose(C, Q @ np.diag(W) @ Q.T):
            #    self.logger.error("Eigen vectors could not be forced to positive")
            #    exit
            #else:
            #    W = Wsparse
            #    self.logger.debug("Eigen vectors are forced to positive")
            Zmean = np.mean(zscores, axis = 0)
            self._W = W
            self._Q = Q
            self._Zmean = Zmean
        self.comm.barrier()


    def compute(self):
        self.WQ_mpiwrap()
        if self.rank == 0: self.logger.debug("Start MPI calculation")
        if self.mpi:
            self.mpicompute()
        else:
            pvals, qscores = self.slavejob(self._W, self._Q, self._Zmean, self._niter)
            self._pvals = pvals.reshape(self._niter, self._W.shape[0])
            self._qscores = qscores
        if self.rank == 0:
            self.logger.debug("Null JPA-scores calculated. Writing to file.")
            self.write_qscores()
        return
