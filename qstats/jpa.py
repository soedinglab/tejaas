import os
import numpy as np
from scipy import special
from scipy import stats
import ctypes
from utils import mpihelper
from utils.logs import MyLogger

class JPA:

    def __init__(self, x, y, comm, rank, ncore, qcalc):
        self.gt = x
        self.gx = y
        self._pvals = None
        self._qscores = None
        self.rank = rank
        self.comm = comm
        self.ncore = ncore
        self.qcalc = qcalc
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

    # @scores.setter
    # def scores(self, scores):
    #     self._qscores = scores

    def apply_mask(self, masks, snpmasks):
        newqscores = list()
        for j, mask in enumerate(masks):
            for i in snpmasks[j]:
                pvals_crop = np.delete(self.pvals[i,:], mask, axis = 0)
                if len(mask) > 0:
                    self._pvals[i,:][np.array(mask)] = 1
                # print("Iterating on SNP {:d} with {:d}/{:d} pvals".format(i, len(pvals_crop), len(self.pvals[i,:])))
                newqscore = self.jpascore(pvals_crop)
                newqscores.append(newqscore)
        self._qscores = np.array(newqscores)

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


    def clinreg(self, geno, expr, nrow):
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
        xsize = x.shape[0]
        nsnps = geno.shape[0]
        nsample = geno.shape[1]
        ngene = expr.shape[0]
        fstat = np.zeros(nsnps * ngene)
        success = cfstat(x, y, nsnps, ngene, nsample, fstat)
        res = 1 - stats.f.cdf(fstat, 1, nsample-2)
        res = res.reshape(nsnps, ngene)
        pvals = np.zeros((nrow, ngene), dtype=np.float64)
        pvals[:nsnps, :] = res
        return pvals


    def slavejob(self, geno, expr, nmax, jpa = True):
        nsnps = geno.shape[0]
        pvals = self.clinreg(geno, expr, nmax)
        if jpa:
            qscores = np.array([self.jpascore(pvals[i,:]) for i in range(nsnps)])
        else:
            qscores = np.zeros(nsnps)
        # return pvals, qscores
        return pvals.reshape(-1,), qscores


    def mpicompute(self):
        if self.rank == 0:
            # this is the master
            # create a list of genotypes for sending to your slaves
            geno = mpihelper.split_genotype(self.gt, self.ncore)
            expr = self.gx
            snp_per_node = [x.shape[0] for x in geno]
            print(snp_per_node)
            # nmax = max(snp_per_node)
        else:
            geno = None
            expr = None
            # nmax = None
            snp_per_node = None
        
        slave_geno = self.comm.scatter(geno, root = 0)
        expr = self.comm.bcast(expr, root = 0)
        nmax = self.comm.scatter(snp_per_node, root = 0)
        # nmax = self.comm.bcast(nmax, root = 0)
        self.comm.barrier()
        
        # ==================================
        # Data sent. Now do the calculations
        # ==================================
        if self.qcalc:
            pvals, qscores = self.slavejob(slave_geno, expr, nmax, jpa = True)
        else:
            pvals, qscores = self.slavejob(slave_geno, expr, nmax, jpa = False)


        recvbuf = None
        if self.rank == 0:
            self.logger.debug("Number of SNPs sent to each slave: " + ", ".join(["{:d}".format(x) for x in snp_per_node])) #print (recvbuf.shape)
            # recvbuf = np.zeros([self.gt.shape[0], self.gx.shape[1]], dtype=np.float64)
            columns = self.gx.shape[0]
            flat_sizes = np.array([rows * columns for rows in snp_per_node])
            total_rows = sum(np.array(snp_per_node))
            recvbuf = np.zeros(sum(flat_sizes), dtype=np.float64)
        else:
            flat_sizes = None


        # self.comm.Gather(pvals, recvbuf, root=0)
        
        self.comm.Gatherv(sendbuf=pvals, recvbuf=(recvbuf, flat_sizes), root = 0)

        if self.rank == 0:
            newarr = recvbuf.reshape(total_rows, columns)
            self._pvals = newarr

        qscores = self.comm.gather(qscores, root = 0)

        if self.rank == 0:
            # self.logger.debug("Shape of gathered pvals array: " + " x ".join(str(x) for x in recvbuf.shape)) #print (recvbuf.shape)
            # self.logger.debug("Number of SNPs: {:d}, Sum of rows in pvals array: {:d}".format(self.gt.shape[0], recvbuf.shape[0] * recvbuf.shape[1]))
            #self.logger.debug("Last 10 unused pvals are " + " ".join("{:g}".format(x) for x in recvbuf[0][-1, -10:]))
            # self._pvals = np.zeros((self.gt.shape[0], self.gx.shape[0]))
            # offset = 0
            # for i in range(self.ncore):
            #     nsnp = snp_per_node[i]
            #     self._pvals[offset:offset+nsnp, :] = recvbuf[i][:nsnp, :]
            #     self.logger.debug("Adding next {:d} SNPs from index {:d}".format(nsnp, offset))
            #     offset += nsnp
            self._qscores = np.concatenate(qscores)
        else:
            assert qscores is None
            assert recvbuf is None

        return
            

    def compute(self):
        if self.mpi:
            self.mpicompute()
        else:
            pvals, qscores = self.slavejob(self.gt, self.gx, self.gt.shape[0])
            self._pvals = pvals
            self._qscores = qscores
        return
