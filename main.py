#/usr/bin/env python

'''
    Distribute job to:
    a. jpa.py
    b. rr.py
    c. jpa-rr.py
    d. optim.py
'''

import numpy as np
#from mpi4py import MPI
import logging
import time
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

from utils.args import Args
from qstats.jpa import JPA
#from qstats.revreg import RevReg
#from qstats.optim import Optim

## Get the version from version.py without importing the package
#exec(compile(open('version.py').read(), 'version.py', 'exec'))
#
## Initialize logger
#logger = logging.getLogger('tejaas')
## set to DEBUG level
#logger.setLevel(10)
## and always display on console
#stderr_log_handler = logging.StreamHandler()
#logger.addHandler(stderr_log_handler)
#
#logger.info('Running TEJAAS v{:s}'.format(locals()['__version__']))
#logger.info('Genotype File: {:s}'.format(args.vcf_file))
#logger.info('Gene expression file: {:s}'.format(args.gx_file))
#logger.info('Method: {:s}'.format(args.method))
#logger.info('using SNPs {:d} to {:d}'.format(args.startsnp, args.endsnp))

nsnps = 100
ngene = 1000
nsample = 338
start = 1
end = 100
geno = np.random.rand(nsnps * nsample).reshape(nsnps, nsample)
expr = np.random.rand(ngene * nsample).reshape(ngene, nsample)

# ==================================================
# Start MPI calculation
# =================================================

start_time = time.time()

MPI.Init()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

jpa = JPA(geno, expr, start, end, mpi = True, rank = rank, comm = comm, ncore = ncore)
jpa.compute()

rr = RR(geno, expr, start, end, mpi = True, rank = rank, comm = comm, ncore = ncore, null = 'perm', sigbeta = np.repeat(0.005, geno.shape[0]))

if rank == 0:
    pvals = jpa.pvals
    scores = jpa.scores
    print (pvals.shape)
    print (scores)
    print ("Job completed in {:g} seconds".format(time.time() - start_time))



MPI.Finalize()

