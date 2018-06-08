#/usr/bin/env python

'''
    Distribute job to:
    a. jpa.py
    b. rr.py
    c. jpa-rr.py
    d. optim.py
'''

import logging
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

# ==================================================
# Start MPI calculation
# =================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncore = comm.Get_size()

if rank == 0:
    # this is the master
    # create a list of genotypes for sending to your slaves

    args = Args()
    out_fileprefix = args.outprefix
    genotype_filename = opts.genotype_filename
    sample_filename = opts.sample_filename
    expression_filename = opts.expression_filename
    startsnp = opts.startsnp
    endsnp = opts.endsnp

    #genotype_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/dosages/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr21.gz"
    #sample_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/dosages/donor_ids.fam"
    #expression_filename = "/usr/users/sbanerj/gwas-eQTL/gtex_data/gtex_wholeblood_normalized.expression.txt"

    start_time = time.time()

    ds = DosageParser(genotype_filename, sample_filename, startsnp, endsnp)
    dosage = ds.dosage
    snpinfo = ds.snpinfo
    donorids = ds.sample_id
    nsnps = ds.nsnps
    nsample = ds.nsample

    sampleids, expression, gene_names = read_expression(expression_filename)
    
    choose_ids = [x for x in sampleids if x in donorids]
    dosage_indices = [i for i, x in enumerate(donorids)  if x in choose_ids]
    exprsn_indices = [i for i, x in enumerate(sampleids) if x in choose_ids]
    
    geno = dosage[:, dosage_indices]
    freq = np.array([x.freq for x in snpinfo])

    print ("Completed data reading")

    # These are the inputs to cpvalcomp
    geno = norm_binom(geno, freq)
    expr = expression[:, exprsn_indices]
    #qcal = qstat.q_calibrate(G = 5000, S = 10000, p0 = 0.01)
    qcal = qstat.q_calibrate()

    print ("Completed Qstat calibration")

    maxsnp = int(nsnps / ncore)
    offset = 0
    data = [None for i in range(ncore)]
    for dest in range(ncore-1):
        start = offset
        end = offset + maxsnp
        data[dest] = geno[start:end, :]
        offset += maxsnp

    data[ncore-1] = geno[offset:nsnps, :]

else:
    data = None
    expr = None
    qcal = None

slave_geno = comm.scatter(data, root = 0)
expr = comm.bcast(expr, root = 0)
qcal = comm.bcast(qcal, root = 0)
comm.barrier()

# ==================================
# Everything sent. Now do the calculations
# ==================================
pvals, qscores, pqvals, gene_indices = cpvalcomp(slave_geno, expr, qcal, shuffle=False)

pvals   = comm.gather(pvals,   root = 0)
qscores = comm.gather(qscores, root = 0)
pqvals  = comm.gather(pqvals,  root = 0)
gene_indices = comm.gather(gene_indices,  root = 0)

if rank == 0:
    pvals = np.vstack(pvals)
    #np.save(out_fileprefix, pvals)
    qscores = np.concatenate(qscores)
    pqvals  = np.concatenate(pqvals)
    gene_indices_list = list()
    for x in gene_indices:
        gene_indices_list += x
    output.write(out_fileprefix, snpinfo, gene_names, pvals, qscores, pqvals, gene_indices_list)
    print ("Job completed in {:g} seconds".format(time.time() - start_time))
else:
    assert gene_indices is None
    assert qscores is None
    assert pvals   is None
    assert pqvals  is None

print ("Done")
