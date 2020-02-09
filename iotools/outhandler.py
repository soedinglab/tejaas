import numpy as np
import itertools
import os

from utils.logs import MyLogger
logger = MyLogger(__name__)

import mpmath
mpmath.mp.dps = 500
def pval(x): return float(1 - 0.5 * (1 + mpmath.erf(x/mpmath.sqrt(2))))

class Outhandler:

    def __init__(self, args, snpinfo, geneinfo):
        self.snpinfo  = snpinfo
        self.geneinfo = geneinfo
        self.args = args
        dirname = os.path.dirname(os.path.realpath(self.args.outprefix))
        logger.debug('Writing result in: {:s}'.format(dirname))
        if not os.path.exists(dirname): os.makedirs(dirname)


    def write_jpa_out(self, jpa):
        fname = self.args.outprefix + "_jpa.txt"
        scores = jpa.scores
        pvals = jpa.jpa_pvals
        with open(fname, "w") as f:
            f.write("snpid\tjpascore\tp-value\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s}\t{:g}\t{:g}\n".format(snp.varid, scores[i], pvals[i]))


    def write_rr_out(self, rr, jpa, snp_select_idx, suffix = "", selected_snps = [], selected_genes = [], write_betas = False):
        mysnpinfo = self.snpinfo
        if len(selected_snps):
            mysnpinfo = [self.snpinfo[int(i)] for i in selected_snps]
        fname = self.args.outprefix + "_rr" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'CHR', 'Pos', 'MAF', 'Q', 'Mu', 'Sigma', 'P'))
            for i, x in enumerate(mysnpinfo):
                if rr.pvals[i] == 0:
                    rr.pvals[i] = pval( (rr.scores[i] - rr.null_mu[i]) / rr.null_sigma[i])
                f.write("{:s}\t{:d}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.chrom, x.bp_pos, x.maf, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]))

        if rr.betas is not None and write_betas:
            betafile = self.args.outprefix + "_betas" + suffix + ".txt"
            with open(betafile, 'w') as outstream:
                gene_names = " ".join([g.ensembl_id for g in self.geneinfo])
                outstream.write(gene_names+"\n")
                np.savetxt(outstream, rr.betas, fmt='%1.4e')

        if len(selected_genes):
            np.savetxt(self.args.outprefix + "_selected_genes" + suffix + ".txt", selected_genes, fmt='%i')
            pvals = jpa.pvals[selected_snps]

        fname = self.args.outprefix + "_gene_snp_list" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for i, sidx in enumerate(snp_select_idx):
                gene_select_idx = np.where(jpa.pvals[i, :] < self.args.pgenecut)[0]
                for gidx in gene_select_idx:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[sidx].varid, jpa.pvals[i, gidx]) )
