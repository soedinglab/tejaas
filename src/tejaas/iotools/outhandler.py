import numpy as np
import itertools
import os

from tejaas.utils.logs import MyLogger
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


    def write_jpa_out(self, jpa, tgjpa, snp_select_idx):
        fname = self.args.outprefix + "_jpa.txt"
        scores = jpa.scores
        pvals = jpa.jpa_pvals
        with open(fname, "w") as f:
            f.write("snpid\tjpascore\tp-value\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s}\t{:g}\t{:g}\n".format(snp.varid, scores[i], pvals[i]))
        
        fname = self.args.outprefix + "_gene_snp_list.txt"
        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for i, sidx in enumerate(snp_select_idx):
                gene_select_idx = np.where(tgjpa.pvals[i, :] < self.args.pgenecut)[0]
                for gidx in gene_select_idx:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[sidx].varid, tgjpa.pvals[i, gidx]) )

    def write_rr_out(self, rr, tgknn, jpa, snp_select_idx, suffix = "", write_betas = False):
        mysnpinfo = self.snpinfo
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
                np.savetxt(outstream, rr.betas[snp_select_idx,:], fmt='%1.4e')

        fname = self.args.outprefix + "_gene_snp_list_knn" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for i, sidx in enumerate(snp_select_idx):
                gene_select_idx = np.where(tgknn.pvals[i, :] < self.args.pgenecut)[0]
                for gidx in gene_select_idx:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[sidx].varid, tgknn.pvals[i, gidx]) )

        fname = self.args.outprefix + "_gene_snp_list" + suffix + ".txt"
        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for i, sidx in enumerate(snp_select_idx):
                gene_select_idx = np.where(jpa.pvals[i, :] < self.args.pgenecut)[0]
                for gidx in gene_select_idx:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[sidx].varid, jpa.pvals[i, gidx]) )

        if jpa.target_fdr is not None:
            fname = self.args.outprefix + "_gene_snp_list_FDR" + suffix + ".txt"
            with open(fname, 'w') as f:
                f.write("geneid\tsnpid\tpval\tadj_pval\n")
                for i, pairs in enumerate(jpa.pass_fdr):
                    f.write( "{:s}\t{:s}\t{:g}\t{:g}\n".format(self.geneinfo[pairs[1]].ensembl_id, self.snpinfo[snp_select_idx[pairs[0]]].varid, pairs[2], jpa.adj_pvals[i]) )


