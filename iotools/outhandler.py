import numpy as np
import itertools
import os

from utils.logs import MyLogger
logger = MyLogger(__name__)

class Outhandler:

    def __init__(self, args, snpinfo, geneinfo):
        self.snpinfo  = snpinfo
        self.geneinfo = geneinfo
        self.args = args
        dirname = os.path.dirname(os.path.realpath(self.args.outprefix))
        logger.debug('Writing result in: {:s}'.format(dirname))
        if not os.path.exists(dirname): os.makedirs(dirname)

    def write_fstat_out(self, cpma):
        fname = self.args.outprefix + "_cpma_fstat.txt"
        fstats = cpma.fstats
        gene_names = " ".join([g.ensembl_id for g in self.geneinfo])
        with open(fname, "w") as f:
            f.write("snpid "+gene_names+"\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s} {:s}\n".format(snp.varid, " ".join(["{:g}".format(i) for i in fstats[i,:]])))

    def write_jpa_out(self, jpa):
        fname = self.args.outprefix + "_jpa.txt"
        scores = jpa.scores
        with open(fname, "w") as f:
            f.write("snpid\tjpascore\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s}\t{:g}\n".format(snp.varid, scores[i]))

    def write_jpa_pvals(self, jpa):
        fname = self.args.outprefix + "_jpa_pvals.txt"
        scores = jpa.scores
        pvals = jpa.pcpma
        with open(fname, "w") as f:
            f.write("snpid\tjpascore\tp-value\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s}\t{:g}\t{:g}\n".format(snp.varid, scores[i], pvals[i]))


    def write_rr_out(self, jpa, rr, prefix = "", selected_snps = [], selected_genes = [], write_betas = False):
        mysnpinfo = self.snpinfo
        if len(selected_snps):
            mysnpinfo = [self.snpinfo[int(i)] for i in selected_snps]
        fname = self.args.outprefix + "_rr" + prefix + ".txt"
        with open(fname, "w") as f:
            f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
            for i, x in enumerate(mysnpinfo):
                f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]))
        if rr.betas is not None and write_betas:
            np.savetxt(self.args.outprefix + "_betas" + prefix + ".txt", rr.betas, fmt='%1.4e')
        select = np.where(rr.pvals < self.args.psnpcut)[0]
        fname = self.args.outprefix + "_gene_snp_list" + prefix + ".txt"
        pvals = jpa.pvals

        if len(selected_genes):
            np.savetxt(self.args.outprefix + "_selected_genes" + prefix + ".txt", selected_genes, fmt='%i')
            pvals = jpa.pvals[selected_snps]

        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for idx in select:
                gene_select = np.where(pvals[idx, :] < self.args.pgenecut)[0]
                for gidx in gene_select:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[idx].varid, pvals[idx][gidx]) )

    # def write_rr_out(self, jpa, rr):
    #     fname = self.args.outprefix + "_rr.txt"
    #     with open(fname, "w") as f:
    #         f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
    #         for i, x in enumerate(self.snpinfo):
    #             f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]))

    #     select = np.where(rr.pvals < self.args.psnpcut)[0]
    #     fname = self.args.outprefix + "_gene_snp_list.txt"
    #     pvals = jpa.pvals
    #     if self.selected is not None and len(self.selected):
    #         pvals = jpa.pvals[self.selected]
    #     with open(fname, "w") as f:
    #         f.write("geneid\tsnpid\tpval\n")
    #         for idx in select:
    #             gene_select = np.where(pvals[idx, :] < self.args.pgenecut)[0]
    #             for gidx in gene_select:
    #                 f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[idx].varid, pvals[idx][gidx]) )
