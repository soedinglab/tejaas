import numpy as np
import itertools
import os

class Outhandler:

    def __init__(self, args, snpinfo, geneinfo):
        self.snpinfo  = snpinfo
        self.geneinfo = geneinfo
        self.args = args
        
        if not os.path.exists(os.path.dirname(self.args.outprefix)):
            os.makedirs(os.path.dirname(self.args.outprefix))

    def write_jpa_out(self, jpa):
        fname = self.args.outprefix + "_jpa.txt"
        with open(fname, "w") as f:
            f.write("snpid\tjpascore\n")
            for i, snp in enumerate(self.snpinfo):
                f.write("{:s}\t{:g}\n".format(snp.varid, jpa.scores[i]))


    def write_rr_out(self, jpa, rr):
        fname = self.args.outprefix + "_rr.txt"
        with open(fname, "w") as f:
            f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}\n".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
            for i, x in enumerate(self.snpinfo):
                f.write("{:s}\t{:d}\t{:g}\t{:g}\t{:g}\t{:g}\n".format(x.varid, x.bp_pos, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]))

        select = np.where(rr.pvals < self.args.psnpcut)[0]
        fname = self.args.outprefix + "_gene_snp_list.txt"
        with open(fname, "w") as f:
            f.write("geneid\tsnpid\tpval\n")
            for idx in select:
                gene_select = np.where(jpa.pvals[idx, :] < self.args.pgenecut)[0]
                for gidx in gene_select:
                    f.write( "{:s}\t{:s}\t{:g}\n".format(self.geneinfo[gidx].ensembl_id, self.snpinfo[idx].varid, jpa.pvals[idx][gidx]) )

        # #indices = rr.pvals < self.args.psnpcut
        # #indices = np.arange(len(rr.pvals))
        # #indices  = list(itertools.compress(range(len(indices)), indices))
        # lines = []
        # for idx in indices:
        #     selected_gene_indices = jpa.pvals[idx,:] < self.args.pgenecut
        #     selected_gene_indices = list(itertools.compress(range(len(selected_gene_indices)), selected_gene_indices))
        #     for gidx in selected_gene_indices:
        #         lines.append(str(self.geneinfo[gidx].ensembl_id) + "\t" + str(self.snpinfo[idx].varid) + "\t" + str(jpa.pvals[idx][gidx]) + "\n")
        # for line in sorted(lines, key=lambda line: line.split("\t")[0]):
        #     f.write(line)
        # f.close()
