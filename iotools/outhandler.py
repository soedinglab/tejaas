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
    def write_jpa_out(self,jpa):
        f = open(self.args.outprefix + "_jpa.txt", "w")
        f.write("snpid\tjpascore")
        outlist = ["\n" + self.snpinfo[i].varid + "\t" + str(jpa.scores[i])for i in range(len(self.snpinfo))]
        f.writelines(outlist)
        f.close()
    def write_rr_out(self,jpa, rr):
        f = open(self.args.outprefix + "_rr.txt", "w")
        f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
        outlist = ["\n{:s}\t {:d}\t{:g}\t{:g}\t{:g}\t{:g}".format(x.varid, x.bp_pos, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]) for i, x in enumerate(self.snpinfo)]
        f.writelines(outlist)
        f.close()

        indices = rr.pvals < self.args.psnpcut
        indices  = list(itertools.compress(range(len(indices)), indices))
        #selected_snps = [snpinfo[i] for i in indices]
        f = open(self.args.outprefix + "_gene_snp_list.txt", "w")
        f.write("geneid\tsnpid\tpval\n")
        lines = []
        for idx in indices:
            selected_gene_indices = jpa.pvals[idx,:] < self.args.pgenecut
            selected_gene_indices = list(itertools.compress(range(len(selected_gene_indices)), selected_gene_indices))
            for gidx in selected_gene_indices:
                lines.append(str(self.geneinfo[gidx].ensembl_id) + "\t" + str(self.snpinfo[idx].varid) + "\t" + str(jpa.pvals[idx][gidx]) + "\n")
        for line in sorted(lines, key=lambda line: line.split("\t")[0]):
            f.write(line)
        f.close()


