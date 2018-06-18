import numpy as np
import itertools

class Outhandler:
    def __init__(self, args, snpinfo, geneinfo):
        self.snpinfo  = snpinfo
        self.geneinfo = geneinfo
    
    def write_jpa_out(jpa):
        f = open(args.outprefix + "_jpa.txt", "w")
        f.write("snpid\tjpascore")
        outlist = ["\n" + self.snpinfo[i].varid + "\t" + str(jpa.scores[i])for i in range(len(self.sninfo))]
        f.writelines(outlist)
        f.close()
    def write_rr_out(jpa, rr):
        f = open(args.outprefix + "_rr.txt", "w")
        f.write("{:s}\t{:s}\t{:s}\t{:s}\t{:s}\t{:s}".format('ID', 'Pos', 'Q', 'Mu', 'Sigma', 'P'))
        outlist = ["\n{:s}\t {:d}\t{:g}\t{:g}\t{:g}\t{:g}".format(x.varid, x.bp_pos, rr.scores[i], rr.null_mu[i], rr.null_sigma[i], rr.pvals[i]) for i, x in enumerate(self.snpinfo)]
        f.writelines(outlist)
        f.close()

        indices = rr.pvals < args.psnpcut
        indices  = list(itertools.compress(range(len(indices)), indices))
        #selected_snps = [snpinfo[i] for i in indices]
        f = open(args.outprefix + "_gene_snp_list.txt", "w")
        f.write("geneid\tsnpid\tpval\n")
        lines = []
        for idx in indices:
            selected_gene_indices = jpa.pvals[idx,:] < args.pgenecut
            selected_gene_indices = list(itertools.compress(range(len(selected_gene_indices)), selected_gene_indices))
            for gidx in selected_gene_indices:
                lines.append(str(self.geneinfo[gidx].ensembl_id) + "\t" + str(self.snpinfo[idx].varid) + "\t" + str(jpa.pvals[idx][gidx]) + "\n")
        for line in sorted(lines, key=lambda line: line.split("\t")[0]):
            f.write(line)
        f.close()


