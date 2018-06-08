#!/usr/bin/env python

def write(filepath, snpinfo, gene_names, pvals, qscore, p_qscore, gene_indices):
    outfile = '{:s}.qscores'.format(filepath)
    with open(outfile, 'w') as mfile:
        for i, snp in enumerate(snpinfo):
            thisline = '{:25s}\t{:g}\t{:g}\n'.format(snp.rsid, qscore[i], p_qscore[i])
            mfile.write(thisline)

    outfile = '{:s}.topgenes'.format(filepath)
    with open(outfile, 'w') as mfile:
        for i, snp in enumerate(snpinfo):
            if len(gene_indices[i]) > 0:
                for j in gene_indices[i]:
                    thisline = '{:25s}\t{:25s}\t{:g}\n'.format(snp.rsid, gene_names[j], pvals[i, j])
                    mfile.write(thisline)
