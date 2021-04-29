#!/usr/bin/env python

''' Unfortunately there is no standard for gtf files.
    So, every version gets different function
'''

import os
import numpy as np
import gzip
from tejaas.utils.containers import GeneInfo


def gencode_v12(filepath, feature = 'gene', trim=False, biotype=['protein_coding'], include_chrom = 0, include_chroms=['{:d}'.format(x + 1) for x in range(22)]):
    annotfile = os.path.realpath(filepath)
    geneinfo = list()
    try:
        with gzip.open(annotfile, 'r') as mfile:
            for line in mfile:
                linesplit = line.decode().strip().split('\t')
                if linesplit[0][0] == '#' or linesplit[2] != feature: continue # skip header

                chrom = linesplit[0][3:]
                if include_chrom > 0:
                    include_chroms = ['{:d}'.format(include_chrom)]
                if chrom not in include_chroms: continue

                # Any particular biotype selected?
                infolist = linesplit[8].split(';')
                if len(biotype) > 0:
                    rowtype = infolist[2].strip().split(' ')[1].replace('"','')
                    if rowtype not in biotype: continue

                # TSS: gene start (0-based coordinates for BED)
                if linesplit[6] == '+':
                    start = np.int64(linesplit[3]) - 1
                    end   = np.int64(linesplit[4])
                elif linesplit[6] == '-':
                    start = np.int64(linesplit[3])  # last base of gene
                    end   = np.int64(linesplit[4]) - 1
                else:
                    raise ValueError('Strand not specified.')

                # For simulation
                if linesplit[1] == 'SIMULATION':
                    start = np.int64(linesplit[3])
                    end   = np.int64(linesplit[4])

                gene_id = infolist[0].strip().split(' ')[1].replace('"','')
                if trim:
                    gene_id = gene_id.split(".")[0]
                gene_name = infolist[4].strip().split(' ')[1].replace('"','')
                this_gene = GeneInfo(name       = gene_name,
                                     ensembl_id = gene_id,
                                     chrom      = int(chrom),
                                     start      = start,
                                     end        = end)

                geneinfo.append(this_gene)
    except IOError as err:
        raise IOError('{:s}: {:s}'.format(annotfile, err.strerror))

    return geneinfo
