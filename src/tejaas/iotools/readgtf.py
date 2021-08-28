#!/usr/bin/env python

''' Unfortunately there is no standard for gtf files.
    So, every version gets different function
'''

import os
import numpy as np
import gzip
import re
from tejaas.utils.containers import GeneInfo

def gencode(filepath, feature = 'gene', trim=False, biotype=['protein_coding', 'lncRNA'], include_chrom = 0, include_chroms=['{:d}'.format(x + 1) for x in range(22)]):
    annotfile = os.path.realpath(filepath)
    geneinfo = list()
    lncRNA_list = ["macro_lncRNA", "non_coding", "bidirectional_promoter_lncRNA", "3prime_overlapping_ncRNA", 
                   "sense_overlapping", "processed_transcript", "sense_intronic", "antisense", "lincRNA", "miRNA"]

    if "lncRNA" in biotype:
        mybiotype = biotype + lncRNA_list
    else:
        mybiotype = biotype
    if annotfile.endswith("affy"):
        geneinfo = affy_exon_chip(annotfile, include_chrom = include_chrom, include_chroms=include_chroms)
        return geneinfo
    else:
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
                    indices_set = False
                    gene_name_ix = None
                    gene_type_ix = None
                    gene_id_ix   = None
                    if not indices_set:
                        for i,e in enumerate(infolist):
                            Id    = e.strip().split(' ')[0]
                            # value = e.strip().split(' ')[1]
                            if Id == "gene_type": gene_type_ix = i
                            if Id == "gene_name": gene_name_ix = i
                            if Id == "gene_id": gene_id_ix = i
                        if (gene_type_ix is None) and (len(mybiotype) > 0):
                            raise ValueError("No 'gene_type' specified in gene annotation file")
                        if (gene_name_ix is not None) and (gene_id_ix is not None):
                            indices_set = True

                    if len(mybiotype) > 0:
                        rowtype = infolist[gene_type_ix].strip().split(' ')[1].replace('"','')
                        if rowtype not in mybiotype: continue
                    gene_name = infolist[gene_name_ix].strip().split(' ')[1].replace('"','')


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
                    this_gene = GeneInfo(name       = gene_name,
                                         ensembl_id = gene_id,
                                         chrom      = int(chrom),
                                         start      = start,
                                         end        = end)

                    geneinfo.append(this_gene)
        except IOError as err:
            raise IOError('{:s}: {:s}'.format(annotfile, err.strerror))

        ensembl_ids = [g.ensembl_id for g in geneinfo]
        try:
            assert( len(ensembl_ids) == len(set(ensembl_ids)) )
        except AssertionError:
            print('GTF annotation file contains non-unique gene_id identifiers')
        return geneinfo

def affy_exon_chip(filepath, include_chrom = 0, include_chroms=['{:d}'.format(x + 1) for x in range(22)]):
    geneinfo = list()
    try:
        with open(filepath, 'r') as mfile:
            next(mfile) # skip header
            for line in mfile:
                linesplit = line.strip().split('\t')
                if linesplit[0][0] == '#' : continue 

                chrom = linesplit[2][3:]
                if include_chrom > 0:
                    include_chroms = ['{:d}'.format(include_chrom)]
                if chrom not in include_chroms: continue

                rowtype = None
                gene_id = None
                gene_name = None
                infolist = [x.strip() for x in linesplit[8].split(';')]
                for info in infolist:
                    if info.startswith('gene_type'):
                        rowtype = info.split(' ')[1].replace('"','')
                    elif info.startswith('gene_id'):
                        gene_id = info.split(' ')[1].replace('"','')
                    elif info.startswith('gene_name'):
                        gene_name = info.split(' ')[1].replace('"','')

                # Any particular biotype selected?
                if 'all' not in biotype:
                    if rowtype not in biotype: continue

                # TSS: gene start (0-based coordinates for BED)
                if linesplit[3] == '+':
                    start = np.int64(linesplit[4]) - 1
                    end   = np.int64(linesplit[5])
                elif linesplit[3] == '-':
                    start = np.int64(linesplit[4])  # last base of gene
                    end   = np.int64(linesplit[5]) - 1
                else:
                    raise ValueError('Strand not specified.')

                gene_name = linesplit[7].split("//")[0].rstrip()
                transcript_cluster_id = linesplit[0]
                this_gene = GeneInfo(name       = gene_name,
                                     ensembl_id = transcript_cluster_id,
                                     chrom      = int(chrom),
                                     start      = start,
                                     end        = end)

                geneinfo.append(this_gene)
    except IOError as err:
        raise IOError('{:s}: {:s}'.format(annotfile, err.strerror))

    return geneinfo
