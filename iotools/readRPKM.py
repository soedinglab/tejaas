#!/usr/bin/env python

''' Read expression data.
    Every row contains gene expression of a donor for all genes.
    First 2 columns are FID and IID. 
    Gene expression from 3rd column onwards.
    Tab-separated.
    ** Tries to find a .geneids file in the same folder to get the gene name
    Reads annotation file for gene information
 
    Usage:
        from iotools.readrpkm import ReadRPKM
        rpkm = ReadRPKM(rpkmfile)
        rpkm.read_annotation(annotationfile)
        expression = rpkm.expression
           ...

   Returns:
        a) expression matrix
        b) identity of the donors
        c) gene names
'''

import numpy as np
import os

class ReadRPKM:


    _read_expression_once = False
    _read_annotation_once = False


    def __init__(self, filepath, dataset):
        self._rpkmfile = os.path.realpath(filepath)
        self._dataset = dataset
        self._genenames = None


    @property
    def expression(self):
        self._run_once()
        return self._gene_expression


    @property
    def donor_ids(self):
        self._run_once()
        return self._donor_ids


    @property
    def gene_names(self):
        self._run_once()
        return self._genenames


    def _run_once(self):
        if self._read_expression_once:
            return
        self._read_expression_once = True
        self._read_gene_expression()
        if self._dataset == "cardiogenics":
            self._read_gene_names()
        if self._dataset == "geuvadis":
            self._read_gene_names()

    def _read_cardiogenics(self):
        expr_list = list()
        donor_list = list()
        with open(self._rpkmfile) as mfile:
            for line in mfile:
                linesplit = line.strip().split()
                donor = linesplit[0].strip()
                expr = np.array([float(x) for x in linesplit[2:]])
                expr_list.append(expr)
                donor_list.append(donor)
        return expr_list, donor_list
    
    def _read_gtex(self):
        expr_list = list()
        donor_list = list()
        gene_list = list()
        with open(self._rpkmfile) as mfile:
            donor_list = mfile.readline().strip().split()[1:]
            for line in mfile:
                linesplit = line.strip().split()
                gene = linesplit[0].strip()
                gene_list.append(gene)
                expr = np.array([float(x) for x in linesplit[1:]])
                expr_list.append(expr)
        expr_list = np.transpose(np.array(expr_list))
        self._genenames = gene_list
        return expr_list, donor_list

    def _read_geuvadis(self):
        expr_list = list()
        donor_list = list()
        with open(self._rpkmfile) as mfile:
            for line in mfile:
                linesplit = line.strip().split()
                donor = linesplit[1].strip()
                expr = np.array([float(x) for x in linesplit[2:]])
                expr_list.append(expr)
                donor_list.append(donor)
        expr_list = np.array(expr_list)
        return expr_list, donor_list


    def _read_gene_expression(self):
        expr_list = list()
        donor_list = list()
        try:
            if self._dataset == "cardiogenics":
                expr_list, donor_list = self._read_cardiogenics()
            if self._dataset == "gtex":
                expr_list, donor_list = self._read_gtex()
            if self._dataset == "geuvadis":
                expr_list, donor_list = self._read_geuvadis()
        except IOError as err:
            raise IOError("{:s}: {:s}".format(self._rpkmfile, err.strerror))
        expression = np.array(expr_list).transpose()
        self._gene_expression = expression
        self._donor_ids = donor_list


    def _read_gene_names(self):
        genefile = os.path.splitext(self._rpkmfile)[0] + ".geneids"
        genenames = list()
        try:
            with open(genefile) as mfile:
                for line in mfile:
                    gene = line.strip().split("\t")[0]
                    genenames.append(gene)
        except IOError as err:
            raise IOError("{:s}: {:s}".format(genefile, err.strerror))
        self._genenames = genenames
