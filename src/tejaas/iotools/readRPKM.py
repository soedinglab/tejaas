#!/usr/bin/env python

''' Read expression data.
    Every row contains gene expression of a donor for all genes.
    First 2 columns are FID and IID. 
    Gene expression from 3rd column onwards.
    Tab-separated.
    Reads annotation file for gene information
 
    Usage:
        from tejaas.iotools.readrpkm import ReadRPKM
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
import pandas as pd
from sklearn.decomposition import PCA

from tejaas.utils.logs import MyLogger

class ReadRPKM:


    _read_expression_once = False
    _read_annotation_once = False


    def __init__(self, filepath, dataset, npca = 0):
        self._rpkmfile = os.path.realpath(filepath)
        self._dataset = dataset
        self._genenames = None
        self._npca = npca
        self.logger = MyLogger(__name__)


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

    def _read_cardiogenics(self):
        gene_list = self._read_gene_names()
        expr_list = list()
        donor_list = list()
        with open(self._rpkmfile) as mfile:
            for line in mfile:
                linesplit = line.strip().split("\t")
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
            donor_list = mfile.readline().strip().split("\t")[1:]
            for line in mfile:
                linesplit = line.strip().split("\t")
                gene = linesplit[0].strip()
                gene_list.append(gene)
                expr = np.array([float(x) for x in linesplit[1:]])
                expr_list.append(expr)
        expr_list = np.transpose(np.array(expr_list))
        return gene_list, expr_list, donor_list

    def _read_geuvadis(self):
        gene_list = self._read_gene_names()
        expr_list = list()
        donor_list = list()
        with open(self._rpkmfile) as mfile:
            for line in mfile:
                linesplit = line.strip().split("\t")
                donor = linesplit[1].strip()
                expr = np.array([float(x) for x in linesplit[2:]])
                expr_list.append(expr)
                donor_list.append(donor)
        expr_list = np.array(expr_list)
        return expr_list, donor_list

    def _normalize_expr(self, Y):
        if isinstance(Y, pd.DataFrame):
            Y_cent = (Y.values - np.mean(Y.values, axis = 1).reshape(-1, 1)) / np.std(Y.values, axis = 1).reshape(-1, 1)
            Y_cent = pd.DataFrame(Y_cent, index=Y.index, columns=Y.columns)
            Y_cent.index.name = Y.index.name
        else:
            Y_cent = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
        return Y_cent

    def _normalize_expr_old(self, Y):
        newY = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
        return newY

    def _quant_normalize_expr(self, Y):
        from sklearn.preprocessing import normalize
        from sklearn.preprocessing import QuantileTransformer

        Y_quant = QuantileTransformer(output_distribution='normal').fit_transform(Y.T).T
        return Y_quant

    def _read_gene_expression(self):
        expr_list = list()
        donor_list = list()
        try:
            if self._dataset == "cardiogenics":
                gene_list, expr_list, donor_list = self._read_cardiogenics()
            if self._dataset == "gtex":
                gene_list, expr_list, donor_list = self._read_gtex()
            if self._dataset == "geuvadis":
                gene_list, expr_list, donor_list = self._read_geuvadis()
        except IOError as err:
            raise IOError("{:s}: {:s}".format(self._rpkmfile, err.strerror))
        expression = np.array(expr_list).transpose()
        normexpr = self._normalize_expr(expression)

        # if self._npca > 0:
        #     ## https://stats.stackexchange.com/questions/229092
        #     nComp = self._npca
        #     self.logger.debug("Using {:d} principal components".format(nComp))
        #     pca = PCA()
        #     pca.fit(expression.T)
        #     expr_pcacorr = np.dot(pca.transform(expression.T)[:, nComp:], pca.components_[nComp:,:]).T
        #     expression = expr_pcacorr
            
        self._genenames = gene_list
        self._gene_expression = normexpr
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
        return genenames
