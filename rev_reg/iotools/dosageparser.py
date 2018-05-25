#!/usr/bin/env python

import numpy as np
import collections
import gzip

SNPINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'ref_allele', 'alt_allele', 'freq']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

class DosageParser:

    _read_dosage_once = False
    _read_sample_once = False


    def __init__(self, dosage_filename, sample_filename, start, end):
        self._dosage_filename = dosage_filename
        self._sample_filename = sample_filename
        self._startsnp = start
        self._endsnp = end


    @property
    def nsample(self):
        self._read_dosage()
        return self._nsample


    @property
    def nsnps(self):
        self._read_dosage()
        return self._nsnps


    @property
    def snpinfo(self):
        self._read_dosage()
        return self._snpinfo


    @property
    def dosage(self):
        self._read_dosage()
        return self._dosage


    @property
    def sample_id(self):
        self._read_sampleid()
        return self._sampleid


    def _read_dosage(self):
        if self._read_dosage_once:
            return
        self._read_dosage_once = True
        _dosagelist = list()
        _snpinfo = list()
        count_poly = 0
        count_comp = 0
        count_maf = 0
        with gzip.open(self._dosage_filename, 'r') as mfile:
            linenum = 0
            for line in mfile:
                linenum += 1
                if linenum >= self._startsnp and linenum <= self._endsnp:
                    linestrip = line.decode().strip().split()
                    if len(linestrip[3]) ==1 and len(linestrip[4]) ==1:
                        if linestrip[4] != SNP_COMPLEMENT[linestrip[3]]:
                            chrm = linestrip[0]
                            rsid = linestrip[1]
                            bpos = linestrip[2]
                            refa = linestrip[3]
                            alta = linestrip[4]
                            freq = float(linestrip[5])
                            #if (freq >= 0.10 and freq <= 0.2) or (freq >= 0.8 and freq <= 0.99):
                            if freq >= 0.10 and freq <=0.90: 
                                snp  = SnpInfo(chrm = chrm,
                                               rsid = rsid,
                                               bp_location = bpos,
                                               ref_allele = refa,
                                               alt_allele = alta,
                                               freq = freq)
                                this_dosage = np.array(linestrip[6:], dtype = float)
                                _dosagelist.append(this_dosage)
                                _snpinfo.append(snp)
                            else:
                                count_maf += 1
                        else:
                            count_comp += 1
                    else:
                        count_poly += 1
        print("\n==================================================")
        print("low maf  poly alleles  complement alleles \n")
        print(count_maf, " ", count_poly," ", count_comp,"\n")
        self._dosage = np.array(_dosagelist)
        self._nsnps = self._dosage.shape[0]
        print("Total SNPs read in : ", count_comp+count_poly+count_maf+self._nsnps)
        print("Total SNPs remained : ", self._nsnps)
        print("==================================================\n")
        self._nsample = self._dosage.shape[1]
        self._snpinfo = _snpinfo


    def _read_sampleid(self):
        if self._read_sample_once:
            return
        self._read_sample_once = True
        samplelist = list()
        with open(self._sample_filename, 'r') as mfile:
            for line in mfile:
                linestrip = line.strip().split()
                samplelist.append(linestrip[1])
        self._sampleid = samplelist
