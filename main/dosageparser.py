#!/usr/bin/env python

import numpy as np
import collections
import gzip

SNPINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'ref_allele', 'alt_allele', 'freq']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()


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
        with gzip.open(self._dosage_filename, 'r') as mfile:
            linenum = 0
            for line in mfile:
                linenum += 1
                if linenum >= self._startsnp and linenum <= self._endsnp:
                    linestrip = line.decode().strip().split()
                    chrm = linestrip[0]
                    rsid = linestrip[1]
                    bpos = linestrip[2]
                    refa = linestrip[3]
                    alta = linestrip[4]
                    freq = float(linestrip[5])
                    snp  = SnpInfo(chrm = chrm,
                                   rsid = rsid,
                                   bp_location = bpos,
                                   ref_allele = refa,
                                   alt_allele = alta,
                                   freq = freq)
                    this_dosage = np.array(linestrip[6:], dtype = float)
                    _dosagelist.append(this_dosage)
                    _snpinfo.append(snp)
        self._dosage = np.array(_dosagelist)
        self._nsnps = self._dosage.shape[0]
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
