#!/usr/bin/env python


import sys
import gzip
import numpy as np
import re

from tejaas.utils.containers import SnpInfo
from tejaas.utils.logs import MyLogger

class Error(Exception):
    pass

class SAMPLE_NUMBER_ERROR(Error):
    pass

class GT_FREQS_NUMBER_ERROR(Error):
    pass

class ReadOxford:

    _read_samples_once = False
    _read_genotype_once = False
    _nloci = 0
    _nsample = 0

    # deafult is GTEx:
    #     - isdosage = True

    def __init__(self, gtfile, samplefile, startsnp=0, endsnp=1e15, isdosage=True):
        self.logger = MyLogger(__name__)
        self._gtfile = gtfile
        self._samplefile = samplefile
        self._startsnp = startsnp
        self._endsnp = endsnp
        self._isdosage = isdosage
        if self._isdosage:
            self._meta_columns = 6
        else:
            self._meta_columns = 5
        self._read_genotypes()

        
    @property
    def nsample(self):
        self._read_samples()
        return self._nsample


    @property
    def samplenames(self):
        self._read_samples()
        return self._samplenames

    
    @property
    def nloci(self):
        return self._nloci


    @property
    def snpinfo(self):
        self._read_genotypes()
        return tuple(self._snpinfo)


    @property
    def dosage(self):
        return tuple(self._dosage)


    @property
    def gtnorm(self):
        return tuple(self._gtnorm)


    @property
    def gtcent(self):
        return tuple(self._gtcent)


    def _read_samples(self):
        if self._read_samples_once:
           return
        self._read_samples_once = True
        with open(self._samplefile, 'r') as samfile:
            sample = 0
            samplenames = list()
            next(samfile)
            next(samfile)
            for line in samfile:
                if re.search('^#', line):
                    continue
                sample += 1
                samplenames.append(line.strip().split()[0])
        self._nsample = sample
        self._samplenames = samplenames


    def _read_dosages(self):
        dosage = list()
        allsnps = list()
        self.logger.info("Started reading genotype.")
        self._nloci = 0
        linenum = 0
        with gzip.open(self._gtfile, 'r') as filereader:
            for snpline in filereader:
                if linenum >= self._startsnp and linenum < self._endsnp:
                    self._nloci += 1
                    mline = snpline.split()

                    if self._isdosage:
                        ngenotypes = len(mline) - self._meta_columns
                    else:
                        ngenotypes = (len(mline) - self._meta_columns) / 3

                    if float(ngenotypes).is_integer():
                        if ngenotypes != self._nsample:
                            self.logger.error('Number of samples differ from genotypes')
                            raise SAMPLE_NUMBER_ERROR;
                    else:
                        self.logger.error('Number of columns in genotype frequencies not divisible by 3')
                        raise GT_FREQS_NUMBER_ERROR;

                    if self._isdosage:
                        snp_dosage = np.array([float(x) for x in mline[self._meta_columns:]])
                    else:
                        gt_freqs = np.array([float(x) for x in mline[self._meta_columns:]])
                        indsAA = np.arange(0,self._nsample) * 3
                        indsAB = indsAA + 1
                        indsBB = indsAB + 1
                        snp_dosage = 2*gt_freqs[indsBB] + gt_freqs[indsAB] # [AA, AB, BB] := [0, 1, 2]

                    maf = sum(snp_dosage) / 2 / len(snp_dosage)
                    try:                               ######## change to get the chrom numberfrom gtfile
                       chrom = int(mline[0])
                    except:
                       chrom = -1
                    this_snp = SnpInfo(chrom      = chrom,
                                       bp_pos     = int(mline[2]),
                                       varid      = mline[1].decode("utf-8"),
                                       ref_allele = mline[3].decode("utf-8"),
                                       alt_allele = mline[4].decode("utf-8"),
                                       maf        = maf)
                    allsnps.append(this_snp)
                    dosage.append(snp_dosage)
                linenum += 1
        return allsnps, np.array(dosage)


    def _read_genotypes(self):
        if self._read_genotype_once:
            return
        self._read_genotype_once = True
        self._read_samples() # otherwise, self._nsample is not set
        allsnps, dosage = self._read_dosages()
        self.logger.info("Found {:d} SNPs of {:d} samples.".format(self._nloci, self._nsample))
        self._dosage = dosage
        self._snpinfo = allsnps
