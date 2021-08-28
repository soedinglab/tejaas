#/usr/bin/env python

''' Parse dosage data from VCF file.
    Usage:
       from tejaas.iotools.readvcf import ReadVCF
       readvcf = ReadVCF(vcf_filepath)
       snpinfo = readvcf.snpinfo
          ...
    Returns:
       a) information of all the variants read
            - chrom
            - pos
            - id
            - reference allele
            - alternate allele
            - minor allele frequency
       b) genotype matrix
       c) identity of the donors
'''

import numpy as np
import gzip
import os
import re
from tejaas.utils.containers import SnpInfo

class ReadVCF:


    _read_genotype_once = False


    def __init__(self, filepath, startsnp, endsnp, mode="DS", samplefile=None):
        self._filepath = filepath
        self._startsnp = startsnp
        self._endsnp = endsnp
        self._mode = mode
        if samplefile is not None:
            self._samplefile = samplefile
        else:
            self._samplefile = ""

    @property
    def snpinfo(self):
        self._run_once()
        return tuple(self._snpinfo)


    @property
    def dosage(self):
        self._run_once()
        return self._dosage


    @property
    def donor_ids(self):
        self._run_once()
        return self._donor_ids


    def _run_once(self):
        if self._read_genotype_once:
            return
        self._read_genotype_once = True
        self._read_samples()
        self._read_dosage()

    def _read_samples(self):
        if os.path.exists(self._samplefile):
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

    def _read_dosage(self):
        dosage = list()
        snpinfo = list()
        linenum = 0
        with gzip.open(self._filepath, 'r') as vcf:
            for line in vcf:
                linestrip = line.decode().strip()
                if linestrip[:2] == '##': continue
                if linestrip[:6] == '#CHROM':
                    linesplit = linestrip.split("\t")
                    donor_ids = linesplit[9:]
                else:
                    if linenum >= self._startsnp and linenum < self._endsnp:
                        linesplit = linestrip.split("\t")
                        if linesplit[0].startswith("chr"):
                            chrom = int(linesplit[0][3:])
                        else:
                            chrom = int(linesplit[0])
                        pos   = int(linesplit[1])
                        varid = linesplit[2]
                        ref   = linesplit[3]
                        alt   = linesplit[4]

                        if self._mode == "DS":
                            if "DS" not in linesplit[8].split(':'):
                                self._mode = "GT"
                            else:
                                dsindx = linesplit[8].split(':').index("DS")
                                ds = [x.split(':')[dsindx] if len(x) > 1 else "." for x in linesplit[9:]]
                                gtindx = linesplit[8].split(':').index("GT")
                                for i, x in enumerate(ds):
                                    if x == ".":
                                        gt = linesplit[9+i].split(':')[gtindx]
                                        if len(gt) == 3 and gt[0] != "." and gt[2] != ".":
                                            ds[i] = float(int(gt[0]) + int(gt[2]))

                        if self._mode == "GT":
                            if "GT" not in linesplit[8].split(':'):
                                print("ERROR: no GT field in VCF file")
                                raise
                            gtindx = linesplit[8].split(':').index("GT")
                            gt = [x.split(':')[gtindx] if len(x) > 1 else "." for x in linesplit[9:]]
                            ds = [ float(int(x[0]) + int(x[2])) if len(x) == 3 and x[0] != "." and x[2] != "." else "." for x in gt ]

                        ds_notna = [float(x) for x in ds if x != "."]
                        freq = sum(ds_notna) / 2 / len(ds_notna)
                        maf = freq
                        snpdosage = [float(x) if x != '.' else 2 * freq for x in ds]
                        #if freq > 0.5:
                        #    maf = 1 - freq
                        #    ref = linesplit[4]
                        #    alt = linesplit[3]
                        #    snpdosage = [2 - x for x in snpdosage]

                        this_snp = SnpInfo(chrom      = chrom,
                                           bp_pos     = pos,
                                           varid      = varid,
                                           ref_allele = ref,
                                           alt_allele = alt,
                                           maf        = maf)

                        dosage.append(snpdosage)
                        snpinfo.append(this_snp)
                    if linenum > self._endsnp: break
                    linenum += 1

        if os.path.exists(self._samplefile):
            common_ids = [x for x in self._samplenames if x in donor_ids]
            print("GT Sample selection: {:d} samples were retained from a total of {:d} samples".format(len(common_ids), len(donor_ids)))
            ix = [donor_ids.index(x) for x in common_ids]
            self._dosage = np.array(dosage)[:,ix]
            self._snpinfo = snpinfo
            self._donor_ids = common_ids
        else:
            self._dosage = np.array(dosage)
            self._snpinfo = snpinfo
            self._donor_ids = donor_ids
