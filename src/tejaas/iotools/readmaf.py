import numpy as np
import pygtrie
import os
from tejaas.utils.logs import MyLogger

logger = MyLogger(__name__)

def load(snpinfo, maf_file = None):

    maf = [x.maf for x in snpinfo]

    if maf_file is not None and os.path.exists(maf_file):
        logger.debug("Updating MAF from user provided data. {:s}".format(maf_file))
        nmafx = 0
        nmafu = 0
        nallx = 0
        t = pygtrie.StringTrie()
        with open(maf_file, "r") as mfile:
            for line in mfile:
                l = line.split()
                try:
                    t[l[2]] = l[3:]
                except:
                    logger.warn("Error reading line {:s}".format(line))
                    pass
        for i, x in enumerate(snpinfo):
            try:
                info = t[x.varid]
                if info[1] == x.alt_allele:
                    nmafu += 1
                    maf[i] = float(info[2])
                elif info[1] == x.ref_allele:
                    nmafu += 1
                    maf[i] = 1 - float(info[2])
                else:
                    nallx += 1
                    logger.debug("Alleles do not match for {:s}".format(x.varid))
            except:
                nmafx += 1
                pass
        logger.debug("No MAF found for {:d} SNPs".format(nmafx))
        logger.debug("Allele mismatch for {:d} SNPs".format(nallx))
        logger.debug("MAF updated for {:d} SNPs".format(nmafu))
        del t

    maf = np.array(maf)

    return maf
