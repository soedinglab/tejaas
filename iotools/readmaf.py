import numpy as np
import pygtrie
from utils.logs import MyLogger

logger = MyLogger(__name__)

def load(snpinfo, nullmodel, maf_file):

    maf = [x.maf for x in snpinfo]

    if nullmodel == "maf":
        logger.debug("Updating MAF from 1000G data")
        nmafx = 0
        nmafu = 0
        nallx = 0
        t = pygtrie.StringTrie()
        with open(maf_file, "r") as mfile:
            for line in mfile:
                l = line.split()
                try:
                    t[l[2]] = l[3:]
                    #t[l[2]] = float(l[5])
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
