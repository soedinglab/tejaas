import numpy as np

def split_genotype(geno, ncore):
    nsnps = geno.shape[0]
    maxsnp = int(nsnps / ncore)
    offset = 0
    offsetlist = [offset]
    data = [None for i in range(ncore)]
    for dest in range(ncore-1):
        start = offset
        end = offset + maxsnp
        data[dest] = geno[start:end, :]
        offset += maxsnp
        offsetlist.append(offset)
    data[ncore-1] = geno[offset:, :]
    return data, offsetlist


def split_genemasks(masks, nsnplist, offsetlist):
    res = None
    if masks is not None:
        res = [masks[offset : nsnp + offset] for offset, nsnp in zip(offsetlist, nsnplist)]
    return res
