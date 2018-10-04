import numpy as np

def split_genotype(geno, ncore):
    nsnps = geno.shape[0]
    maxsnp = int(nsnps / ncore)
    offset = 0
    data = [None for i in range(ncore)]
    for dest in range(ncore-1):
        start = offset
        end = offset + maxsnp
        data[dest] = geno[start:end, :]
        offset += maxsnp
    data[ncore-1] = geno[offset:, :]
    return data
