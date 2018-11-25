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
    data[ncore-1] = geno[offset:nsnps, :]
    return data

def split_1darray(myarray, ncore):
    length = len(myarray)
    max_e = int(length / ncore)
    offset = 0
    data = [None for i in range(ncore)]
    for dest in range(ncore-1):
        start = offset
        end = offset + max_e
        data[dest] = myarray[start:end]
        offset += max_e 
    data[ncore-1] = myarray[offset:]
    return data
