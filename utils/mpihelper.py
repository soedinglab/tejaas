import numpy as np

def split_n(n, ncore):
    offset = 0
    endlist = list()
    startlist = list()
    nmax = int(n / ncore)
    for i in range(ncore-1):
        startlist.append(offset)
        endlist.append(offset + nmax)
        offset += nmax
    startlist.append(offset)
    endlist.append(n)
    return startlist, endlist

# convert this split as in split_genotype using split_n
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


def split_genotype(geno, ncore):
    start, end = split_n(geno.shape[0], ncore)
    data = [geno[i:j, :] for i, j in zip(start, end)]
    return data, start
    

def split_genemasks(masks, nsnplist, offsetlist):
    res = None
    if masks is not None:
        res = [masks[offset : nsnp + offset] for offset, nsnp in zip(offsetlist, nsnplist)]
    return res


def split_maskcomp(masks, ncore):
    #nsnps = np.array([len(x.apply2) for x in masks])
    #start, end, _ = fuzzysplit(nsnps, ncore)
    nmask = len(masks)
    start, end = split_n(nmask, ncore)
    masksplit = [masks[i:j] for i,j in zip(start, end)]
    return masksplit


def fuzzysplit(narr, ncore):
    # A fuzzy splitting algorithm, to keep similar number of snps in each core.
    # Probably useless because SVD takes more time than calculation of each SNP.
    # Kept for legacy reason and backup.
    # Use this like:
    #   start, end, _ = fuzzysplit(nsnps, ncore)
    #   masksplit = [masks[i:j] for i,j in zip(start, end)]
    # It tries to distribute the remaining number of SNPs equally in all cores.
    # When confused, it decided a boundary based on the number of SNPs in the 2 adjacent masks.

    ncum = np.cumsum(narr)
    startlist = list()
    endlist = list()
    startlist.append(0)
    ntot = narr.shape[0]
    nrem = narr.shape[0]
    crem = ncore
    h = ncum[-1] / crem
    w = 0
    wold = 0
    for i, n in enumerate(narr):
        w += n
        nrem -= 1
        #print("Element {:d} [{:d}]. Sum = {:g}, Allowed = {:g}".format(i, n, w, h))
        #print("Remaining elements: {:d}. Remaining cores: {:d}".format(nrem, crem))
        if w >= h or nrem < crem:
            #if w >= h : print("Overflow add")
            #if nrem <= crem: print ("Constraint add")
            endix = 0
            _mw = 0
            if i > 0:
                excess = w - h
                deficit = h - wold
                endix = i
                _mw = 0
                #print ("Excess: {:g}, Deficit: {:g}".format(excess, deficit), i - 1, startlist[-1])
                if excess >= deficit and (i-1) >= startlist[-1]:
                    endix = i - 1
                    _mw = n
            w = _mw
            if nrem > 0:
                #print("Adding up to index {:d} in core {:d}. Carrying over: {:g}".format(endix, len(startlist), w))
                endlist.append(endix + 1)
                startlist.append(endix + 1)
            crem -= 1
            if crem > 0:
                h = (ncum[-1] - ncum[endix]) / crem
        wold = w
        #print("wold = {:g}".format(wold))
        
    #print("Adding up to index {:d} in core {:d}.".format(i, len(startlist)))
    endlist.append(narr.shape[0])
    arrsplit = [narr[i:j] for i,j in zip(startlist, endlist)]
    return startlist, endlist, arrsplit
