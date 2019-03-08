import numpy as np
from utils.logs import MyLogger

logger = MyLogger(__name__)

def load(qnull_file):
    qnull = list()
    with open(qnull_file, 'r') as mfile:
        for line in mfile:
            l = line.strip().split()
            q = float(l[0].strip())
            qnull.append(q)
    return np.array(qnull)
