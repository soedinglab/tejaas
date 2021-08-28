import numpy as np
from tejaas.utils.logs import MyLogger

logger = MyLogger(__name__)

def load(qnull_file):
    qnull = list()
    with open(qnull_file, 'r') as mfile:
        for line in mfile:
            l = line.strip().split()
            q = float(l[0].strip())
            qnull.append(q)
    qnull = np.array(qnull)
    qmod = qnull[np.isfinite(qnull)]
    logger.debug("Read {:d} null Q-scores".format(qmod.shape[0]))
    return qmod
