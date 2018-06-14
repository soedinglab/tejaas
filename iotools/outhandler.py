import numpy as np
import sqlite3
import itertools

class DBwriter:
    def __init__(self, snpids, geneids, pvals):
        self.geneids  = geneids                                    #G
        self.snpids   = snpids                                     #I
        self.pvals     = pvals                                   #I X G
    
    def write_new_table(self):
        try:
            db      = sqlite3.connect('testdb')
            cursor  = db.cursor() 
            cursor.execute(''' CREATE TABLE IF NOT EXISTS test (id INTEGER PRIMARY KEY, snpid TEXT, geneid TEXT, pval REAL)''')
            ixg     = list(itertools.product(self.snpids, self.geneids))
            tuples  = [ixg[i] + (pvals[i],) for i in range(len(pvals))]
            cursor.executemany(''' INSERT INTO test(snpid, geneid, pval) VALUES(?,?,?)''', tuples)
            db.commit()
        except Exception as e:
            db.rollback()
            print(e)
        finally:
            db.close()
        


class rrOutWriter:
    def __init__(self, snpids, pvals, rscores, mu, sigma):
        self.snpids     = snpids
        self.pvals      = pvals
        self.rscores    = rscores
        self.mu         = mu
        self.sigma      = sigma
        self.outfile    = "rr_out.txt"
    def writetxt(self):
        f = open(self.outfile, "w")
        f.write("snpid\trscore\tpval\tmu\tsigma")
        outlist = ["\n" + snpids[i] + "\t" + str(rscores[i]) + "\t" + str(pvals[i]) + "\t" + str(mu[i]) + "\t" + str(sigma[i]) for i in range(len(snpids))]
        f.writelines(outlist)
        f.close()
            
snpids = ["snp1","snp2","snp3"]
rscores = [12, 13, 14]
pvals = [0.1,0.01, 0.001]
mu = [0,0,0]
sigma = [1,1,1]
rrwriter = rrOutWriter(snpids, pvals, rscores, mu, sigma)
rrwriter.writetxt()

