import sqlite3
import itertools

class DBwriter:
    def __init__(self, snpinfo, geneinfo, pvals,dbname, tablename):
        self.geneids  = [x.ensembl_id for x in geneinfo]                                    #G
        self.snpids   = [x.varid for x in snpinfo]            #I
        self.pvals  = pvals                                   #I X G
        self.dbname = dbname
        self.tablename = tablename
    def write_new_table(self):
        try:
            db      = sqlite3.connect(dbname)
            cursor  = db.cursor()
            cursor.execute(''' CREATE TABLE IF NOT EXISTS ? (id INTEGER PRIMARY KEY, snpid TEXT, geneid TEXT, pval REAL)''',(self.tablename,))
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
    def __init__(self, snpinfo, pvals, rscores, mu, sigma, outfile):
        self.snpids     = [x.varid for x in snpinfo]
        self.pvals      = pvals
        self.rscores    = rscores
        self.mu         = mu
        self.sigma      = sigma
        self.outfile    = outfile
    def write(self):
        f = open(self.outfile, "w")
        f.write("snpid\trscore\tpval\tmu\tsigma")
        outlist = ["\n" + str(self.snpids[i]) + "\t" + str(self.rscores[i]) + "\t" + str(self.pvals[i]) + "\t" + str(self.mu[i]) + "\t" + str(self.sigma[i]) for i in range(len(self.snpids))]
        f.writelines(outlist)
        f.close()

class jpaOutWriter:
    def __init__(self, snpinfo, jpascores, outfile):
        self.snpids     = [x.varid for x in snpinfo]
        self.jpascores    = jpascores
        self.outfile    = outfile
    def write(self):
        f = open(self.outfile, "w")
        f.write("snpid\tjpascore")
        outlist = ["\n" + self.snpids[i] + "\t" + str(self.jpascores[i])for i in range(len(self.snpids))]
        f.writelines(outlist)
        f.close()

