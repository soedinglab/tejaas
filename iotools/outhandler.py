import sqlite3
import itertools
import numpy as np

class DBwriter:
    def __init__(self, snpinfo, geneinfo, pvals,dbname):
        self.snpinfo  = snpinfo
        self.geneinfo = geneinfo
        self.geneids  = [x.ensembl_id for x in geneinfo]                                    #G
        self.snpids   = [x.varid for x in snpinfo]            #I
        self.pvals  = pvals                                   #I X G
        self.dbname = dbname
    def write(self):
        try:
            db      = sqlite3.connect(self.dbname)
            cursor  = db.cursor()
            cursor.execute(''' CREATE TABLE IF NOT EXISTS pvals (id INTEGER PRIMARY KEY, snpid TEXT, geneid TEXT, pval REAL)''')
            ixg     = list(itertools.product(self.snpids, self.geneids))
            tuples  = [ixg[i*(self.pvals.shape[1])+ j] + (self.pvals[i][j],) for i in range(self.pvals.shape[0]) for j in range(self.pvals.shape[1])]
            cursor.executemany(''' INSERT INTO pvals (snpid, geneid, pval) VALUES(?,?,?)''', tuples)
            
            cursor.execute(''' CREATE TABLE IF NOT EXISTS snpinfo (id INTEGER PRIMARY KEY, varid TEXT, chrom INTEGER, bp_pos INTEGER, ref_allele TEXT, alt_allele TEXT, maf REAL)''')
            snp_tuples = list((str(i.varid), i.chrom, i.bp_pos, i.ref_allele, i.alt_allele, i.maf) for i in self.snpinfo)
            cursor.executemany(''' INSERT INTO snpinfo (varid, chrom, bp_pos, ref_allele, alt_allele, maf) VALUES(?,?,?,?,?,?)''', snp_tuples)

            cursor.execute(''' CREATE TABLE IF NOT EXISTS geneinfo (id INTEGER PRIMARY KEY, ensembl_id TEXT, name TEXT, chrom INTEGER, start INTEGER, end INTEGER)''')
            gene_tuples = list((i.ensembl_id, i.name, i.chrom, i.start, i.end) for i in self.geneinfo)
            cursor.executemany(''' INSERT INTO geneinfo (ensembl_id, name, chrom, start, end) VALUES(?,?,?,?,?)''', gene_tuples)

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

