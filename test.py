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
        

snpids  = ["snp1","snp2","snp3","snp4"]
geneids = ["gene1", "gene2", "gene3"]
pvals =  [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12]

dbwriter = DBwriter(snpids, geneids, pvals)
dbwriter.write_new_table()

db      = sqlite3.connect('testdb')
cursor  = db.cursor()
cursor.execute('''SELECT snpid, geneid, pval FROM test WHERE snpid = ?''',("snp1",))  #argument for the ? mark needs to be given in tuple
for row in cursor:
    # row['name'] returns the name column in the query, row['email'] returns email column.
    print(row)


###########Do whatever staff you want###################

cursor.execute('''DROP TABLE IF EXISTS test''')   ### To remove the table
db.close()


