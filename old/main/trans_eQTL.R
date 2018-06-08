library(argparser)
library(Rmpi)

p <- arg_parser("Get Values from job_submit.sh")
p <- add_argument(p, "GENOTYPE", help="")
p <- add_argument(p, "GENEEXP", help="")
p <- add_argument(p, "DONORIDS", help="")
p <- add_argument(p, "STARTSNP", help="")
p <- add_argument(p, "ENDSNP", help="")
p <- add_argument(p, "P_QVAL_CUTOFF", help="")
p <- add_argument(p, "PVAL_CUTOFF", help="")
p <- add_argument(p, "OUTPUT_QSTAT", help="")
p <- add_argument(p, "OUTPUT_TOPGENES", help="")
argv <- parse_args(p)

TRANSEQTL <- argv$TRANSEQTL
GENOTYPE <- argv$GENOTYPE
GENEEXP <- argv$GENEEXP
STARTSNP <- argv$STARTSNP
ENDSNP <- argv$ENDSNP
DONORIDS <- argv$DONORIDS
P_QVAL_CUTOFF <- argv$P_QVAL_CUTOFF
PVAL_CUTOFF <- argv$PVAL_CUTOFF
OUTPUT_QSTAT <- argv$OUTPUT_QSTAT
OUTPUT_TOPGENES <- argv$OUTPUT_TOPGENES

startsnp = as.numeric(STARTSNP)
endsnp = as.numeric(ENDSNP)
p_qval_cutoff = as.numeric(P_QVAL_CUTOFF) 
pval_cutoff = as.numeric(PVAL_CUTOFF)
expression = read.table(GENEEXP, header = T, row.names = 1)
donor_ids = read.table(DONORIDS)[,2]
donor_ids = gsub('-', '.', donor_ids)
genotype = read.table(GENOTYPE, header = F, row.names = 2, skip = startsnp, nrows = (endsnp - startsnp))[, 6: (length(donor_ids) + 5)]
colnames(genotype) = donor_ids
output_qstat = OUTPUT_QSTAT
output_topgenes = OUTPUT_TOPGENES


for (rsid in rownames(genotype)) {
  X = genotype[rsid ,colnames(expression)]
  Y = (as.numeric(X)-mean(as.numeric(X)))/(var(as.numeric(X)))
  genotype[rsid ,colnames(expression)] = Y
}

genotype = genotype[rownames(genotype) ,colnames(expression)]
expression = expression[rownames(expression),]
normalized_geneexp = expression

X_vector<- as.vector(t(genotype))
Y_vector<- as.vector(t(expression))

print ("File reading complete.")


Qstat <- function(pvals) {
  pvals <- sort(pvals)
  G <- length(pvals) 
  kmax = min(100,G)
  digammaG1 <- digamma(G+1)
  z <- - (log(pvals[1:kmax]) - (digamma(1:kmax) - digammaG1) )
  sumz <- cumsum(z)	
  max(sumz) 
}

calibrateQtest <- function(G=20000,S=1E5) {
  # Generate negative samples for Q statistic
  s10 = max(round(S/100),100)
  qneg <- array(0, dim=S)	
  for (s in 1:S) {
    qneg[s] <- Qstat( runif(G) )
    if ((s %% s10)==0) { print(c("Negative samples:", s, "out of", S)) }
  }
  # Compute the paramters for the exponential extrapolation
  # P-value( Q) = p0 * exp( -lam * Q)
  qneg_sorted <- sort(qneg, decreasing=TRUE)
  n1 <- 0.001*S  # fit exponential in highest q=0.001 quantile
  n2 <- n1/4	
  if (1) {  # NEEDS TO BE TESTED!!
    # do ML in highest 0.001 quantile
    lam = n1/ (sum (qneg_sorted[1:n1]))  
    p0 = n1 / S
  } else {  
    # Do a shitty fitting at two grid points, n1 and n2
    lam <- -log(n1/n2)/(qneg_sorted[n1] - qneg_sorted[n2]) 
    p0 <- n1/S * exp(lam*qneg_sorted[n1])
  }
  Qcal = list()	
  Qcal$ecdf <- ecdf(qneg)	# empirical cumulative distribution function
  Qcal$p0 <- p0
  Qcal$lam <- lam
  Qcal$Qmax <- qneg_sorted[n1] # beyond Qmax we extrapolate the p-value, below we take the empirical p-value
  return(Qcal)
}

Q.test <- function(pvals, Qcal) {
  Qscore <- Qstat(pvals)
  if(Qscore>Qcal$Qmax) {
    res <- 1- Qcal$p0 * exp( -Qcal$lam * Qscore)	
  } else {
    res <- 1 - Qcal$ecdf(Qscore)
  }	
  return(res)
}

Qcal = calibrateQtest()


# The C-function wrapper which returns the array of F-statistic
cfunc <- function(xarr, yarr, num_gene, num_sample) {
  xsize <- length(xarr)
  num_snps <- as.integer(xsize / num_sample)
  fsize <- num_snps * num_gene
  cres <- .C("fit",
             as.double(xarr),
             as.double(yarr),
             as.integer(num_snps),
             as.integer(num_gene),
             as.integer(num_sample),
             double(fsize) )
  fstat <- cres[[6]]
  pval <- 1 - pf(fstat, 1, num_sample - 2)
}

rfunc <- function(xarr, yarr, num_gene, num_sample) {
  xsize <- length(xarr)
  num_snps <- as.integer(xsize / num_sample)
  res <- rep(NA, num_snps * num_gene)
  for (i in 1:num_snps) {
    xstart <- (i - 1) * num_sample + 1
    xend <- i * num_sample
    X <- xarr[xstart : xend]
    for (j in 1:num_gene) {
      ystart <- (j - 1) * num_sample + 1
      yend <- j * num_sample
      Y <- yarr[ystart : yend]
      model <- lm(Y ~ X)
      #print(anova(model))
      #print(summary(model))
      pval <- summary(model)$coefficients[2,4]
      pos <- (i - 1) * num_gene + j
      res[pos] <- pval
    }
  }
  #print(res)
}



# Get the genotype and expression
nsample <- dim(expression)[2]
ngene <- dim(expression)[1]
nsnps <- dim(genotype)[1]

Res_qval_pqval = matrix(0, nsnps, 3)


for (i in 1:nsnps){
  Res_qval_pqval[i,1] = rownames(genotype)[i]
}

#genotype   = as.numeric(X_vector) #runif(nsnps * nsample, 0, 1)
#expression = as.numeric(Y_vector) #runif(ngene * nsample, 0, 1)
genotype = X_vector
expression = Y_vector


# without parallelization, it happens on the master node.
# keep this for debugging.
#start.time <- Sys.time()
#rfunc(genotype, expression, ngene, nsample)
#time.taken <- Sys.time() - start.time
#print(time.taken)

# with C code, no parallelization. only for master node
# keep this for debugging
#start.time <- Sys.time()
#pvals <- cfunc(genotype, expression, ngene, nsample)
#time.taken <- Sys.time() - start.time
#print(time.taken)

ncore <- mpi.universe.size() - 1
nmax <- (as.integer(nsnps / ncore) + 1) * nsample
gsplit = split(genotype, ceiling(seq_along(genotype)/nmax))

# spawn the slaves, send the objects, and get the work done
mpi.spawn.Rslaves(nslaves=ncore)
mpi.scatter.Robj2slave(gsplit)
mpi.bcast.Robj2slave(expression)
mpi.bcast.Robj2slave(nsample)
mpi.bcast.Robj2slave(ngene)
mpi.bcast.Robj2slave(cfunc)
mpi.bcast.cmd( dyn.load("/usr/users/sbanerj/nagial/trans-eqtl/scripts/linear_regression.so") )
final_res <- mpi.remote.exec(cfunc(gsplit, expression, ngene, nsample))
pvals <- unlist(final_res, use.names = FALSE)
mpi.close.Rslaves(dellog = FALSE)

print ("pvals calculated")


Qval = list()
Pval_Qval = list()
for (i in 1:nsnps){
  Qval[i] = Qstat(pvals[ (i-1)*ngene + 1: i*ngene ]) 
  Pval_Qval[i] = Q.test(pvals[ (i-1)*ngene + 1: i*ngene ], Qcal)
  Res_qval_pqval[i,2] = Qstat(pvals[ (i-1)*ngene + 1: i*ngene ])
  Res_qval_pqval[i,3] = Q.test(pvals[ (i-1)*ngene + 1: i*ngene ], Qcal)
}

print ("Qstat calculated")

Res_pval_gene_cutoff = matrix(0, 1, 3)
tempmatrix = matrix(0, 1, 3)
for (i in 1:nsnps){
    if (as.numeric(Res_qval_pqval[i,3]) < p_qval_cutoff){
      for(j in 1:ngene){
	  if (pvals[(i-1)*ngene+j] < pval_cutoff){
	     tempmatrix[1,1] = Res_qval_pqval[i,1]
	     tempmatrix[1,2] = colnames(normalized_geneexp)[j] #change to expression afteryou link to tmp.r
	     tempmatrix[1,3] = pvals[(i-1)*ngene+j]
      	     #Res_qval_pqval[i,1]
	     #colnames(normalized_geneexp)[j]
	     #pvals[(i-1)*ngene+j]	
	     Res_pval_gene_cutoff = rbind(Res_pval_gene_cutoff, tempmatrix)			  
            }
	}
    }

}

print ("Cutoff applied")

Res_pval_gene_cutoff = Res_pval_gene_cutoff[-1,]
write.table(Res_qval_pqval, output_qstat, sep="\t", row.names = F, col.names = F)
write.table(Res_pval_gene_cutoff, output_topgenes, sep="\t", row.names = F, col.names = F)

# quit mpi
mpi.quit()

#sink("outfile.qstat.txt")
#for (rsid in rownames(genotype)){
#   cat("rsid Qval[1] Pval_Qval[1]")
#   cat("\n")
#}
#sink()
#  # split the genotype fo sending to slaves
#  ncore <- mpi.universe.size() - 1
#  cat(sprintf("%d\n", ncore))
#  nmax <- (as.integer(nsnps / ncore) + 1) * nsample
#  gsplit = split(genotype, ceiling(seq_along(genotype)/nmax))
#  
#  # spawn the slaves, send the objects, and get the work done
#  mpi.spawn.Rslaves(nslaves=ncore)
#  mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
#  mpi.scatter.Robj2slave(gsplit)
#  mpi.bcast.Robj2slave(expression)
#  mpi.bcast.Robj2slave(nsample)
#  mpi.bcast.Robj2slave(ngene)
#  mpi.bcast.Robj2slave(cfunc)
#  mpi.bcast.cmd( dyn.load("/mnt/work/transeqtl_gtex/transeqtl_pipeline/linear_regression.so") )
#  final_res <- mpi.remote.exec(cfunc(gsplit, expression, ngene, nsample))
#  mpi.close.Rslaves(dellog = FALSE)

# print the final result
#print(unlist(final_res, use.names = FALSE))


#mpi.bcast.cmd( id <- mpi.comm.rank() )
##mpi.bcast.cmd( ns <- mpi.comm.size() )
##mpi.bcast.cmd( host <- mpi.get.processor.name() )
#mpi.bcast.cmd( dyn.load("linear_regression.so") )
##mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
#mpi.remote.exec(ls(.GlobalEnv))
#mpi.quit()

