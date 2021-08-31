#!/bin/sh

OUTDIR="$1"
NCORE="$2" # number of cores for running Tejaas
DATA_DIR="${OUTDIR}/data" # directory for downloading data and output results
EXMPL_DIR=$( pwd ) # this code is always run from the current working directory

#====================
# DO NOT CHANGE BELOW
#====================

#FTP_URL="http://wwwuser.gwdg.de/~compbiol/tejaas/example"
GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
RUN_PATH="tejaas"
GENOFILE="${DATA_DIR}/genotype.vcf.gz"
GXPRFILE="${DATA_DIR}/expression.txt"
GENEINFO="${DATA_DIR}/gencode.v19.annotation.gtf.gz"
OUTPREFIX_RR="${DATA_DIR}//result"
OUTPREFIX_JPA="${DATA_DIR}/result_fr"
NULLFILE="${DATA_DIR}/result_jpanull.txt"
INCSTRING=0:110
CHROM=22

if [ ! -d ${DATA_DIR} ]; then mkdir -p ${DATA_DIR}; fi

if [ ! -f ${GENOFILE} ]; then cp ${EXMPL_DIR}/genotype.vcf.gz ${GENOFILE} ; fi
if [ ! -f ${GXPRFILE} ]; then tar -zxf ${EXMPL_DIR}/expression.txt.tar.gz -C ${DATA_DIR}/ ; fi
if [ ! -f ${GENEINFO} ]; then wget -P ${DATA_DIR}/ ${GTF_URL} ; fi
if [   -f ${NULLFILE} ]; then rm -f ${NULLFILE} ; fi

mpirun -n ${NCORE} ${RUN_PATH} --vcf ${GENOFILE} --include-SNPs ${INCSTRING} --chrom ${CHROM} \
                               --gx ${GXPRFILE} --gxcorr ${GXPRFILE} --gtf ${GENEINFO} --trim --outprefix ${OUTPREFIX_RR} \
                               --method rr --prior-sigma 0.1 --knn 30 --cismask --null perm --psnpthres 0.1 --pgenethres 0.1 

#mpirun -n ${NCORE} ${RUN_PATH} --vcf ${GENOFILE} --include-SNPs ${INCSTRING} --chrom ${CHROM} \
#                               --gx ${GXPRFILE} --gxfmt gtex --gtf ${GENEINFO} --trim --outprefix ${OUTPREFIX_JPA} \
#                               --knn 0 --method jpa --jpanull ${NULLFILE} --jpanull-iter 1000 --seed 100 \
#                               --psnpthres 0.1 --pgenethres 0.1
