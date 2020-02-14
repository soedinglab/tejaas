#!/bin/sh

NCORE="4" # number of cores for running Tejaas
DATA_DIR="data" # directory for downloading data and output results

#====================
# DO NOT CHANGE BELOW
#====================

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load intel/mpi/64/2017/2.174

FTP_URL="http://wwwuser.gwdg.de/~compbiol/tejaas/example"
RUN_PATH="../bin/tejaas"
GENOFILE="${DATA_DIR}/GEUVADIS.chr22.vcf.gz"
GXPRFILE="${DATA_DIR}/GEUVADIS.expression.mod.txt"
GENEINFO="${DATA_DIR}/gencode.v19.annotation.gtf.gz"
OUTPREFIX="data/result"
INCSTRING=0:2000
CHROM=22

if [ ! -d ${DATA_DIR} ]; then mkdir -p ${DATA_DIR}; fi

for FPATH in ${GENOFILE} ${GXPRFILE} ${GENEINFO}; do
    if [ ! -f ${FPATH} ]; then
        FNAME=$( basename ${FPATH} )
        wget -P ${DATA_DIR}/ ${FTP_URL}/${FNAME}
    fi
done

mpirun -n ${NCORE} ${RUN_PATH} --vcf ${GENOFILE} --include-SNPs ${INCSTRING} --chrom ${CHROM} \
                               --gx ${GXPRFILE} --gxcorr ${GXPRFILE} --gtf ${GENEINFO} --outprefix ${OUTPREFIX} \
                               --method rr --prior-sigma 0.1 --knn 30 --cismask --null perm --psnpthres 0.001 --pgenethres 0.001
