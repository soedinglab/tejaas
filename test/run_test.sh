#!/bin/sh

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load intel/mpi/64/2017/2.174

RUN_PATH=../bin/tejaas
GENOFILE=data/GEUVADIS.chr22.vcf.gz
# SAMPFILE=data/
GXPRFILE=data/GEUVADIS.expression.mod.txt
GENEINFO=/cbscratch/franco/datasets/GENCODE/gencode.v19.annotation.gtf.gz

TJMETHOD=jpa-rr
NULLMODL=perm
OUTPRFIX=data/result
INCSTRNG=0:2000
SNPTHRES=0.001
GENTHRES=0.001
SBETA=0.1
CHROM=22

mpirun -n 8 ${RUN_PATH} --vcf          ${GENOFILE} \
                        --gx           ${GXPRFILE} \
                        --gtf          ${GENEINFO} \
                        --method       ${TJMETHOD} \
                        --null         ${NULLMODL} \
                        --outprefix    ${OUTPRFIX} \
                        --include-SNPs ${INCSTRNG} \
                        --psnpthres    ${SNPTHRES} \
                        --pgenethres   ${GENTHRES} \
                        --prior-sigma  ${SBETA} \
                        --chrom        ${CHROM} \
                        --dosage                \
                        --knn 30                \
                        --cismask               

# --fam          ${SAMPFILE} \