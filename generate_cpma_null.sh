#!/bin/sh


module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load intel/mpi/64/2017/2.174

# source $HOME/miniconda3/envs/py36/bin/activate py36

RUN_PATH=/usr/users/fsimone/tejaas/bin/tejaas_cpma
GENOFILE=/cbscratch/franco/datasets/gtex/genotypes/dosages_allsamples/GTEx_all_chrs_dosages.gz
SAMPFILE=/cbscratch/franco/datasets/gtex/gtex.sample
GXPRFILE=/cbscratch/franco/datasets/gtex/expression/norm_lmcorrected/gtex.normalized.expression.lmcorrected.ms.txt.gencode_filtered
GENEINFO=/cbscratch/franco/datasets/GENCODE/gencode.v19.annotation.gtf.gz

TJMETHOD=jpa
NULLMODL=perm
OUTPRFIX=/cbscratch/franco/tejaas_output/tests/cpma_null_test
INCSTRNG=1:10000
SNPTHRES=0.001
GENTHRES=0.001
SBETA=0.01
# CHROM=22
MAFFILE=None
EXTRAFLAGS="--dosage"

if [ "${NULLMODL}" = "maf" ]; then
    EXTRAFLAGS="$EXTRAFLAGS --maf-file ${MAFFILE}"
fi


mpirun -n 8 ${RUN_PATH} --oxf          ${GENOFILE} \
                        --fam          ${SAMPFILE} \
                        --gx           ${GXPRFILE} \
                        --gtf          ${GENEINFO} \
                        --method       ${TJMETHOD} \
                        --null         ${NULLMODL} \
                        --outprefix    ${OUTPRFIX} \
                        --psnpthres    ${SNPTHRES} \
                        --pgenethres   ${GENTHRES} \
                        --prior-sigma  ${SBETA} \
                        ${EXTRAFLAGS}
                        
                        # --include-SNPs ${INCSTRNG} \
                        # --chrom        ${CHROM} \

