#!/bin/sh
GENOFILE=${HOME}/gwas-eQTL/gtex_data/dosages/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr1.gz
SAMPFILE=${HOME}/gwas-eQTL/gtex_data/dosages/donor_ids.fam
GXPRFILE=/scratch/sbanerj/trans-eQTL/data/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt
OUTPRFIX=trial
INCSTRNG=0:20000
RUN_PATH=${HOME}/trans-eQTL/codebase/tejaas/bin/tejaas
SNPTHRES=0.001
GENTHRES=0.001
GENEINFO=${HOME}/gwas-eQTL/gtex_data/gencode.v19.annotation.gtf.gz

#source activate py35
mpirun -n 8 ${RUN_PATH} --oxf          ${GENOFILE} \
                        --fam          ${SAMPFILE} \
                        --gx           ${GXPRFILE} \
                        --gtf          ${GENEINFO} \
                        --method       rr \
                        --null         perm \
                        --outprefix    ${OUTPRFIX} \
                        --include-SNPs ${INCSTRNG} \
                        --psnpthres    ${SNPTHRES} \
                        --pgenethres   ${GENTHRES} \
                        --prior-sigma  0.005
