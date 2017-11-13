#!/bin/bash

INPUTDIR="${HOME}/gwas-eQTL/gtex_data"
SCRIPTDIR="${HOME}/trans-eQTL/codebase/tejaas/main"
OUTDIRBASE="${HOME}/trans-eQTL/gtex_peer"
JOBSUBDIR="${OUTDIRBASE}/jobsub"
CWD=`pwd`

P_QVAL_CUTOFF="0.05"
PVAL_CUTOFF="0.05"
MAX_NSNP=10000

if [ ! -d ${JOBSUBDIR} ]; then
    mkdir -p ${JOBSUBDIR}
fi
    
for j in {1..22}; do

    OUTDIR="${OUTDIRBASE}/chr${j}"
    
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi

    ## do not change below
    
    TRANSEQTL="${SCRIPTDIR}/trans_eqtl.py"
    GENOTYPE="${INPUTDIR}/dosages/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${j}.gz"
    GENEEXP="${INPUTDIR}/gtex_wholeblood_normalized_PEER_residuals.txt"
    DONORIDS="${INPUTDIR}/dosages/donor_ids.fam"
    TOTALSNPS=`zcat ${GENOTYPE} | wc -l`
    
    NJOBS=$(echo $(( TOTALSNPS/MAX_NSNP )))
    
    for (( i=0; i <= ${NJOBS}; i++ )); do 
        INDEX=`echo $i | awk '{printf "%03d", $1}'`
        STARTSNP=$(( MAX_NSNP * i + 1 ))
        ENDSNP=$(( MAX_NSNP * (i+1) ))
        if [ $ENDSNP -gt $TOTALSNPS ]; then
            ENDSNP=${TOTALSNPS}
        fi
        OUTPUT_QSTAT="${OUTDIR}/chr${j}_${INDEX}"
        JOBPREFIX="chr${j}"
        JOBNAME="${JOBPREFIX}_${INDEX}"
        # create the job submission file
        sed "s|_JOBNAME_|${JOBNAME}|g;
           15s|_trn_eqtl|${TRANSEQTL}|g;
           16s|_genotype|${GENOTYPE}|g;
           17s|_gene_exp|${GENEEXP}|g;
           18s|_donor_id|${DONORIDS}|g;
           19s|_startsnp|${STARTSNP}|g;
           20s|_endsnp__|${ENDSNP}|g;
           21s|_pqvalcut|${P_QVAL_CUTOFF}|g;
           22s|_pvalcut_|${PVAL_CUTOFF}|g;
           23s|_qstatout|${OUTPUT_QSTAT}|g;" master.jobsub > ${JOBSUBDIR}/${JOBNAME}.bsub
    
        # Submit the job
        cd ${JOBSUBDIR}
        bsub < ${JOBNAME}.bsub
        cd ${CWD}
    done

done
