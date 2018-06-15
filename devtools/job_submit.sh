#!/bin/bash

INPUTDIR="${HOME}/gwas-eQTL/gtex_data"
OUTDIRBASE="${HOME}/trans-eQTL/gtex_nopeer_lmcorr_rr_perm_sb005"
JOBSUBDIR="${OUTDIRBASE}/jobsub"
TEJAAS="${HOME}/trans-eQTL/codebase/tejaas/bin/tejaas"
CWD=`pwd`

MAX_NSNP=100000
SNPCUT=0.05
GENCUT=0.05

if [ ! -d ${JOBSUBDIR} ]; then
    mkdir -p ${JOBSUBDIR}
fi
    
for j in {1..2}; do

    echo "Submitting jobs for Chromosome ${j}."

    GTFILE="${INPUTDIR}/dosages/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr${j}.gz"
    GXFILE="${INPUTDIR}/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"
    DONORS="${INPUTDIR}/dosages/donor_ids.fam"
    GENINF="${INPUTDIR}/gencode.v19.annotation.gtf"
    JOBPREFIX="chr${j}"

    ## do not change below

    OUTDIR="${OUTDIRBASE}/chr${j}"
    
    if [ ! -d ${OUTDIR} ]; then
        mkdir -p ${OUTDIR}
    fi   

    TOTALSNPS=`zcat ${GTFILE} | wc -l`
    NJOBS=$(echo $(( TOTALSNPS/MAX_NSNP )))
    
    for (( i=0; i <= ${NJOBS}; i++ )); do 
        INDEX=`echo $i | awk '{printf "%03d", $1}'`
        JOBNAME="${JOBPREFIX}_${INDEX}"

        STARTSNP=$(( MAX_NSNP * i + 1 ))
        ENDSNP=$(( MAX_NSNP * (i+1) ))
        if [ $ENDSNP -gt $TOTALSNPS ]; then
            ENDSNP=${TOTALSNPS}
        fi
        INCSNP="${STARTSNP}:${ENDSNP}"

        OUTPRF="${OUTDIR}/chunk${INDEX}"

        # create the job submission file
        sed "s|_JOBNAME|${JOBNAME}|g;
             s|_GTFILE_|${GTFILE}|g;
             s|_FAMFIL_|${DONORS}|g;
             s|_GXFILE_|${GXFILE}|g;
             s|_GENINF_|${GENINF}|g;
             s|_PREFIX_|${OUTPRF}|g;
             s|_ST_END_|${INCSNP}|g;
             s|_TEJAAS_|${TEJAAS}|g;
             s|_SNPCUT_|${SNPCUT}|g;
             s|_GENCUT_|${GENCUT}|g;" master.jobsub > ${JOBSUBDIR}/${JOBNAME}.bsub

        # Submit the job
        cd ${JOBSUBDIR}
        #bsub < ${JOBNAME}.bsub
        cd ${CWD}
    done

done
