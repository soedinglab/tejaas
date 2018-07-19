#!/bin/sh
#BSUB -J _JOBNAME
#BSUB -q mpi
#BSUB -W 48:00
#BSUB -n 8
#BSUB -o _JOBNAME.o%J
#BSUB -e _JOBNAME.e%J
#BSUB -a intelmpi
#BSUB -R "span[hosts=1]"

module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load openmpi/intel/64/1.10.7

GENOFILE=_GTFILE_
SAMPFILE=_FAMFIL_
GXPRFILE=_GXFILE_
OUTPRFIX=_PREFIX_
INCSTRNG=_ST_END_
RUN_PATH=_TEJAAS_
SNPTHRES=_SNPCUT_
GENTHRES=_GENCUT_
GENEINFO=_GENINF_

#source activate py35
mpirun -n 8 ${RUN_PATH} --oxf          ${GENOFILE} \
                        --fam          ${SAMPFILE} \
                        --gx           ${GXPRFILE} \
                        --gtf          ${GENEINFO} \
                        --method       jpa-rr \
                        --null         perm \
                        --outprefix    ${OUTPRFIX} \
                        --include-SNPs ${INCSTRNG} \
                        --psnpthres    ${SNPTHRES} \
                        --pgenethres   ${GENTHRES} \
                        --prior-sigma  0.0006	\
			--isdosage	\
                        --oxf-columns  6

