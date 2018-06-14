# TEJAAS

**T**rans-**E**qtls from **J**oint **A**ssociation **A**nalysi**S**

A tool for finding trans-eQTLs (and cis-eQTLs).

## Description

TEJAAS is a command line tool that performs an association analysis of all genes with every SNP to find trans-eQTLs.
First, it performs a linear regression of every SNP with the gene expression of all `G` genes 
and combine all the `G` p-values to obtain a JPA-score which gives the significance of having multiple genes being associated with a single SNP.
Second, it performs a reverse regression, which is a L<sub>2</sub>-regularized multiple regression with every SNP as the target variable 
and all genes as exaplanatory variables.

## Installation

TEJAAS is written in python and C. To run TEJAAS, you will need
- python version 3.4 or higher,
- Python packages for scientific computing NumPy and SciPy
- C compiler
- Intel MKL library
- any flavor of MPI linked to the Intel MKL library (e.g. OpenMPI)
- Python package mpi4py linked to MPI and MKL

### (Example) Instructions for GWDG cluster at MPIBPC
1. (Optional) Uninstall mpi4py if you already have it installed
2. Load the modules
```
module load intel/compiler/64/2017/17.0.2
module load intel/mkl/64/2017/2.174
module load openmpi/intel/64/1.10.7
```
3. Install mpi4py (so that it will now link to openmpi and mkl)
```
conda install mpi4py
```
4. Make the `C` libraries
```
cd lib
make clean
make -f Makefile.gwdg
cd ..
```
5. Update file path of TEJAAS executable in job submission file located in `devtools/test.jobsub`
6. Submit job
```
cd devtools
bsub < test.jobsub
```
7. (Optional) If you are using interactive node, you can submit the job directly using
```
EXPRFILE="/scratch/sbanerj/trans-eQTL/data/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"
mpirun -n 8 bin/tejaas --gx ${EXPRFILE} --simulate --test --method rr --null maf
```
