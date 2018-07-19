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
- Python packages for scientific computing [NumPy](http://www.numpy.org/) and [SciPy](https://www.scipy.org/)
- Python package for Trie data structure, [Pygtrie](https://pypi.org/project/pygtrie/)
- C compiler
- Intel MKL library
- any flavor of MPI linked to the Intel MKL library (e.g. OpenMPI)
- Python package mpi4py linked to MPI and MKL

You can find examples of getting started here:
- [Example 1 (GWDG Cluster)](https://github.com/soedinglab/tejaas/wiki/GWDG-Cluster)
- [Example 2 (Minion)](https://github.com/soedinglab/tejaas/wiki/Minion2)

## Example Command

For quick check (dev only):
```
EXPRFILE="/scratch/sbanerj/trans-eQTL/data/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"
mpirun -n 8 bin/tejaas --gx ${EXPRFILE} --simulate --method rr --null perm --simparams 0.1 0.1 1000 0 0.001 --test
```
