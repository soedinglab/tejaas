#!/bin/bash

GENO=$1
echo $1
python reverse_regression_old.py --genotype $1  --sample data/donor_ids.fam --expression data/expressions/GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt --output out --start 0 --end 100000  --transgeno data/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_filteredTEJAAS_transeqtls.gz --optimize 0 --sigbeta 0.001636581825407

#python reverse_regression_old.py --genotype $1  --sample data/donor_ids.fam --expression PCA_KNN_Correction/main/out  --output out --start 0 --end 100000  --transgeno data/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_filteredTEJAAS_transeqtls.gz --optimize 0 --sigbeta 0.001636581825407
