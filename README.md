# Tejaas: Discover trans-eQTLs

## Description

Tejaas is a command line tool to find trans-eQTLs from eQTL data.
It is released under the GNU General Public License version 3.

Tejaas is based on the hypothesis that a trans-eQTL should regulate the expression levels of multiple genes.
In brief, it implements two statistical methods to find trans-eQTLs:
- **RR-score (Reverse Regression)**: It performs a multiple linear regression with L<sub>2</sub>-regularization 
using expression levels of all genes to explain the genotype of a candidate SNP.
In contrast to conventional methods, the direction of the regression is reversed, with the gene expressions as explanatory variables.
RR-score is a statistic which estimates whether more genes are required to explain the allele counts of a SNP than expected by chance.
- **JPA-score (Joint P-value Analysis)**: It evaluates the distribution of p-values of the pairwise linear association of a candidate SNP with all available gene expression levels. 
Any null SNP (no trans-effect) will have a uniform distribution of p-values,
while a trans-eQTL will be associated with more genes than expected by chance, leading to overdispersion near zero. 
The JPA-score is a statistic which estimates whether the distribution of p-values is significantly overdispersed near zero.

Additionally, it also implements a non-linear unsupervised confounder correction using k-nearest neighbors called **KNN correction**.

## Dependencies

- Python version 3.4 or higher,
- Intel MKL library
- C compiler
- Python libraries:
  - [NumPy](http://www.numpy.org/) / array operations
  - [SciPy](https://www.scipy.org/) / optimization and other special functions
  - [statsmodel](http://www.statsmodels.org/stable/index.html) / used for ECDF calculation in JPA-score
  - [Pygtrie](https://pypi.org/project/pygtrie/) / used for reading MAF file in RR-score / maf null
  - [mpi4py](https://mpi4py.readthedocs.io/en/stable/) / linked to MPI and MKL for python parallelization
  - [scikit-learn](https://scikit-learn.org/stable/index.html) / used for PCA decomposition in KNN correction

Optional:
- any flavor of MPI linked to the Intel MKL library (e.g. OpenMPI)

You can find examples of getting started here:
- [Example 1 (GWDG Cluster)](https://github.com/soedinglab/tejaas/wiki/GWDG-Cluster)
- [Example 2 (Minion)](https://github.com/soedinglab/tejaas/wiki/Minion2)

## Installation
1. Clone this repository.
2. Compile the C libraries provided in the `lib` subdirectory. Some example makefiles are provided within the `lib` subdirectory.
```
cd lib
make all -f Makefile
```
3. Run Tejaas!
```
bin/tejaas [OPTIONS]
```
See below for valid options or try `bin/tejaas --help`.

#### Run an example to check installation
An example script `test/run_test.sh` is provided to check the installation.
Open the script in your favorite editor and modify the variables `NCORE` and `DATA_DIR`.
The script will download some example input files in the `DATA_DIR` directory and run Tejaas on `NCORE` cores.
The output will be created in `DATA_DIR`. 
Check if the output matches with the results provided in the `test/gold` subdirectory.
```
cd test
./run_test.sh
```

## Input Files
- Gene expression file
- Genotype file
  - VCF
  - Oxford
  - Dosage
- GENCODE file
- Population minor allele frequency

## Tejaas [OPTIONS]

Option&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Argument | Description | Priority | Default value
:--- | :--- | :---        |:---      | :--
`--vcf`       | `FILEPATH` | Input genotype file in vcf.gz format | Required (vcf or oxf) | --
`--oxf`       | `FILEPATH` | Input genotype file in Oxford format | Required (vcf or oxf) | --
`--dosage`    |            | Flag for reading dosage files. The file is specified with the `--oxf` option, e.g. `--oxf FILEPATH --dosage` | Optional | `False` 
`--fam`       | `FILEPATH` | Input fam file for samples names of Oxford genotype | Optional | -- 
`--chrom`     | `INT`      | Chromosome number of the genotype file | Required | -- 
`--include-SNPs` | `START:END` | Colon-separated index of SNPs to be included | Optional | -- 
`--gx`        | `FILEPATH` | Input gene expression file for trans-eQTL discovery | Required | --
`--gxcorr`    | `FILEPATH` | Input gene expression file for target gene discovery | Optional | `--gx` file
`--gxfmt`     | `OPTION`   | Input gene expression file format (see format details below). Supported options: `gtex`, `cardiogenics`, `geuvadis` | Optional | `gtex`
`--gtf`       | `FILEPATH` | Input GTF file from GENCODE to read gene Ensembl IDs. Used for selecting biotypes and getting genomic locations. | Required | --
`--trim`      |            | Flag to trim version number from GENCODE Ensembl IDs | Optional | `False`
`--biotype`   | `OPTION`   | Which biotypes to select from the GTF file. Supported options: `protein_coding`, `lncRNA`. | Optional | `protein_coding lncRNA`
`--outprefix` | `STRING`   | Full path to output file names. The extensions are generated by Tejaas. | Optional | `out`  
`--method`    | `OPTION`   | Name of method to run. Supported options: `jpa` or `rr` | Optional | `rr`
`--null`      | `OPTION`   | Null model to use for RR-score. Supported options: `perm` or  `maf` | Optional | `perm` 
`--cismask`   |            | Flag to mask cis-Genes within a window  for each candidate SNP. Gene positions are obtained from the GENCODE annotation file.  | Optional | `False`
`--window`    | `FLOAT`    | Window (number of base pairs) used for masking cis genes | Optional | 1e6 
`--prior-sigma`| `FLOAT`   | Standard deviation of the normal prior for reverse multiple linear regression | Optional | 0.1
`--knn`       | `INT`      | Number of neighbours for KNN (use 0 if you do not want to use KNN) | Optional | 0
`--psnpthres` | `FLOAT`    | Target genes will be reported only for trans-eQTLs below this threshold p-value for RR/JPA-score | Optional | 0.0001
`--pgenethres`| `FLOAT`    | Target genes will be reported only if their association with trans-eQTLs are below this threshold p-value | Optional | 0.05
`--jpanull`   | `FILEPATH` | File containing list of null model JPA-scores | Optional | -- 
`--maf-file`  | `FILEPATH` | Read minor allele frequency (MAF) of SNPs from this file, e.g. to read population MAF for `maf` null (see documentation for file format) | Optional | -- 
`--shuffle`   |            | Flag to randomly shuffle the genotypes to obtain a null distribution | Optional | `False`
`--shuffle-with` | `FILEPATH` | Shuffle the genotypes in the same order of donor IDs specified in `FILEPATH` | Optional | -- 
`--test`      |            | Flag to do test run | Optional | -- 

## Usage Examples

1. For quick start or installation check, run Tejaas with all default options:
```
bin/tejaas --vcf ${VCFFILE} --chrom ${CHRM} --gx ${GXFILE} --gtf ${GTFFILE} --cismask --outprefix ${OUTPREFIX}
```
This will create RR-scores at &gamma;=0.1 and masking all genes within 1Mb of each SNP. The p-values will be computed from the permuted null model.
Default format for the gene expression is the same as the GTEx format, and default gtf file is the [GENCODE v26](https://www.gencodegenes.org/human/release_26.html) release.
For target gene discovery, it will use the same file as used for trans-eQTL discovery.

2. Example of running Tejaas RR-score.
We recommend using the `perm` null model for calcuting p-values from the RR-score
and a separate confounder-corrected gene expression file for target gene discovery.
In this example, RR-score is calculated for first 1000 SNPs excluding the first 20 `--include-SNPs 21:1000`.
KNN correction is performed with 20 nearest neighbors `--knn 20`.
All cis-genes within +2MB and -2Mb are masked during analysis `--cismask --window 2e6`.
RR-score calculation uses a prior normal distribution with standard deviation of 0.05 `--prior-sigma 0.05`.
The output reports target genes only for SNPs with p-value < 1e-6 `--psnpthres 0.000001`.
Here, `GXFILE` is the raw gene expression file, `GXCORRFILE` is the confounder-corrected gene expression file,
`VCFFILE` is the genotype file in `.vcf.gz` format and `GTFFILE` is the GENCODE annotation file.
```
mpirun -n 8 bin/tejaas --vcf ${VCFFILE} --chrom ${CHRM} --include-SNPs 21:1000 --gx ${GXFILE} --gxcorr ${GXCORRFILE} \
                       --gxfmt gtex --gtf ${GTFFILE} --trim  --outprefix ${OUTPREFIX} \
                       --cismask --window 2e6 --psnpthres 0.000001 \
                       --knn 20 --method rr --null perm --prior-sigma 0.05
```

3. Example of running JPA-score with no KNN correction. 
Empirical p-values are calculated from the null scores loaded from `NULLFILE` specified by the `--jpanull` option.
If `NULLFILE` does not exist, then it will create `100000` null scores and write them in the `NULLFILE` before calculating JPA-scores.
If `--jpanull` option is not used, then p-values for the JPA-scores are calculated from an analytical construction of null model.
```
mpirun -n 8 bin/tejaas --vcf ${VCFFILE} --chrom ${CHRM} --include-SNPs 1:100 \
                       --gx ${GXFILE} --gxfmt gtex --gtf ${GTFFILE} --outprefix ${OUTPREFIX} \
                       --knn 0 --method jpa --jpanull ${NULLFILE}
```

4. Example of parallelizing job submission.
```
NMAX=20000 # number of SNPs per job
for CHRM in $( seq 1 22 ); do
    VCFFILE="file_path_here_${CHRM}.vcf.gz"
    NTOT=$( calculate_no_of_SNPs_in_this_chromosome )
    NJOB=$( echo $(( (NTOT + NMAX - 1)/NMAX )) )
    for (( i=0; i < ${NJOB}; i++ )); do
        STARTSNP=$(( NMAX * i + 1 ))
        ENDSNP=$(( NMAX * (i + 1) ))
        if [ ${ENDSNP} -gt ${NTOT} ]; then
            ENDSNP=${NTOT}
        fi
        mpirun -n 8 bin/tejaas --vcf ${VCFFILE} --chrom ${CHRM} --include-SNPs ${STARTSNP}:${ENDSNP} --gx ${GXFILE} --gxcorr ${GXCORRFILE} \
                               --gxfmt gtex --gtf ${GTFFILE} --trim  --outprefix ${OUTPREFIX} \
                               --cismask --psnpthres 0.000001 --knn 20 --method rr --null perm --prior-sigma 0.05
    done
done
```

## Contributors
<a href="https://github.com/soedinglab/tejaas/graphs/contributors">
  <img src="https://contributors-img.web.app/image?repo=soedinglab/tejaas" />
</a>
