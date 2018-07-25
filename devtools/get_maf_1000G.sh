#!/bin/bash

for chrom in {1..22}; do
    filename="ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    outfile="1000G_phase3_v5a_20130502_snpinfo_EUR_chr${chrom}.txt"

    if [ -f ${outfile} ]; then
        rm -f ${outfile}
    fi

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${filename}
    zcat ${filename} | sed '/^#/ d' | awk '{print $1, $2, $3, $4, $5, $8}' > ALL.chr${chrom}.snpinfo.txt
    cat ALL.chr${chrom}.snpinfo.txt | awk '{snp="False"; eur_af=0;
                                            split($6, infoarr, ";"); 
                                            for (i in infoarr) {split(infoarr[i], m, "="); 
                                                                if (m[1] == "VT") {if (m[2] == "SNP"){snp="True"; } };
                                                                if (m[1] == "EUR_AF") {eur_af = m[2];};
                                                               };
                                            if (snp == "True" && eur_af !~ /.*,.*$/ && eur_af > 0 && eur_af < 1) {print $1, $2, $3, $4, $5, eur_af;}; 
                                           }' >> ${outfile}
done
