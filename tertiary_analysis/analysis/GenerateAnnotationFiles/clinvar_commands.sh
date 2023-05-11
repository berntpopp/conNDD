#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/clinvar/

## clinvar.vcf.gz was dowloaded from the offical Clinvar website
## "chr" was added to each record for compatibility; "|" was replaced with "," for compatibility with SnpSift

wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20200127.vcf.gz
gunzip clinvar_20200127.vcf.gz

(cat clinvar_20200127.vcf | grep "^#" | grep -v "^#CHROM" | sed 's/GRCh38/HG38/g'; \
cat clinvar_20200127.vcf | grep "^#CHROM"; \
cat clinvar_20200127.vcf | grep -v "^#" | grep -v "^NW_009646201.1" | sed 's/^/chr/g' | sed 's/chrMT/chrM/g' | sed 's/|/,/g') > clinvar_20200127.hg38.changed.vcf

## normalize indels and split multiallelic 
gatk LeftAlignAndTrimVariants \
-R ../../../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-V clinvar_20200127.hg38.changed.vcf \
-O clinvar_20200127.hg38.changed.split.vcf \
--max-indel-length 5000 \
--max-leading-bases 5000 \
--split-multi-allelics

bcftools sort -O v clinvar_20200127.hg38.changed.split.vcf -o clinvar_20200127.hg38.vcf

bgzip clinvar_20200127.hg38.vcf
tabix -f -p vcf clinvar_20200127.hg38.vcf.gz

# clean up
rm -rf *.vcf *.idx
