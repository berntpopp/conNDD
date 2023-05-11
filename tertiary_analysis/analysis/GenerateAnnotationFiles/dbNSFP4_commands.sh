#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/dbNSFP/

## Manually download dbNSFP database from https://drive.google.com/open?id=1BNLEdIc4CjCeOa7V7Z8n8P8RHqUaF5GZ on 2020-01-10

## Uncompress
unzip dbNSFP4.0a.zip

# Create a single file version
(zcat dbNSFP4.0a_variant.chr1.gz | head -n 1 ; zcat dbNSFP4.0a_variant.chr1.gz dbNSFP4.0a_variant.chr2.gz dbNSFP4.0a_variant.chr3.gz dbNSFP4.0a_variant.chr4.gz dbNSFP4.0a_variant.chr5.gz dbNSFP4.0a_variant.chr6.gz dbNSFP4.0a_variant.chr7.gz dbNSFP4.0a_variant.chr8.gz dbNSFP4.0a_variant.chr9.gz dbNSFP4.0a_variant.chr10.gz dbNSFP4.0a_variant.chr11.gz dbNSFP4.0a_variant.chr12.gz dbNSFP4.0a_variant.chr13.gz dbNSFP4.0a_variant.chr14.gz dbNSFP4.0a_variant.chr15.gz dbNSFP4.0a_variant.chr16.gz dbNSFP4.0a_variant.chr17.gz dbNSFP4.0a_variant.chr18.gz dbNSFP4.0a_variant.chr19.gz dbNSFP4.0a_variant.chr20.gz dbNSFP4.0a_variant.chr21.gz dbNSFP4.0a_variant.chr22.gz dbNSFP4.0a_variant.chrX.gz dbNSFP4.0a_variant.chrY.gz dbNSFP4.0a_variant.chrM.gz | grep -v "^#" ) > dbNSFP4.0a.txt

# Compress using block-gzip algorithm
bgzip dbNSFP4.0a.txt

# Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFP4.0a.txt.gz


####################
## for hg19 version 
# Replace coordinates by columns 7 and 8 (hg19 coordinates) and sort by those coordinates
(zcat dbNSFP4.0a.txt.gz | head -n 1 | sed -e 's/hg19/hg38/g'; zcat dbNSFP4.0a.txt.gz | sed 1,1d | awk '{FS="\t"; OFS="\t"} $8 == "." { next; } $8 != "." { s=$1; $1=$8; $8=s; t=$2; $2=$9; $9=t; print; }' | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n) > dbNSFP4.0a_hg19-sorted.txt

# Compress using block-gzip algorithm
nohup bgzip dbNSFP4.0a_hg19-sorted.txt &

# Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFP4.0a_hg19-sorted.txt.gz