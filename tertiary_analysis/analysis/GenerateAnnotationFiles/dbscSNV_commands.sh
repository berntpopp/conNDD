#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/dbscSNV/

## manually download dbscSNV database from https://drive.google.com/file/d/0B60wROKy6OqcZkw2bWt2TGU5NDA/view on 2020-01-11

## Uncompress
unzip dbscSNV.zip

(echo "##fileformat=VCFv4.2"; echo "##fileDate=2020-01-11"; echo "##reference=HG19"; echo "##dbscSNV_version=v.1.1"; echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; echo "##INFO=<ID=dbscSNV_ada_score,Number=A,Type=Float,Description=\"Field 'ada_score' from dbscSNV\">"; echo "##INFO=<ID=dbscSNV_rf_score,Number=A,Type=Float,Description=\"Field 'rf_score' from dbscSNV\">"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; cat dbscSNV.chr1 dbscSNV.chr2 dbscSNV.chr3 dbscSNV.chr4 dbscSNV.chr5 dbscSNV.chr6 dbscSNV.chr7 dbscSNV.chr8 dbscSNV.chr9 dbscSNV.chr10 dbscSNV.chr11 dbscSNV.chr12 dbscSNV.chr13 dbscSNV.chr14 dbscSNV.chr15 dbscSNV.chr16 dbscSNV.chr17 dbscSNV.chr18 dbscSNV.chr19 dbscSNV.chr20 dbscSNV.chr21 dbscSNV.chr22 dbscSNV.chrX dbscSNV.chrY | grep -v "chr" | \
cut -f 1,2,3,4,15,16 | awk 'BEGIN{OFS="\t";OFS="\t";} {print "chr"$1,$2,".",$3,$4,"100","PASS","dbscSNV_ada_score="$5";dbscSNV_rf_score="$6}' | sort -k1,1d -k2,2n) \
 > dbscSNV.hg19.vcf

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz

CrossMap.py vcf hg19ToHg38.over.chain dbscSNV.hg19.vcf ../../../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa dbscSNV.liftover.hg38.vcf

(grep "#" dbscSNV.liftover.hg38.vcf; grep -v "#" dbscSNV.liftover.hg38.vcf | grep -E "chr[0-9A-Za-z]+_") > dbscSNV.hg38.otherchrom
(grep "##" dbscSNV.liftover.hg38.vcf; grep "#" dbscSNV.liftover.hg38.vcf | grep -v "##"; grep -v "#" dbscSNV.liftover.hg38.vcf | grep -Ev "chr[0-9A-Za-z]+_") > dbscSNV.hg38.canonicalchrom.vcf

bcftools sort -O v dbscSNV.hg38.canonicalchrom.vcf -o dbscSNV.hg38.vcf

bgzip dbscSNV.hg38.vcf
tabix -f -p vcf dbscSNV.hg38.vcf.gz

# clean up
rm -rf *.vcf *.zip *.chr[0-9] *.chr[A-Z] *.chr[0-9][0-9]
