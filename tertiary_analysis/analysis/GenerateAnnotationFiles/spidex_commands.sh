#!/bin/sh

##2200-01-10
#########

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/spidex/

# download files
wget http://www.openbioinformatics.org/annovar/download/IlvUMvrpPT/hg19_spidex.zip
unzip hg19_spidex.zip

(echo "##fileformat=VCFv4.2"; echo "##fileDate=2020-01-10"; echo "##reference=HG19"; echo "##SPIDEX_version=v.1.0"; echo "##INFO=<ID=spidex_dpsi_max_tissue,Number=A,Type=Float,Description=\"Field 'dpsi_max_tissue' from hg19_spidex.txt\">"; echo "##INFO=<ID=spidex_dpsi_zscore,Number=A,Type=Float,Description=\"Field 'dpsi_zscore' from hg19_spidex.txt\">"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; grep -v "#" hg19_spidex.txt | awk 'BEGIN{OFS="\t";OFS="\t";} {print "chr"$1,$2,".",$4,$5,"100","PASS","spidex_dpsi_max_tissue="$6";spidex_dpsi_zscore="$7}') > hg19_spidex.vcf

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip hg19ToHg38.over.chain.gz

CrossMap.py vcf hg19ToHg38.over.chain hg19_spidex.vcf ../../../ref/GRCh38_full_analysis_set_plus_decoy_hla.fa hg19_spidex.hg38.vcf

(grep "#" hg19_spidex.hg38.vcf; grep -v "#" hg19_spidex.hg38.vcf | grep -E "chr[0-9A-Za-z]+_") > hg19_spidex.hg38.otherchrom

(grep "##" hg19_spidex.hg38.vcf; echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; grep "#" hg19_spidex.hg38.vcf | grep -v "##"; grep -v "#" hg19_spidex.hg38.vcf | grep -Ev "chr[0-9A-Za-z]+_") > hg19_spidex.hg38.canonicalchrom.vcf

bcftools sort -O v hg19_spidex.hg38.canonicalchrom.vcf -o spidex.hg38.vcf

bgzip spidex.hg38.vcf
tabix -f -p vcf spidex.hg38.vcf.gz

# clean up
rm *.vcf *.zip *.docx *.idx *.txt