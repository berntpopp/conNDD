#!/bin/sh

#############
## 2020-01-27

# go to annovar installation folder
cd "$(config_get annovar_folder)"

# download files
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz

tar -xf annovar.latest.tar.gz

# use annotate_variation to downlaod gene4denovo201907 database
./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gene4denovo201907 humandb/

cat humandb/hg38_gene4denovo201907.txt | grep -v "#" > humandb/hg38_gene4denovo201907.noheader.txt

cat humandb/hg38_gene4denovo201907.noheader.txt | awk '{print "samtools faidx "$(config_get work_directory)"/analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa "$1":"$2-1"-"$2-1}' | sh - | paste - - | cut -f2 > humandb/hg38_gene4denovo201907.noheader.previousbase.txt

cat humandb/hg38_gene4denovo201907.noheader.txt | cut -f1,2,3,4,5,6,7,8,9,11 > humandb/hg38_gene4denovo201907.noheader.noline10.txt

(echo "##fileformat=VCFv4.2"; echo "##fileDate=2020-01-27"; echo "##reference=HG38"; echo "##hg38_gene4denovo_version=201907"; echo "##INFO=<ID=gene4denovo_PubmedID,Number=A,Type=Float,Description=\"Field 'Pubmed ID' from ANNOVAR hg38_gene4denovo201907 txt file\">"; echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; paste  humandb/hg38_gene4denovo201907.noheader.noline10.txt humandb/hg38_gene4denovo201907.noheader.previousbase.txt | awk 'BEGIN{OFS="\t";OFS="\t";} { if ( $4 ~ /-/ || $5 ~ /-/) {gsub("-","",$4); gsub("-","",$5); print $1,$2-1,".",$11$4,$11$5,"100","PASS","gene4denovo_PubmedID="$10} else {print $1,$2,".",$4,$5,"100","PASS","gene4denovo_PubmedID="$10} }') > humandb/hg38_gene4denovo201907.noheader.vcf

bcftools norm --check-ref x -f "$(config_get work_directory)"/analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa -o humandb/hg38_gene4denovo201907.norm.vcf humandb/hg38_gene4denovo201907.noheader.vcf

bgzip humandb/hg38_gene4denovo201907.norm.vcf
tabix -f -p vcf humandb/hg38_gene4denovo201907.norm.vcf.gz
