#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/cadd/

# download files
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz.tbi
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz.tbi

#############
# merge, convert to VCF (awk, normalize to reference (bcftools)

(echo "##fileformat=VCFv4.2"; echo "##fileDate=2020-01-27"; echo "##reference=HG38"; echo "##CADD_version=1.5"; echo "##INFO=<ID=CADDv15_RawScore,Number=A,Type=Float,Description=\"Field 'RawScore' from CADD tsv files\">"; echo "##INFO=<ID=CADDv15_PHRED,Number=A,Type=Float,Description=\"Field 'PHRED' from CADD tsv files\">"; echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; zcat whole_genome_SNVs.tsv.gz InDels.tsv.gz ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz | grep -v "#" | awk 'BEGIN{OFS="\t";OFS="\t";} {print "chr"$1,$2,".",$3,$4,"100","PASS","CADDv15_RawScore="$5";CADDv15_PHRED="$6}') > CADDv15_all.hg38.canonicalchrom.vcf

bcftools sort -m 24000M -O v CADDv15_all.hg38.canonicalchrom.vcf -o CADDv15_all.hg38.vcf

bgzip CADDv15_all.hg38.vcf
tabix -f -p vcf CADDv15_all.hg38.vcf.gz

# clean up
rm CADDv15_all.hg38.canonicalchrom.vcf
rm ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz
rm ALL.TOPMed_freeze5_hg38_dbSNP.tsv.gz.tbi
rm InDels.tsv.gz.tbi
rm InDels.tsv.gz
rm whole_genome_SNVs.tsv.gz.tbi
rm whole_genome_SNVs.tsv.gz

#########################################
# subset to target for faster annotation:
sort -k 1,1 -k2,2n -k3,2n  "$(config_get work_directory)"/analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.bed > S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed

bcftools view --threads 24 -O z -l 1 -T S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed CADDv15_all.hg38.vcf.gz -o CADDv15_all.hg38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz

tabix -f -p vcf CADDv15_all.hg38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz