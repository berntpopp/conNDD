#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg38/gnomAD/

# download files
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz.tbi
wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz
wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz.tbi

#########################################
##subset to target for faster annotation:
sort -k 1,1 -k2,2n -k3,2n  "$(config_get work_directory)"/analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.bed > S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed

bcftools view --threads 24 -O z -l 1 -T S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz -o gnomad.exomes.r2.1.1.sites.liftover_grch38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz

bcftools view --threads 24 -O z -l 1 -T S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed gnomad.genomes.r3.0.sites.vcf.bgz -o gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz

tabix -f -p vcf gnomad.exomes.r2.1.1.sites.liftover_grch38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz
tabix -f -p vcf gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz

#############
## SNPsift throughs an error due to empty description field, these commands corect it

zcat gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz | sed 's/Description=\"\"/Description=\"x\"/g' | bgzip -c > gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.corrected.vcf.gz

$ tabix -f -p vcf gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.corrected.vcf.gz