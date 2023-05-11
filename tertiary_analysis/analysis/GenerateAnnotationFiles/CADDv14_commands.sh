#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
$ cd analysis/annotation/hg19/cadd/

# download files
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/InDels.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/InDels.tsv.gz.tbi
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/gnomad.genomes.r2.0.1.sites.tsv.gz.tbi

#############
# merge, convert to VCF (awk, normalize to reference (bcftools)

(echo "##fileformat=VCFv4.2"; echo "##fileDate=2020-01-27"; echo "##reference=hg19"; echo "##CADD_version=1.5"; echo "##INFO=<ID=CADDv14_RawScore,Number=A,Type=Float,Description=\"Field 'RawScore' from CADD tsv files\">"; echo "##INFO=<ID=CADDv14_PHRED,Number=A,Type=Float,Description=\"Field 'PHRED' from CADD tsv files\">"; echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; zcat whole_genome_SNVs.tsv.gz InDels.tsv.gz gnomad.genomes.r2.0.1.sites.tsv.gz | grep -v "#" | awk 'BEGIN{OFS="\t";OFS="\t";} {print "chr"$1,$2,".",$3,$4,"100","PASS","CADDv14_RawScore="$5";CADDv14_PHRED="$6}') > CADDv14_all.hg19.canonicalchrom.vcf

bcftools sort -m 24000M -O v CADDv14_all.hg19.canonicalchrom.vcf -o CADDv14_all.hg19.vcf

bgzip -@ 24 CADDv14_all.hg19.vcf
tabix -f -p vcf CADDv14_all.hg19.vcf.gz

# clean up
rm CADDv14_all.hg19.canonicalchrom.vcf
rm gnomad.genomes.r2.0.1.sites.tsv.gz
rm gnomad.genomes.r2.0.1.sites.tsv.gz.tbi
rm InDels.tsv.gz.tbi
rm InDels.tsv.gz
rm whole_genome_SNVs.tsv.gz.tbi
rm whole_genome_SNVs.tsv.gz


#########################################
# subset to target for faster annotation:
sort -k 1,1 -k2,2n -k3,2n /catstor/user/bernt/GCKDpanel/analysis/targets/OID46243_hg19_02Jul2018_primary_targets.addedlines.plusminus50.bed > OID46243_hg19_02Jul2018_primary_targets.addedlines.plusminus50.sorted.bed

bcftools view --threads 24 -O z -l 1 -TOID46243_hg19_02Jul2018_primary_targets.addedlines.plusminus50.sorted.bed CADDv14_all.hg19.vcf.gz -o CADDv14_all.hg19.OID46243_hg19_02Jul2018_primary_targets.hg19.vcf.gz

tabix -f -p vcf CADDv14_all.hg19.OID46243_hg19_02Jul2018_primary_targets.hg19.vcf.gz