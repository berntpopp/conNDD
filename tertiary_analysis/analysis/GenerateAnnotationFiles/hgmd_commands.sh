#!/bin/sh

#############
## 2020-02-11

## Genome_Trax_HGMD_2019.3.tar.gz was copied from locally

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
cd analysis/annotation/hg19/hgmd/

tar xvzf Genome_Trax_HGMD_2019.3.tar.gz

(grep "^##" hgmd-hg19.vcf | sed 's/Description=""/Description="na"/g'; \
echo "##contig=<ID=chr1>"; echo "##contig=<ID=chr2>"; echo "##contig=<ID=chr3>"; echo "##contig=<ID=chr4>"; echo "##contig=<ID=chr5>"; echo "##contig=<ID=chr6>"; echo "##contig=<ID=chr7>"; echo "##contig=<ID=chr8>"; echo "##contig=<ID=chr9>"; echo "##contig=<ID=chr10>"; echo "##contig=<ID=chr11>"; echo "##contig=<ID=chr12>"; echo "##contig=<ID=chr13>"; echo "##contig=<ID=chr14>"; echo "##contig=<ID=chr15>"; echo "##contig=<ID=chr16>"; echo "##contig=<ID=chr17>"; echo "##contig=<ID=chr18>"; echo "##contig=<ID=chr19>"; echo "##contig=<ID=chr20>"; echo "##contig=<ID=chr21>"; echo "##contig=<ID=chr22>"; echo "##contig=<ID=chrX>"; echo "##contig=<ID=chrY>"; \
grep "^#" hgmd-hg19.vcf | grep -v "^##"; \
grep -v "^#" hgmd-hg19.vcf | sed 's/^/chr/g') > hgmd.2019-3.hg19.changed.vcf

# normalize indels and split multiallelic 
gatk LeftAlignAndTrimVariants \
-R /storage/genomes/hg19/hg19.fa \
-V hgmd.2019-3.hg19.changed.vcf \
-O hgmd.2019-3.hg19.changed.split.vcf \
--max-indel-length 5000 \
--max-leading-bases 5000 \
--split-multi-allelics

bcftools sort -O v hgmd.2019-3.hg19.changed.split.vcf -o hgmd.2019-3.hg19.vcf

bgzip hgmd.2019-3.hg19.vcf
tabix -f -p vcf hgmd.2019-3.hg19.vcf.gz

# clean up
rm -rf *.vcf *.bed *.gff