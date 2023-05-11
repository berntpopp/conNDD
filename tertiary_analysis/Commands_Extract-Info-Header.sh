#!/bin/sh
## 2020-01-27:


#########################################
# go to workdirektory
cd "$(config_get work_directory)"
#########################################

# extract RG tag from BAMs
parallel 'samtools view -H {} | grep "^@RG" | sed "s/^@RG/{/}/g"' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.RG.raw.txt

# extract PG tag from BAMs
parallel 'samtools view -H {} | grep "^@PG" | sed "s/^@PG/{/}/g"' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.PG.raw.txt

# extract SQ tag from BAMs
parallel 'samtools view -H {} | grep "^@SQ" | sed "s/^@SQ/{/}/g"' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.SQ.raw.txt

# extract PairedReadCount from BAM
parallel -j 12 'printf {/}"\t"; samtools view -c -f 1 {} | tr "\n" "\t"; printf "\n"' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.PairedReadCount.txt

# extract ReadLengthMax from BAM
parallel -j 12 'printf {/}"\t"; samtools view {} | head -100000 | cut -f 10 | awk "{print length($1)}" | sort -n -r | head -1' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.ReadLengthMax.txt

# extract ReadLengthMin from BAM
parallel -j 12 'printf {/}"\t"; samtools view {} | head -100000 | cut -f 10 | awk "{print length($1)}" | sort -n | head -1' ::: data/seqs/BAMs/*/*.bam > data/seqs/BAMS-additional-and-extra.ReadLengthMin.txt