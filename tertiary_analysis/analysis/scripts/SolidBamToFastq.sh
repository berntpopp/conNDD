#!/bin/sh
# SolidBamToFastq.sh
# A shell script to convert aligned color-space BAM files from the SOLID platform to gziped FASTQ files
# Written by: Bernt Popp
# Last updated on: 2019-12-23
## example usage: `nohup sh SolidBamToFastq.sh bam_filename fq_filename &`
# -------------------------------------------------------
# Set vars
bam_filename=$1											# first command line argument equal to input BAM filename
fq_filename=$2											# second command line argument equal to output FASSTQ filename
input_dir=${3:-"data/seqs/BAMs"}						# optional argument to set input directory
output_dir=${4:-"data/seqs/FQs"}						# optional argument to set output directory
samtools view -u -F 2048 -F 256 $input_dir/$bam_filename | samtools sort -@ 4 -m 5G -n | samtools fastq -T CS,CQ - | paste - - - - | awk -v OFS="\n" -v ORS="\n" '{gsub("CQ:Z:","",$0); gsub("CS:Z:","",$0); {print $1,$2,$5,$3}}' | gzip > $output_dir/$fq_filename