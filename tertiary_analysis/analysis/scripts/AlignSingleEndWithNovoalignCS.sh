#!/bin/sh
# AlignSingleEndWithNovoalignCS.sh
# A shell script to align single-end color-space FASTQ files using novoalignCS
# Written by: Bernt Popp
# Last updated on: 2019-12-30
## example usage: `nohup sh AlignSingleEndWithNovoalignCS.sh input_filename_1 dna_id run_id platform_id nCS_reference &`
## make execuatable chmod a+x AlignSingleEndWithNovoalignCS.sh
# -------------------------------------------------------
# Set vars
input_filename_1=$1																							# first command line argument equal to input BAM filename
dna_id=$2																									# second command line argument equal to DNA identifier
run_id=$3																									# third command line argument equal to run identifier
platform_id=$4																								# fourth command line argument equal to platform identifier
nCS_reference=$5																							# fifth command line argument argument to set reference
#########################################
## perform alignment
novoalignCS -r Random -o SAM "@RG\tID:$dna_id\tPL:$platform_id\tSM:$dna_id\tLB:$run_id-$dna_id" -d $nCS_reference -f $input_filename_1 2>$input_filename_1.novoalignCS.log | samtools view -Sb - | samtools sort -@ 4 -m 5G - -o $input_filename_1.bam
rm $input_filename_1