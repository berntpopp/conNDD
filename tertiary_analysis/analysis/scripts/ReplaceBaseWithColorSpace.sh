#!/bin/sh
# ReplaceBaseWithColorSpace.sh
# A shell script to convert tab delimited files with base and color space information into color space fastq files
# Written by: Bernt Popp
# Last updated on: 2019-12-30
## example usage: `nohup sh ReplaceBaseWithColorSpace.sh input_filename &`
# -------------------------------------------------------
# Set vars
input_filename=$1																					# first command line argument equal to input filename
#########################################
## read the file, convert to 4 column tab and split into N lines
cat $input_filename | awk -v OFS="\n" -v ORS="\n" '{gsub("CQ:Z:","",$0); gsub("CS:Z:","",$0); {print $1,$2,$5,$3}}' > $input_filename-fq
wc -l $input_filename > $input_filename.ReplaceBaseWithColorSpace.log
rm $input_filename