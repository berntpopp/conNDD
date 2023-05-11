#!/bin/sh
# SplitFastqNLines.sh
# A shell script to convert FASTQ files into tab delimited files and split them into multiple with N lines
# Written by: Bernt Popp
# Last updated on: 2019-12-30
## example usage: `nohup sh SplitFastqNLines.sh input_filename line_number &`
# -------------------------------------------------------
# Set vars
input_filename=$1																					# first command line argument equal to input filename
line_number=$2																						# second command line argument equal number of lines
#########################################
## read the file, convert to 4 column tab and split into N lines
cat $input_filename | paste - - - - | awk "NR%$line_number==1{x=\"$input_filename\"\"-n\"++i;}{print > x}"
wc -l $input_filename > $input_filename.SplitFastqNLines.log
rm $input_filename