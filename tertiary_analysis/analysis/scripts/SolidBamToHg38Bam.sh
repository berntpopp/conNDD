#!/bin/sh
# SolidBamToHg38Bam.sh
# A shell script to align previously aligned BAM files from the SOLID platform to the hg38 reference using novoalignCS aligner
# Written by: Bernt Popp
# Last updated on: 2019-12-30
## example usage: `nohup sh SolidBamToHg38Bam.sh bam_input_filename bam_output_filename dna_id run_id platform_id &`
# -------------------------------------------------------
# Set Path to other scripts
PATH=analysis/scripts/:$PATH
# Set vars
bam_input_filename=$1																						# first command line argument equal to input BAM filename
bam_output_filename=$2																						# second command line argument equal to output BAM filename
dna_id=$3																									# third command line argument equal to DNA identifier
run_id=$4																									# fourth command line argument equal to run identifier
platform_id=$5																								# fifth command line argument equal to platform identifier
paired_bool=$6																								# sixth command line argument equal to boolean value indication whether the BAM file is paired or not
split_line_number=${7:-1000000}																				# optional argument to set the number of lines to split the file
input_dir=${8:-"data/seqs/BAMs"}																			# optional argument to set input directory
output_dir=${9:-"results/BAMs/hg38"}																		# optional argument to set output directory
nCS_reference=${10:-"analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.ncx"}								# optional argument to set reference
#########################################
## extract an unmapped fastq from aligned BAM
if [[ "$paired_bool" == "yes" ]]; then
	samtools sort -@ 10 -m 5G -n $input_dir/$bam_input_filename | samtools fastq -F 0x900 -T CS,CQ -1 $output_dir/$bam_input_filename.p1 -2 $output_dir/$bam_input_filename.p2 -0 /dev/null -s $output_dir/$bam_input_filename.s1 -n -
	parallel -j20 --delay 1 "bash SplitFastqNLines.sh {} $split_line_number" ::: $output_dir/$bam_input_filename.*
	parallel -j20 --delay 1 "bash ReplaceBaseWithColorSpace.sh {}" ::: $output_dir/$bam_input_filename.*-n*
	parallel -j20 --delay 20 "bash AlignPairedEndWithNovoalignCS.sh {1} {2} $dna_id $run_id $platform_id $nCS_reference; echo {1}.bam >> {1.}.merge.list " ::: $output_dir/$bam_input_filename.p1-n*-fq :::+ $output_dir/$bam_input_filename.p2-n*-fq
	parallel -j20 --delay 20 "bash AlignSingleEndWithNovoalignCS.sh {} $dna_id $run_id $platform_id $nCS_reference; echo {}.bam >> {.}.merge.list " ::: $output_dir/$bam_input_filename.s1-n*-fq
	samtools merge -@ 20 -c -p -b $output_dir/$bam_input_filename.merge.list $output_dir/$bam_output_filename
	samtools index $output_dir/$bam_output_filename
	parallel -j20 "rm {}" :::: $output_dir/$bam_input_filename.merge.list
elif [[ "$paired_bool" == "no" ]]; then
	samtools sort -@ 10 -m 5G -n $input_dir/$bam_input_filename | samtools fastq -F 0x900 -T CS,CQ -n - > $output_dir/$bam_input_filename.s1
	parallel -j20 --delay 1 "bash SplitFastqNLines.sh {} $split_line_number" ::: $output_dir/$bam_input_filename.*
	parallel -j20 --delay 1 "bash ReplaceBaseWithColorSpace.sh {}" ::: $output_dir/$bam_input_filename.*-n*
	parallel -j20 --delay 20 "bash AlignSingleEndWithNovoalignCS.sh {} $dna_id $run_id $platform_id $nCS_reference; echo {}.bam >> {.}.merge.list " ::: $output_dir/$bam_input_filename.s1-n*-fq
	samtools merge -@ 20 -c -p -b $output_dir/$bam_input_filename.merge.list $output_dir/$bam_output_filename
	samtools index $output_dir/$bam_output_filename
	parallel -j20 "rm {}" :::: $output_dir/$bam_input_filename.merge.list
fi
grep "." $output_dir/$bam_input_filename*.SplitFastqNLines.log $output_dir/$bam_input_filename*.ReplaceBaseWithColorSpace.log $output_dir/$bam_input_filename*.novoalignCS.log > $output_dir/_logs/$bam_input_filename.logs
rm $output_dir/$bam_input_filename*.SplitFastqNLines.log $output_dir/$bam_input_filename*.ReplaceBaseWithColorSpace.log $output_dir/$bam_input_filename*.novoalignCS.log