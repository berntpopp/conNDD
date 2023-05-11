#!/bin/sh
# IlluminaBamToHg38Bam.sh
# A shell script to align previously aligned BAM files from the Illumina platform to the hg38 reference using bwa-mem aligner
# Written by: Bernt Popp
# Last updated on: 2019-12-25
## example usage: `nohup sh IlluminaBamToHg38Bam.sh bam_input_filename bam_output_filename dna_id run_id platform_id &`
# -------------------------------------------------------
# Set vars
bam_input_filename=$1																						# first command line argument equal to input BAM filename
bam_output_filename=$2																						# second command line argument equal to output BAM filename
dna_id=$3																									# third command line argument equal to DNA identifier
run_id=$4																									# fourth command line argument equal to run identifier
platform_id=$5																								# fifth command line argument equal to platform identifier
input_dir=${6:-"data/seqs/BAMs"}																			# optional argument to set input directory
output_dir=${7:-"results/BAMs/hg38"}																		# optional argument to set output directory
bwa_reference=${8:-"analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"}								# optional argument to set reference
#########################################
## Generate an unmapped BAM from aligned BAM (https://gatkforums.broadinstitute.org/gatk/discussion/6484)
picard -Xmx8G RevertSam \
I=$input_dir/$bam_input_filename \
O=$output_dir/$bam_input_filename.revertsam \
SANITIZE=true \
MAX_DISCARD_FRACTION=0.005 \
ATTRIBUTE_TO_CLEAR=XT \
ATTRIBUTE_TO_CLEAR=XN \
ATTRIBUTE_TO_CLEAR=AS \
ATTRIBUTE_TO_CLEAR=OC \
ATTRIBUTE_TO_CLEAR=OP \
SORT_ORDER=queryname \
RESTORE_ORIGINAL_QUALITIES=true \
REMOVE_DUPLICATE_INFORMATION=true \
REMOVE_ALIGNMENT_INFORMATION=true 2>$output_dir/logs/$bam_input_filename.RevertSam.log
#########################################
## Mark adapter sequences using MarkIlluminaAdapters (https://software.broadinstitute.org/gatk/documentation/article?id=6483#step2)
picard -Xmx8G MarkIlluminaAdapters \
I=$output_dir/$bam_input_filename.revertsam \
O=$output_dir/$bam_input_filename.markilluminaadapters \
M=$output_dir/logs/$bam_input_filename.markilluminaadapters_metrics.txt \
TMP_DIR=_tmp/ 2>$output_dir/logs/$bam_input_filename.MarkIlluminaAdapters.log
#########################################
## Pipe SamToFastq, BWA-MEM and MergeBamAlignment to generate a clean BAM (https://software.broadinstitute.org/gatk/documentation/article?id=6483#step3D)
## bwa commands and reference building also based on https://www.biorxiv.org/content/10.1101/868570v1.supplementary-material
picard -Xmx8G SamToFastq \
I=$output_dir/$bam_input_filename.markilluminaadapters \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=_tmp/ 2>$output_dir/logs/$bam_input_filename.SamToFastq.log | \
bwa mem -M -t 12 -p $bwa_reference /dev/stdin 2>$output_dir/logs/$bam_output_filename.bwa-mem.log | \
picard -Xmx16G MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=$output_dir/$bam_input_filename.revertsam \
OUTPUT=$output_dir/$bam_output_filename \
R=$bwa_reference CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=_tmp/ 2>$output_dir/logs/$bam_output_filename.MergeBamAlignment.log
#########################################
# Remove intermediate files
rm $output_dir/$bam_input_filename.revertsam
rm $output_dir/$bam_input_filename.markilluminaadapters