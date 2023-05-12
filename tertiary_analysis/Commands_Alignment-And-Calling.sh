#!/bin/sh

#########################################
# initialize anaconda
export PATH="$(config_get conda_path)"
#########################################


#########################################
# go to workdirektory
cd "$(config_get work_directory)"
#########################################


#########################################
# make project folders

mkdir data/ data/seqs/ data/seqs/BAMs/ data/seqs/BAMs/additional/ data/seqs/BAMs/extra/
mkdir analysis/ analysis/ref/ analysis/annotation/ analysis/scripts/ analysis/targets/ analysis/gatk/ analysis/genelists/
mkdir results/ results/VCFs/ results/VCFs/logs/ results/VCFs/annotation/ results/VCFs/annotation/stats/ results/VCFs/annotation/logs/ results/VCFs/filtering/ results/VCFs/filtering/logs/
mkdir results/BAMs/ results/BAMs/hg38/ results/BAMs/hg38/logs/ results/BAMs/hg38/additional/ results/BAMs/hg38/additional/logs/ results/BAMs/hg38/extra/ results/BAMs/hg38/extra/logs/
mkdir results/QC/ results/QC/qualimap/ results/QC/qualimap/mergedbed/ results/QC/qualimap/mergedbed/multibamqc/ results/QC/qualimap/logs/ results/QC/fastQC/ results/QC/fastQC/logs/
mkdir results/CNV/ results/CNV/CNVkit/ results/CNV/CNVkit/half-vs-half/ results/CNV/CNVkit/one-vs-all/ results/CNV/CNVkit/half-vs-half/calls/ results/CNV/CNVkit/half-vs-half/logs/ results/CNV/CNVkit/one-vs-all/calls/ results/CNV/CNVkit/one-vs-all/logs/
mkdir _tmp/

#########################################


#########################################
# Copy BAM data for all the runs, index them and calcuate checksums

scp -r "$(config_get bam_folder_storage_server)"{114951_Ampliseq_Exome,AmandaKing_103222.ontarget,AmandaKing_31261.ontarget,AmandaKing_42197.ontarget,AmandaKing_42236.ontarget,AmandaKing_51638.ontarget,AmandaKing_53175.ontarget,AmandaKing_76705.ontarget,BieneMinion_141020,BieneMinion_155541,BieneMinion_177896,BieneMinion_178223,BieneMinion_178226,Bonanza_53788.ontarget,Bonanza_53950.ontarget,Bonanza_53969.ontarget,Bonanza_57379.ontarget,Bonanza_57464.ontarget,Bonanza_57498.ontarget,Bonanza_58219.ontarget,Bonanza_62799.ontarget,Booker_SY005.ontarget,CaptainFuture_100097,CaptainFuture_143787,CaptainFuture_143812,CaptainFuture_143813,CaptainFuture_148007,CaptainFuture_151894,CaptainFuture_161929,CaptainFuture_173796,CaptainFuture_173797,CaptainFuture_173798,CaptainFuture_173799,CaptainFuture_173800,CaptainFuture_64910,Cops_132998.ontarget,DornRoesChen_127630,DornRoesChen_134073,DornRoesChen_134074,DornRoesChen_141019,DornRoesChen_145465,DrSnuggles_58440,DrSnuggles_79717,DukeBrothers_79717.ontarget,FastFurious_124939,FastFurious_148944,FastFurious_149202,FastFurious_53988,FastFurious_53992,FastFurious_56057,Fitz_30488.ontarget,Fitz_53854.ontarget,Fitz_54021.ontarget,Fitz_56081.ontarget,Fitz_58119.ontarget,Fitz_62657.ontarget,Fitz_75933.ontarget,FixFoxi_103222,FixFoxi_103247,FixFoxi_103248,FixFoxi_115318,FixFoxi_115319,FixFoxi_127470,FixFoxi_127637,FixFoxi_139414,FixFoxi_139415,FixFoxi_140223,FixFoxi_141520,FixFoxi_151833,FixFoxi_151834,FixFoxi_152050,FixFoxi_152051,FixFoxi_155540,FixFoxi_155541,FixFoxi_177839,FixFoxi_177852,FixFoxi_177859,FixFoxi_177896,FixFoxi_178221,FixFoxi_178223,FixFoxi_178225,FixFoxi_178226,FixFoxi_179763,FixFoxi_179785,FixFoxi_181379,FixFoxi_31261,FixFoxi_39447,FixFoxi_39588,FixFoxi_51638,FixFoxi_53162,FixFoxi_53163,FixFoxi_53777,FixFoxi_54004,FixFoxi_58676,FixFoxi_62829,FixFoxi_68969,FixFoxi_68970,FixFoxi_76922,Harry_SY217.ontarget,Higgins_53859.ontarget,hutch_SY155.ontarget,LolaRennt_151883,LolaRennt_151927,LolaRennt_79588,MikeHammer_113915.ontarget,MikeHammer_122565.ontarget,MikeHammer_125977.ontarget,MikeHammer_125984.ontarget,MikeHammer_79438.ontarget,MollyMill_30699.ontarget,MollyMill_57718.ontarget,MollyMill_72755.ontarget,Monk_SY005.ontarget,Paired_SY005.ontarget,PerryMason_113066.ontarget,PerryMason_114951.ontarget,PerryMason_53985.ontarget,PerryMason_72333.ontarget,rami2_53647.ontarget,rami2_53840.ontarget,rami2_53857.ontarget,rami_SY142.ontarget,rami_SY148.ontarget,rami_SY195.ontarget,rami_SY196.ontarget,RemingtonSteele_53801.ontarget,RemingtonSteele_53963.ontarget,RemingtonSteele_53979.ontarget,Rockford01_SY125.ontarget,Rockford01_SY138.ontarget,SherlockHolmes_114193.ontarget,SherlockHolmes_114194.ontarget,SherlockHolmes_52019.ontarget,SherlockHolmes_58226.ontarget,SherlockHolmes_58266.ontarget,SherlockHolmes_62408.ontarget,SpeedyGonzales_151758,SpeedyGonzales_151769,SpeedyGonzales_151914,SpeedyGonzales_152049,SpeedyGonzales_57425,TequilaBonetti_31758.ontarget,TequilaBonetti_32458.ontarget,TequilaBonetti_44462.ontarget,TequilaBonetti_51065.ontarget,TequilaBonetti_53816.ontarget,TequilaBonetti_54692.ontarget,TequilaBonetti_58441.ontarget,TequilaBonetti_75638.ontarget,TomUndJerry_115321,TomUndJerry_172395,TomUndJerry_172396,TomUndJerry_172397,TomUndJerry_53755,TomUndJerry_53783,TomUndJerry_53904,TomUndJerry_53972,TomUndJerry_54003,TomUndJerry_57375,TomUndJerry_57459,TomUndJerry_57506,TomUndJerry_58433,TomUndJerry_58434,TomUndJerry_58438,TomUndJerry_58444,TomUndJerry_69003}.bam data/seqs/BAMs/

# additional BAM files
scp -r "$(config_get bam_folder_storage_server)"{Fitz_53852.ontarget,MikeHammer_111393.ontarget,SpeedyGonzales_151843,SpeedyGonzales_151861,MRNet_62372,MRNet_62373,RemingtonSteele_53792.ontarget,PerryMason_53994.ontarget,hutch_SY308.ontarget}.bam data/seqs/BAMs/additional/

# extra BAM files
scp -r "$(config_get bam_folder_storage_server)"{KillBill_62430,FastFurious_124940,MikeHammer_124940.ontarget,Monk_SY085.ontarget,MRNet_62370,Paired_SY085.ontarget,PerryMason_100951.ontarget,KillBill_163731,rami2_53602.ontarget,rami_SY113.ontarget,rami_SY211.ontarget,RemingtonSteele_53790.ontarget,TomUndJerry_172398}.bam data/seqs/BAMs/extra/

# calcluate md5 checksums
parallel -j8 "md5sum data/seqs/BAMs/{/.}.bam >> data/seqs/BAMS.md5sum.txt" ::: data/seqs/BAMs/*.bam
parallel -j8 "md5sum data/seqs/BAMs/additional/{/.}.bam >> data/seqs/BAMS.md5sum.txt" ::: data/seqs/BAMs/additional/*.bam
parallel -j8 "md5sum data/seqs/BAMs/extra/{/.}.bam >> data/seqs/BAMS.md5sum.txt" ::: data/seqs/BAMs/extra/*.bam

# index BAM files
parallel -j8 "samtools index data/seqs/BAMs/{/.}.bam" ::: data/seqs/BAMs/*.bam
parallel -j8 "samtools index data/seqs/BAMs/additional/{/.}.bam" ::: data/seqs/BAMs/additional/*.bam
parallel -j8 "samtools index data/seqs/BAMs/extra/{/.}.bam" ::: data/seqs/BAMs/extra/*.bam

#########################################


#########################################
# download reference and index

wget ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -P analysis/ref/
wget ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt -P analysis/ref/
wget ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai -P analysis/ref/

picard CreateSequenceDictionary R=analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa O=analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict

bwa index analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

cat analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | grep -v "_" | grep -v "-" | grep -v "chrEBV" | cut -f 1 > analysis/ref/GRCh38.ID-list.txt

seqkit grep -f analysis/ref/GRCh38.ID-list.txt analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa > analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa

picard CreateSequenceDictionary R=analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa O=analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.dict

samtools faidx analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa

novoindex -c analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.ncx analysis/ref/GRCh38_full_analysis_set_no_decoy_hla.fa
#########################################


#########################################
# convert SOLID color-space BAMs to FASTQ and align to hg38 using novoalignCS

parallel --joblog _tmp/SolidBamToHg38Bam.log -j8 --delay 20 --header : --colsep '\t' bash analysis/scripts/SolidBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} {paired} :::: data/seqs/SOLID.BAMs.list

# additional BAM files
parallel --joblog _tmp/SolidBamToHg38Bam.additional.log -j4 --delay 20 --header : --colsep '\t' bash analysis/scripts/SolidBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} {paired} 1000000 data/seqs/BAMs/additional results/BAMs/hg38/additional :::: data/seqs/SOLID.BAMs-additional.list

# extra BAM files
parallel --joblog _tmp/SolidBamToHg38Bam.extra.log -j4 --delay 20 --header : --colsep '\t' bash analysis/scripts/SolidBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} {paired} 1000000 data/seqs/BAMs/extra results/BAMs/hg38/extra :::: data/seqs/SOLID.BAMs-extra.list
#########################################


#########################################
# revert Illumina BAMs and align to hg38 using bwa-mem

parallel --joblog _tmp/IlluminaBamToHg38Bam.log -j14 --header : --colsep '\t' sh analysis/scripts/IlluminaBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} :::: data/seqs/ILLUMINA.BAMs.list

# additional BAM files
parallel --joblog _tmp/IlluminaBamToHg38Bam.additional.log -j4 --header : --colsep '\t' sh analysis/scripts/IlluminaBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} data/seqs/BAMs/additional results/BAMs/hg38/additional :::: data/seqs/ILLUMINA.BAMs-additional.list

# extra BAM files
parallel --joblog _tmp/IlluminaBamToHg38Bam.extra.log -j5 --header : --colsep '\t' sh analysis/scripts/IlluminaBamToHg38Bam.sh {FilenameBAM} {FilenameRealignedBAM} {DNAid} {Run} {PL_attribute} data/seqs/BAMs/extra results/BAMs/hg38/extra :::: data/seqs/ILLUMINA.BAMs-extra.list

#########################################


#########################################
# mark duplicates with picard

mkdir results/BAMs/hg38/MarkDups/ results/BAMs/hg38/MarkDups/logs/

parallel --joblog _tmp/ILLUMINA_MarkDuplicatesWithMateCigar.log -j40 --delay 20 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/ILLUMINA.BAMs.list

parallel --joblog _tmp/SOLID_MarkDuplicatesWithMateCigar.log -j70 --delay 20 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/SOLID.BAMs.list

nohup parallel -j20 --delay 10 'samtools index -@ 12 {}' ::: results/BAMs/hg38/MarkDups/*.bam &
parallel -j8 "md5sum results/BAMs/hg38/MarkDups/{/.}.bam >> results/BAMs/hg38/MarkDups/BAMS.md5sum.txt" ::: results/BAMs/hg38/MarkDups/*.bam


# additional BAM files
mkdir results/BAMs/hg38/additional/MarkDups/ results/BAMs/hg38/additional/MarkDups/logs/

parallel --joblog _tmp/ILLUMINA_MarkDuplicatesWithMateCigar-additional.log -j10 --delay 10 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/additional/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/additional/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/additional/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/ILLUMINA.BAMs-additional.list

parallel --joblog _tmp/SOLID_MarkDuplicatesWithMateCigar-additional.log -j10 --delay 10 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/additional/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/additional/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/additional/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/SOLID.BAMs-additional.list

nohup parallel -j20 --delay 10 'samtools index -@ 12 {}' ::: results/BAMs/hg38/additional/MarkDups/*.bam &
parallel -j4 "md5sum results/BAMs/hg38/additional/MarkDups/{/.}.bam >> results/BAMs/hg38/additional/MarkDups/BAMS.md5sum.txt" ::: results/BAMs/hg38/additional/MarkDups/*.bam


# extra BAM files
mkdir results/BAMs/hg38/extra/MarkDups/ results/BAMs/hg38/extra/MarkDups/logs/

parallel --joblog _tmp/ILLUMINA_MarkDuplicatesWithMateCigar-extra.log -j10 --delay 10 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/extra/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/extra/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/extra/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/ILLUMINA.BAMs-extra.list

parallel --joblog _tmp/SOLID_MarkDuplicatesWithMateCigar-extra.log -j10 --delay 10 --header : --colsep '\t' 'picard MarkDuplicatesWithMateCigar I=results/BAMs/hg38/extra/{FilenameRealignedBAM.}.bam O=results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam M=results/BAMs/hg38/extra/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.txt 2>results/BAMs/hg38/extra/MarkDups/logs/{FilenameRealignedBAM.}.MarkDuplicatesWithMateCigar.log' :::: data/seqs/SOLID.BAMs-extra.list

nohup parallel -j20 --delay 10 'samtools index -@ 12 {}' ::: results/BAMs/hg38/extra/MarkDups/*.bam &
parallel -j4 "md5sum results/BAMs/hg38/extra/MarkDups/{/.}.bam >> results/BAMs/hg38/extra/MarkDups/BAMS.md5sum.txt" ::: results/BAMs/hg38/extra/MarkDups/*.bam
#########################################



#########################################
# quality control on BAM files using `qualimap_v2.2.2c` (installed on 2019-11-07 `conda install -c bioconda qualimap`)

parallel --joblog _tmp/ILLUMINA_qualimap.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/mergedbed/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/ILLUMINA.BAMs.list

parallel --joblog _tmp/SOLID_qualimap.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/mergedbed/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/SOLID.BAMs.list

cat data/seqs/ILLUMINA.BAMs.list data/seqs/SOLID.BAMs.list | grep -v "^FilenameBAM" | awk 'BEGIN{FS="\t"; OFS="\t";} {print $3,"results/QC/qualimap/mergedbed/"gensub(".bam", "", "g", $13),$2}' > results/QC/qualimap/conNDD.multibamqc.list
qualimap multi-bamqc -d results/QC/qualimap/conNDD.multibamqc.list -outdir results/QC/qualimap/mergedbed/multibamqc/


# additional BAM files
mkdir results/QC/qualimap/additional/ results/QC/qualimap/additional/multibamqc/

parallel --joblog _tmp/ILLUMINA_qualimap-additional.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/additional/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/ILLUMINA.BAMs-additional.list

parallel --joblog _tmp/SOLID_qualimap-additional.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/additional/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/SOLID.BAMs-additional.list


# extra BAM files
mkdir results/QC/qualimap/extra/ results/QC/qualimap/extra/multibamqc/

parallel --joblog _tmp/ILLUMINA_qualimap-extra.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/extra/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/ILLUMINA.BAMs-extra.list

parallel --joblog _tmp/SOLID_qualimap-extra.log -j20 --delay 40 --header : --colsep '\t' 'qualimap bamqc --java-mem-size=8G -c -nt 8 -outdir results/QC/qualimap/extra/{FilenameRealignedBAM.}/ -outformat PDF:HTML -gff analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed -bam results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam >& results/QC/qualimap/logs/{FilenameRealignedBAM.}.MarkDups.bam.qualimap.log' :::: data/seqs/SOLID.BAMs-extra.list


# grep quality summary information
grep ">= 10X" results/QC/qualimap/mergedbed/*/genome_results.txt | sed 's/results\/QC\/qualimap\/mergedbed\///g' | sed 's/\/genome_results.txt\://g' | sed 's/There is a //g' | sed 's/ of reference with a coverageData >= 10X//g' | sed 's/\.nCS\.hg38/\.nCS\.hg38\.MarkDups\.bam/g' | sed 's/\.bwa\.hg38/\.bwa\.hg38\.MarkDups\.bam/g' > results/QC/qualimap/conNDD.multibamqc.10X.txt

grep ">= 10X" results/QC/qualimap/extra/*/genome_results.txt | sed 's/results\/QC\/qualimap\/extra\///g' | sed 's/\/genome_results.txt\://g' | sed 's/There is a //g' | sed 's/ of reference with a coverageData >= 10X//g' | sed 's/\.nCS\.hg38/\.nCS\.hg38\.MarkDups\.bam/g' | sed 's/\.bwa\.hg38/\.bwa\.hg38\.MarkDups\.bam/g' > results/QC/qualimap/conNDD.multibamqc.extra.10X.txt

#########################################



#########################################
# quality control on BAM files using `fastQC`

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/ILLUMINA.BAMs.list
parallel --joblog _tmp/ILLUMINA_fastqc.log -j20 --delay 20 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/ILLUMINA.BAMs.list

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/SOLID.BAMs.list
parallel --joblog _tmp/SOLID_fastqc.log -j20 --delay 20 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/SOLID.BAMs.list


# additional BAM files
mkdir results/QC/fastQC/additional/ results/QC/fastQC/additional/logs/

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/additional/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/ILLUMINA.BAMs-additional.list
parallel --joblog _tmp/ILLUMINA_fastqc-additional.log -j20 --delay 10 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/additional/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/additional/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/ILLUMINA.BAMs-additional.list

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/additional/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/SOLID.BAMs-additional.list
parallel --joblog _tmp/SOLID_fastqc-additional.log -j20 --delay 10 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/additional/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/additional/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/additional/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/SOLID.BAMs-additional.list


# extra BAM files
mkdir results/QC/fastQC/extra/ results/QC/fastQC/extra/logs/

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/extra/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/ILLUMINA.BAMs-extra.list
parallel --joblog _tmp/ILLUMINA_fastqc-extra.log -j20 --delay 10 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/extra/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/extra/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/ILLUMINA.BAMs-extra.list

parallel -j20 --header : --colsep '\t' 'mkdir results/QC/fastQC/extra/{FilenameRealignedBAM.}.MarkDups.bam/' :::: data/seqs/SOLID.BAMs-extra.list
parallel --joblog _tmp/SOLID_fastqc-extra.log -j20 --delay 10 --header : --colsep '\t' 'fastqc -t 8 results/BAMs/hg38/extra/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam -o results/QC/fastQC/extra/{FilenameRealignedBAM.}.MarkDups.bam/ 2> results/QC/fastQC/extra/logs/{FilenameRealignedBAM.}.MarkDups.bam.fastQC.log' :::: data/seqs/SOLID.BAMs-extra.list

#########################################



#########################################
# call variants from exomes with GATK HaplotyeCaller

parallel --joblog _tmp/ILLUMINA_gatkHC.log -j40 --delay 20 --header : --colsep '\t' 'gatk HaplotypeCaller \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-ERC GVCF \
-I results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam \
-O results/VCFs/{FilenameRealignedBAM.}.MarkDups.g.vcf 2> results/VCFs/logs/{FilenameRealignedBAM.}.MarkDups.gatkHC.log' :::: data/seqs/ILLUMINA.BAMs.list

parallel --joblog _tmp/SOLID_gatkHC.log -j40 --delay 20 --header : --colsep '\t' 'gatk HaplotypeCaller \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-ERC GVCF \
-I results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam \
-O results/VCFs/{FilenameRealignedBAM.}.MarkDups.g.vcf 2> results/VCFs/logs/{FilenameRealignedBAM.}.MarkDups.gatkHC.log' :::: data/seqs/SOLID.BAMs.list
#########################################



#########################################
# GenotypeGVCFs

##$ cat data/seqs/ILLUMINA.BAMs.list data/seqs/SOLID.BAMs.list | grep -v "^FilenameBAM" | awk 'BEGIN{FS="\t"; OFS="\t";} {print $3,"results/VCFs/"gensub(".bam", ".MarkDups.g.vcf", "g", $13)}' > results/VCFs/conNDDcohort.sample_map
## had to make this in excel

## downloaded scattered calling intervals on 2020-01-08
gsutil cp -r gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ analysis/targets/
cat analysis/targets/scattered_calling_intervals/*/scattered.interval_list | grep -v "^@" > analysis/targets/resources_broad_hg38.scattered_calling_intervals.merged.bed

mkdir results/VCFs/conNDDdatabase/

parallel --joblog _tmp/GenomicsDBImport.log -j15 --delay 20 --colsep '\t' 'gatk GenomicsDBImport --genomicsdb-workspace-path results/VCFs/conNDDdatabase/conNDDdatabase_{1}-{2}-{3}/ --sample-name-map results/VCFs/conNDDcohort.sample_map -L {1}:{2}-{3} 2> results/VCFs/conNDDdatabase/conNDDdatabase_{1}-{2}-{3}.GenomicsDBImport.log' :::: analysis/targets/resources_broad_hg38.scattered_calling_intervals.merged.bed

mkdir results/VCFs/JointCallCohort/ results/VCFs/JointCallCohort/logs/

parallel -j80 --delay 5 --colsep '\t' 'gatk GenotypeGVCFs \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-V gendb://results/VCFs/conNDDdatabase/conNDDdatabase_{1}-{2}-{3}/ \
-G StandardAnnotation \
-L {1}:{2}-{3} \
-O results/VCFs/JointCallCohort/conNDDcohort.hc-joint.{1}-{2}-{3}.vcf 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.{1}-{2}-{3}.log' :::: analysis/targets/resources_broad_hg38.scattered_calling_intervals.merged.bed

parallel -j1 --colsep '\t' ' echo results/VCFs/JointCallCohort/conNDDcohort.hc-joint.{1}-{2}-{3}.vcf ' :::: analysis/targets/resources_broad_hg38.scattered_calling_intervals.merged.bed > results/VCFs/JointCallCohort/conNDDcohort.JointCallCohort.vcf.list

picard MergeVcfs \
I=results/VCFs/JointCallCohort/conNDDcohort.JointCallCohort.vcf.list \
D=analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict \
O=results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.MergeVcf.vcf.gz.log

tabix -f -p vcf results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz
#########################################



#########################################
# recalibrate variants scores using VQSR
# based on: https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-fc.wdl
# based on: https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.hg38.wgs.inputs.json
# based on: https://gatk.broadinstitute.org/hc/en-us/articles/360036727711-VariantRecalibrator)
# then filter variants and prepare for annotation

# Subset to SNPs-only callset with SelectVariants
# for all variants
gatk SelectVariants \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz \
--select-type-to-include SNP \
-O results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.snps.vcf.gz 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.MergeVcf.vcf.gz.SelectVariants-SNP.log

# only on target
gatk SelectVariants \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz \
--select-type-to-include SNP \
--intervals analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.bed \
-O results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.vcf.gz 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.MergeVcf.vcf.gz.SelectVariants-SNP-ontarget.log

# Subset to indels-only callset with SelectVariants
# for all variants
gatk SelectVariants \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz \
--select-type-to-exclude SNP \
-O results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.indels.vcf.gz 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.MergeVcf.vcf.gz.SelectVariants-indels.log

# only on target
gatk SelectVariants \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.vcf.gz \
--select-type-to-exclude SNP \
--intervals analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.bed \
-O results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.vcf.gz 2> results/VCFs/JointCallCohort/logs/conNDDcohort.hc-joint.MergeVcf.vcf.gz.SelectVariants-indels-ontarget.log

mkdir results/VCFs/VariantRecalibrator/ results/VCFs/VariantRecalibrator/logs/ results/VCFs/VariantRecalibrator/model/

# download resoure bundle files
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz analysis/gatk/1000G_omni2.5.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi analysis/gatk/1000G_omni2.5.hg38.vcf.gz.tbi
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz analysis/gatk/hapmap_3.3.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi analysis/gatk/hapmap_3.3.hg38.vcf.gz.tbi
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz analysis/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi analysis/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf analysis/gatk/Homo_sapiens_assembly38.dbsnp138.vcf
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx analysis/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.idx
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz analysis/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi analysis/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz analysis/gatk/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi analysis/gatk/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi

# VariantRecalibrator for SNPs
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.vcf.gz \
-O results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recall.vcf.gz \
--tranches-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.VQSR.snps-ontarget.vcf.gz.tranches \
--trust-all-polymorphic \
-tranche 100.00 -tranche 99.95 -tranche 99.90 -tranche 99.80 -tranche 99.60 -tranche 99.50 -tranche 99.40 -tranche 99.30 -tranche 99.00 -tranche 98.00 -tranche 97.00 -tranche 90.00 \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
-mode SNP \
--output-model results/VCFs/VariantRecalibrator/model/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.vcf.gz.model \
--max-gaussians 6 \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 analysis/gatk/hapmap_3.3.hg38.vcf.gz \
--resource:omni,known=false,training=true,truth=false,prior=12.0 analysis/gatk/1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 analysis/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 analysis/gatk/Homo_sapiens_assembly38.dbsnp138.vcf \
--rscript-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recall.vcf.gz.output.plots.R 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recall.log

# VariantRecalibrator for INDELs
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.vcf.gz \
-O results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recall.vcf.gz \
--tranches-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.VQSR.indels-ontarget.vcf.gz.tranches \
--trust-all-polymorphic \
-tranche 100.00 -tranche 99.95 -tranche 99.90 -tranche 99.50 -tranche 99.00 -tranche 97.00 -tranche 96.00 -tranche 95.00 -tranche 94.00 -tranche 93.40 -tranche 93.00 -tranche 92.00 -tranche 91.00 -tranche 90.00 \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
-mode INDEL \
--max-gaussians 4 \
--resource:mills,known=false,training=true,truth=true,prior=12 analysis/gatk/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10 analysis/gatk/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2 analysis/gatk/Homo_sapiens_assembly38.dbsnp138.vcf \
--rscript-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recall.vcf.gz.output.plots.R 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recall.log

# Apply calibration for SNPs
gatk --java-options "-Xmx10g -Xms10g" ApplyVQSR \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.vcf.gz \
-O results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recalibrated.vcf.gz \
--recal-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recall.vcf.gz \
--tranches-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.VQSR.snps-ontarget.vcf.gz.tranches \
--truth-sensitivity-filter-level 99.50 \
--create-output-variant-index true \
-mode SNP 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.ApplyVQSR.log

# Apply calibration for INDELS
gatk --java-options "-Xmx10g -Xms10g" ApplyVQSR \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-V results/VCFs/JointCallCohort/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.vcf.gz \
-O results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recalibrated.vcf.gz \
--recal-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recall.vcf.gz \
--tranches-file results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.VQSR.indels-ontarget.vcf.gz.tranches \
--truth-sensitivity-filter-level 99.50 \
--create-output-variant-index true \
-mode INDEL 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.ApplyVQSR.log

# join the filtered snps and indels
picard MergeVcfs \
I=results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.snps-ontarget.recalibrated.vcf.gz \
I=results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.indels-ontarget.recalibrated.vcf.gz \
D=analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.dict \
O=results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.vcf.gz 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.vcf.gz.MergeVcfs.log

# normalize indels and split multiallelic 
gatk LeftAlignAndTrimVariants \
-R analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-V results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.vcf.gz \
-O results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.vcf.gz \
--split-multi-allelics 2> results/VCFs/VariantRecalibrator/logs/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.vcf.gz.LeftAlignAndTrimVariants.log

tabix -f -p vcf results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.vcf.gz
#########################################


#########################################
# annotation of DNA variants

snpEff -Xmx16g hg38 -lof -noInteraction -noMotif -noNextProt -spliceRegionIntronMax 12 results/VCFs/VariantRecalibrator/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.vcf.gz -stats results/VCFs/annotation/stats/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.vcf.gz.html 2> results/VCFs/annotation/logs/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.vcf.gz.snpEff.log > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.vcf

SnpSift -Xmx4g dbnsfp -db analysis/annotation/hg38/dbNSFP/dbNSFP4.0a.txt.gz -n -f dbNSFP_pos_1_coor_,genename results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.vcf > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/spidex/spidex.hg38.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.vcf -info spidex_dpsi_max_tissue,spidex_dpsi_zscore > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/dbscSNV/dbscSNV.hg38.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.vcf -info dbscSNV_ada_score,dbscSNV_rf_score > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/cadd/CADDv15_all.hg38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.vcf -info CADDv15_PHRED > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.vcf -info nhomalt,AC,AN,AF -name gADe211_ > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/gnomAD/gnomad.genomes.r3.0.sites.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.corrected.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.vcf -info nhomalt,AC,AN,AF -name gADg30_ > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/bravo/bravo-dbsnp-all.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.vcf -info Hom,AC,AN,AF -name bravo_ > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/clinvar/clinvar_20200127.hg38.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.vcf -info CLNSIG -name clinvar_20200127_ > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.clinvar.vcf

SnpSift -Xmx4g annotate analysis/annotation/hg38/hgmd/hgmd.2019-3.hg38.vcf.gz results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.clinvar.vcf -info variant_type -name hgmd_20193_ > results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.clinvar.hgmd.vcf

zcat results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.clinvar.hgmd.vcf.gz | grep "#CHROM" > results/VCFs/annotation/conNDDcohort.header.txt

parallel -j8 bgzip -@ 24 {} ::: results/VCFs/annotation/*.vcf
parallel -j8 tabix -f -p vcf {} ::: results/VCFs/annotation/*.vcf.gz

zcat results/VCFs/annotation/conNDDcohort.hc-joint.MergeVcf.ontarget.recalibrated.split.ann.dbNSFP40a.spidex.dbscSNV.CADDv15.gADe211.gADg30.bravo.clinvar.hgmd.vcf.gz | sed 's/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t100097\t103222\t103247\t103248\t103282\t108121\t113066\t113915\t114193\t114194\t114951\t115318\t115319\t115321\t122565\t124939\t125977\t125984\t127470\t127630\t127637\t132998\t134073\t134074\t134296\t139414\t139415\t140223\t141019\t141020\t141520\t143787\t143812\t143813\t145465\t147640\t148007\t148090\t148944\t149202\t149206\t150817\t151079\t151084\t151151\t151161\t151162\t151165\t151173\t151192\t151733\t151736\t151739\t151758\t151769\t151788\t151831\t151833\t151834\t151883\t151885\t151894\t151914\t151927\t151942\t152049\t152050\t152051\t155540\t155541\t161929\t162117\t163722\t163733\t163741\t163742\t172395\t172396\t172397\t173796\t173797\t173798\t173799\t173800\t177839\t177852\t177859\t177896\t178221\t178223\t178225\t178226\t179763\t179785\t181379\t30488\t30699\t31261\t31758\t32458\t39446\t39447\t39588\t42197\t42236\t44462\t51065\t51638\t52019\t53162\t53163\t53647\t53755\t53769\t53777\t53783\t53788\t53801\t53808\t53816\t53840\t53854\t53857\t53859\t53904\t53950\t53963\t53969\t53972\t53979\t53985\t53988\t53992\t54003\t54004\t54021\t54093\t54692\t56057\t56081\t57375\t57379\t57425\t57459\t57464\t57498\t57506\t57718\t58119\t58219\t58226\t58266\t58302\t58433\t58434\t58438\t58440\t58441\t58444\t58676\t62408\t62425\t62657\t62799\t62828\t62829\t62920\t64910\t68969\t68970\t69003\t72333\t72755\t75638\t75933\t76922\t79438\t79588\t79717\tSY005\tSY125\tSY138\tSY142\tSY148\tSY155\tSY195\tSY196\tSY217/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts100097\ts103222\ts103247\ts103248\ts103282\ts108121\ts113066\ts113915\ts114193\ts114194\ts114951\ts115318\ts115319\ts115321\ts122565\ts124939\ts125977\ts125984\ts127470\ts127630\ts127637\ts132998\ts134073\ts134074\ts134296\ts139414\ts139415\ts140223\ts141019\ts141020\ts141520\ts143787\ts143812\ts143813\ts145465\ts147640\ts148007\ts148090\ts148944\ts149202\ts149206\ts150817\ts151079\ts151084\ts151151\ts151161\ts151162\ts151165\ts151173\ts151192\ts151733\ts151736\ts151739\ts151758\ts151769\ts151788\ts151831\ts151833\ts151834\ts151883\ts151885\ts151894\ts151914\ts151927\ts151942\ts152049\ts152050\ts152051\ts155540\ts155541\ts161929\ts162117\ts163722\ts163733\ts163741\ts163742\ts172395\ts172396\ts172397\ts173796\ts173797\ts173798\ts173799\ts173800\ts177839\ts177852\ts177859\ts177896\ts178221\ts178223\ts178225\ts178226\ts179763\ts179785\ts181379\ts30488\ts30699\ts31261\ts31758\ts32458\ts39446\ts39447\ts39588\ts42197\ts42236\ts44462\ts51065\ts51638\ts52019\ts53162\ts53163\ts53647\ts53755\ts53769\ts53777\ts53783\ts53788\ts53801\ts53808\ts53816\ts53840\ts53854\ts53857\ts53859\ts53904\ts53950\ts53963\ts53969\ts53972\ts53979\ts53985\ts53988\ts53992\ts54003\ts54004\ts54021\ts54093\ts54692\ts56057\ts56081\ts57375\ts57379\ts57425\ts57459\ts57464\ts57498\ts57506\ts57718\ts58119\ts58219\ts58226\ts58266\ts58302\ts58433\ts58434\ts58438\ts58440\ts58441\ts58444\ts58676\ts62408\ts62425\ts62657\ts62799\ts62828\ts62829\ts62920\ts64910\ts68969\ts68970\ts69003\ts72333\ts72755\ts75638\ts75933\ts76922\ts79438\ts79588\ts79717\tsSY005\tsSY125\tsSY138\tsSY142\tsSY148\tsSY155\tsSY195\tsSY196\tsSY217/g' | bgzip -@ 24 > results/VCFs/annotation/conNDDcohort.annoated.vcf.gz

tabix -f -p vcf results/VCFs/annotation/conNDDcohort.annoated.vcf.gz

#########################################


#########################################
# VCF filtering commands
# using custom scripts

###############
# homozygous
parallel --joblog _tmp/FilterHomozygous_single.log -j80 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterHomozygous_single.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA NA {analysis} {project_FamilyID} {IndividualID} NA' :::: results/VCFs/filtering/filtering.single.list :::: results/VCFs/filtering/analysis.list

parallel --joblog _tmp/FilterHomozygous_trio.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterHomozygous_trio.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} s{Mother_DNA_Identifier} s{Father_DNA_Identifier} NA {analysis} {project_FamilyID} {IndividualID} NA' :::: results/VCFs/filtering/filtering.trio.list :::: results/VCFs/filtering/analysis.list

parallel --joblog _tmp/FilterHomozygous_duo.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterHomozygous_duo.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA s{Other_DNA_Identifier} {analysis} {project_FamilyID} {IndividualID} NA' :::: results/VCFs/filtering/filtering.duo.list :::: results/VCFs/filtering/analysis.list

###############
# dominant
parallel --joblog _tmp/FilterDominant_single.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterDominant_single.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA NA {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.single.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhADaXaOther.list

parallel --joblog _tmp/FilterDominant_trio.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterDominant_trio.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} s{Mother_DNA_Identifier} s{Father_DNA_Identifier} NA {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.trio.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhADaXaOther.list

parallel --joblog _tmp/FilterDominant_duo.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterDominant_duo.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA s{Other_DNA_Identifier} {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.duo.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhADaXaOther.list

###############
# recessive
parallel --joblog _tmp/FilterRecessive_single.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterRecessive_single.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA NA {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.single.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhAR.list

parallel --joblog _tmp/FilterRecessive_trio.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterRecessive_trio.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} s{Mother_DNA_Identifier} s{Father_DNA_Identifier} NA {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.trio.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhAR.list

parallel --joblog _tmp/FilterRecessive_duo.log -j30 --delay 1 --header : --colsep '\t' 'bash analysis/scripts/filter/FilterRecessive_single.sh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz results/VCFs/filtering s{DNA_Identifier} NA NA s{Other_DNA_Identifier} {analysis} {project_FamilyID} {IndividualID} {genelist}' :::: results/VCFs/filtering/filtering.duo.list :::: results/VCFs/filtering/analysis.list :::: analysis/genelists/genelists_inhAR.list

###############
# tar and zip for upload and processing in R
tar -zcvf filtering_2020-02-07.tar.gz results/VCFs/filtering/

#########################################



#########################################
# run CNVkit for copy number detection

###############
# one-vs-all analysis

# SOLID
parallel -j8 --header : --colsep ';' 'cnvkit.py batch {File} --normal {=Control_Files uq(); =} --drop-low-coverage -p 20 -m amplicon \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/one-vs-all/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}.cnb \
--output-dir results/CNV/CNVkit/one-vs-all/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}/ \
--diagram --scatter 2> results/CNV/CNVkit/one-vs-all/logs/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}.log' :::: results/CNV/CNVkit/CNVkit.one-vs-all.worklist_SOLID.2020-03-29.csv

cnvkit.py heatmap $(ls results/CNV/CNVkit/one-vs-all/*SOLID/*.cns) -x f -d -o results/CNV/CNVkit/one-vs-all/conNDD.SOLID.CNVkit.heatmap.pdf

# ILLUMINA
parallel -j8 --header : --colsep ';' 'cnvkit.py batch {File} --normal {=Control_Files uq(); =} --drop-low-coverage -p 20 -m hybrid \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/one-vs-all/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}.cnb \
--output-dir results/CNV/CNVkit/one-vs-all/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}/ \
--diagram --scatter 2> results/CNV/CNVkit/one-vs-all/logs/conNDD.{Cluster}.{Sample}.{Group}.{PL_attribute}.log' :::: results/CNV/CNVkit/CNVkit.one-vs-all.worklist_ILLUMINA.2020-03-29.csv

cnvkit.py heatmap $(ls results/CNV/CNVkit/one-vs-all/*ILLUMINA/*.cns) -x f -d -o results/CNV/CNVkit/one-vs-all/conNDD.ILLUMINA.CNVkit.heatmap.pdf

# call and aggregate CNVkit results for Excel
parallel -j48 'cnvkit.py call -m threshold -t=-2.0,-0.41,0.32,0.80,1.1 {} -o results/CNV/CNVkit/one-vs-all/calls/{/.}.call.cns' ::: results/CNV/CNVkit/one-vs-all/*/*.cns

grep '.*' results/CNV/CNVkit/one-vs-all/calls/*.call.cns | grep -v "chromosome" | sed 's/results\/CNV\/CNVkit\/one-vs-all\/calls\///g' | sed 's/\.call\.cns\:/\t/g' | sed 's/\.bwa\.hg38\.MarkDups//g' | sed 's/\.nCS\.hg38\.MarkDups//g' | awk 'BEGIN {OFS="\t"; FS="\t"}; {if($7!=2){split($1,a,"_"); print(a[1],a[2],a[3],a[4],$2,$3,$4,$5,$6,$7,$8,$9,$10)}}' > results/CNV/CNVkit/one-vs-all/calls/conNDD.CNVkit.exome-aberrations.txt

###############


###############
## half-vs-half analysis

###############
# general functions and files preparation
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# generate accessable region file for hg38
cnvkit.py access -o analysis/targets/access-5k-mappable.hg38.bed analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# download refflat file for hg38 from UCSC for annotation
wget -O analysis/targets/refFlat.hg38.txt.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz
gunzip analysis/targets/refFlat.hg38.txt.gz

#########
# SOLID
parallel --header : --colsep '\t' 'echo results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam' :::: data/seqs/SOLID.BAMs.list | shuf --random-source=<(get_seeded_random 42) > results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list
split --lines=34 results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list.

nohup cnvkit.py batch \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list.aa | tr '\n' ' ') \
--normal \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list.ab | tr '\n' ' ') \
--drop-low-coverage \
-p 80 \
-m amplicon \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/half-vs-half/conNDD.SOLID.CNVkit.cohort.ab.cnb \
--output-dir results/CNV/CNVkit/half-vs-half/conNDD.SOLID.CNVkit.cohort.aa/ \
--diagram --scatter 2> results/CNV/CNVkit/half-vs-half/logs/conNDD.SOLID.CNVkit.cohort.aa.log &

nohup cnvkit.py batch \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list.ab | tr '\n' ' ') \
--normal \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.SOLID.CNVkit.cohort.list.aa | tr '\n' ' ') \
--drop-low-coverage \
-p 80 \
-m amplicon \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/half-vs-half/conNDD.SOLID.CNVkit.cohort.aa.cnb \
--output-dir results/CNV/CNVkit/half-vs-half/conNDD.SOLID.CNVkit.cohort.ab/ \
--diagram --scatter 2> results/CNV/CNVkit/half-vs-half/logs/conNDD.SOLID.CNVkit.cohort.ab.log &

cnvkit.py heatmap $(ls results/CNV/CNVkit/half-vs-half/*SOLID*/*.cns) -x f -d -o results/CNV/CNVkit/half-vs-half/conNDD.SOLID.CNVkit.heatmap.pdf

#########
# ILLUMINA
parallel --header : --colsep '\t' 'echo results/BAMs/hg38/MarkDups/{FilenameRealignedBAM.}.MarkDups.bam' :::: data/seqs/ILLUMINA.BAMs.list | shuf --random-source=<(get_seeded_random 42) > results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list
split --lines=67 results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list.

nohup cnvkit.py batch \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list.aa | tr '\n' ' ') \
--normal \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list.ab | tr '\n' ' ') \
--drop-low-coverage \
-p 80 \
-m hybrid \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--access analysis/targets/access-5k-mappable.hg38.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/half-vs-half/conNDD.ILLUMINA.CNVkit.cohort.ab.cnb \
--output-dir results/CNV/CNVkit/half-vs-half/conNDD.ILLUMINA.CNVkit.cohort.aa/ \
--diagram --scatter 2> results/CNV/CNVkit/half-vs-half/logs/conNDD.ILLUMINA.CNVkit.cohort.aa.log &

nohup cnvkit.py batch \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list.ab | tr '\n' ' ') \
--normal \
$(cat results/CNV/CNVkit/conNDD.half-vs-half.ILLUMINA.CNVkit.cohort.list.aa | tr '\n' ' ') \
--drop-low-coverage \
-p 80 \
-m hybrid \
--targets analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.bed \
--access analysis/targets/access-5k-mappable.hg38.bed \
--fasta analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--annotate analysis/targets/refFlat.hg38.txt \
--output-reference results/CNV/CNVkit/half-vs-half/conNDD.ILLUMINA.CNVkit.cohort.aa.cnb \
--output-dir results/CNV/CNVkit/half-vs-half/conNDD.ILLUMINA.CNVkit.cohort.ab/ \
--diagram --scatter 2> results/CNV/CNVkit/half-vs-half/logs/conNDD.ILLUMINA.CNVkit.cohort.ab.log &

cnvkit.py heatmap $(ls results/CNV/CNVkit/half-vs-half/*ILLUMINA*/*.cns) -x f -d -o results/CNV/CNVkit/half-vs-half/conNDD.ILLUMINA.CNVkit.heatmap.pdf

# call and aggregate CNVkit results for Excel
parallel -j48 'cnvkit.py call -m threshold -t=-2.0,-0.41,0.32,0.80,1.1 {} -o results/CNV/CNVkit/half-vs-half/calls/{/.}.call.cns' ::: results/CNV/CNVkit/half-vs-half/*/*.cns

grep '.*' results/CNV/CNVkit/half-vs-half/calls/*.call.cns | grep -v "chromosome" | sed 's/results\/CNV\/CNVkit\/calls\///g' | sed 's/\.call\.cns\:/\t/g' | sed 's/\.bwa\.hg38\.MarkDups//g' | sed 's/\.nCS\.hg38\.MarkDups//g' | awk 'BEGIN {OFS="\t"; FS="\t"}; {if($7!=2){split($1,a,"_"); print(a[1],a[2],a[3],a[4],$2,$3,$4,$5,$6,$7,$8,$9,$10)}}' > results/CNV/CNVkit/half-vs-half/calls/conNDD.CNVkit.exome-aberrations.txt
######################################