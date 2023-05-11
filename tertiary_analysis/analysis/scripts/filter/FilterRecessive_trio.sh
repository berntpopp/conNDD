#!/bin/sh
# FilterHomozygous.sh
# A shell script to filter the cohort BAM for homozygous variants
# Written by: Bernt Popp
# Last updated on: 2020-02-02
## example usage: `bash FilterRecessive_single.sh vcf_file output_dir gene_list sample_id mother_id father_id other_id analysis familyid individualid`
# -------------------------------------------------------
# Set vars
vcf_file=$1																							# 1. command line argument equal to input VCF filename
output_dir=$2																						# 2. command line argument equal to output directory
sample_id=$3																						# 3. command line argument equal to DNA of the sample as in the VCF
mother_id=$4																						# 4. command line argument equal to DNA of the mother as in the VCF
father_id=$5																						# 5. command line argument equal to DNA of the father as in the VCF
other_id=$6																							# 6. command line argument equal to DNA of the other family individual as in the VCF
analysis=$7																							# 7. argument to set the type of filtering analysis
family_id=$8																						# 8. argument to set the family identifier
individual_id=$9																					# 9. argument to set the individual identifier
gene_list=${10}																						# 10. command line argument equal to genelist to use for filtering
#########################################
## start filtering
if [[ "$analysis" == "A45" ]]; then
	zcat $vcf_file | \
	SnpSift -Xmx8g filter --set analysis/genelists/$gene_list"_2020-01-11.txt" \
	" (ANN[ANY].GENE in SET[0]) & isVariant(GEN[$sample_id]) & \
	( (AF < 0.02) & ((QUAL >= 250) | (FILTER = 'PASS')) & (GEN[$sample_id].DP[0] >= 4) & \
	(! (clinvar_20200127_CLNSIG = 'Benign') & ! (clinvar_20200127_CLNSIG = 'Likely_benign') & ! (clinvar_20200127_CLNSIG = 'Benign/Likely_benign')) & \
	((gADe211_AF < 0.01) | (na gADe211_AC)) & \
	((gADe211_nhomalt = 0) | (na gADe211_nhomalt)) & \
	((gADg30_AF < 0.01) | (na gADg30_AC)) & \
	((gADg30_nhomalt = 0) | (na gADg30_nhomalt)) & \
	((bravo_AF < 0.01) | (na bravo_AC)) & \
	((bravo_Hom = 0) | (na bravo_Hom)) ) & \
	((clinvar_20200127_CLNSIG = 'Likely_pathogenic') | (clinvar_20200127_CLNSIG = 'Pathogenic') | (hgmd_20193_variant_type = 'DM')) " | \
	SnpSift -Xmx8g extractFields -s "," -e "na" - \
	CHROM POS ID REF ALT dbNSFP_hg19_chr dbNSFP_hg19_pos_1_based_ FILTER QUAL AC AF AN ANN[0].GENE ANN[0].FEATUREID ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P \
	gADe211_AC[0] gADe211_AN[0] gADe211_AF[0] \
	gADg30_AC[0] gADg30_AN[0] gADg30_AF[0] \
	bravo_AC[0] bravo_AN[0] bravo_AF[0] \
	CADDv15_PHRED[0] dbNSFP_MetaLR_score[0] dbNSFP_MetaSVM_score[0] dbNSFP_M_CAP_score[0] dbNSFP_REVEL_score[0] dbNSFP_MVP_score[0] dbNSFP_MPC_score[0] dbNSFP_PrimateAI_score[0] \
	dbNSFP_SIFT_score[0] dbNSFP_Polyphen2_HVAR_score[0] dbNSFP_MutationTaster_score[0] \
	spidex_dpsi_zscore[0] dbscSNV_ada_score[0] dbscSNV_rf_score[0] clinvar_20200127_CLNSIG[0] hgmd_20193_variant_type[0] \
	GEN[$sample_id].GT GEN[$sample_id].AD[0] GEN[$sample_id].DP[0] GEN[$mother_id].GT GEN[$father_id].GT | awk -F "\t" "BEGIN {OFS = FS} {gsub("/$sample_id/",\"index\",\$0); gsub("/$mother_id/",\"mother\",\$0); gsub("/$father_id/",\"father\",\$0); print \$0,\"na\"}" | sed 's/GEN\[father\].GT\tna/GEN\[father\]\.GT\tGEN\[other\]\.GT/g' 2>$output_dir/logs/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".log" > $output_dir/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".txt"
elif [[ "$analysis" == "LGD" ]]; then
	zcat $vcf_file | \
	SnpSift -Xmx8g filter --set analysis/genelists/$gene_list"_2020-01-11.txt" \
	" (ANN[ANY].GENE in SET[0]) & isVariant(GEN[$sample_id]) & \
	( (AF < 0.02) & ((QUAL >= 250) | (FILTER = 'PASS')) & (GEN[$sample_id].DP[0] >= 4) & \
	(! (clinvar_20200127_CLNSIG = 'Benign') & ! (clinvar_20200127_CLNSIG = 'Likely_benign') & ! (clinvar_20200127_CLNSIG = 'Benign/Likely_benign')) & \
	((gADe211_AF < 0.01) | (na gADe211_AC)) & \
	((gADe211_nhomalt = 0) | (na gADe211_nhomalt)) & \
	((gADg30_AF < 0.01) | (na gADg30_AC)) & \
	((gADg30_nhomalt = 0) | (na gADg30_nhomalt)) & \
	((bravo_AF < 0.01) | (na bravo_AC)) & \
	((bravo_Hom = 0) | (na bravo_Hom)) ) & \
	((ANN[ANY].IMPACT has 'HIGH') | (LOF[ANY].PERC > 0.5) | (NMD[ANY].PERC > 0.5) | (ANN[ANY].EFFECT =~ 'inframe') | (ANN[ANY].EFFECT =~ 'stop') | (ANN[ANY].EFFECT =~ 'initiator')) " | \
	SnpSift -Xmx8g extractFields -s "," -e "na" - \
	CHROM POS ID REF ALT dbNSFP_hg19_chr dbNSFP_hg19_pos_1_based_ FILTER QUAL AC AF AN ANN[0].GENE ANN[0].FEATUREID ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P \
	gADe211_AC[0] gADe211_AN[0] gADe211_AF[0] \
	gADg30_AC[0] gADg30_AN[0] gADg30_AF[0] \
	bravo_AC[0] bravo_AN[0] bravo_AF[0] \
	CADDv15_PHRED[0] dbNSFP_MetaLR_score[0] dbNSFP_MetaSVM_score[0] dbNSFP_M_CAP_score[0] dbNSFP_REVEL_score[0] dbNSFP_MVP_score[0] dbNSFP_MPC_score[0] dbNSFP_PrimateAI_score[0] \
	dbNSFP_SIFT_score[0] dbNSFP_Polyphen2_HVAR_score[0] dbNSFP_MutationTaster_score[0] \
	spidex_dpsi_zscore[0] dbscSNV_ada_score[0] dbscSNV_rf_score[0] clinvar_20200127_CLNSIG[0] hgmd_20193_variant_type[0] \
	GEN[$sample_id].GT GEN[$sample_id].AD[0] GEN[$sample_id].DP[0] GEN[$mother_id].GT GEN[$father_id].GT | awk -F "\t" "BEGIN {OFS = FS} {gsub("/$sample_id/",\"index\",\$0); gsub("/$mother_id/",\"mother\",\$0); gsub("/$father_id/",\"father\",\$0); print \$0,\"na\"}" | sed 's/GEN\[father\].GT\tna/GEN\[father\]\.GT\tGEN\[other\]\.GT/g' 2>$output_dir/logs/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".log" > $output_dir/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".txt"
elif [[ "$analysis" == "Missense" ]]; then
	zcat $vcf_file | \
	SnpSift -Xmx8g filter --set analysis/genelists/$gene_list"_2020-01-11.txt" \
	" (ANN[ANY].GENE in SET[0]) & isVariant(GEN[$sample_id]) & \
	( (AF < 0.02) & ((QUAL >= 250) | (FILTER = 'PASS')) & (GEN[$sample_id].DP[0] >= 4) & \
	(! (clinvar_20200127_CLNSIG = 'Benign') & ! (clinvar_20200127_CLNSIG = 'Likely_benign') & ! (clinvar_20200127_CLNSIG = 'Benign/Likely_benign')) & \
	((gADe211_AF < 0.01) | (na gADe211_AC)) & \
	((gADe211_nhomalt = 0) | (na gADe211_nhomalt)) & \
	((gADg30_AF < 0.01) | (na gADg30_AC)) & \
	((gADg30_nhomalt = 0) | (na gADg30_nhomalt)) & \
	((bravo_AF < 0.01) | (na bravo_AC)) & \
	((bravo_Hom = 0) | (na bravo_Hom)) ) & \
	((CADDv15_PHRED >= 15)) & ((ANN[*] has 'missense_variant')) " | \
	SnpSift -Xmx8g extractFields -s "," -e "na" - \
	CHROM POS ID REF ALT dbNSFP_hg19_chr dbNSFP_hg19_pos_1_based_ FILTER QUAL AC AF AN ANN[0].GENE ANN[0].FEATUREID ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P \
	gADe211_AC[0] gADe211_AN[0] gADe211_AF[0] \
	gADg30_AC[0] gADg30_AN[0] gADg30_AF[0] \
	bravo_AC[0] bravo_AN[0] bravo_AF[0] \
	CADDv15_PHRED[0] dbNSFP_MetaLR_score[0] dbNSFP_MetaSVM_score[0] dbNSFP_M_CAP_score[0] dbNSFP_REVEL_score[0] dbNSFP_MVP_score[0] dbNSFP_MPC_score[0] dbNSFP_PrimateAI_score[0] \
	dbNSFP_SIFT_score[0] dbNSFP_Polyphen2_HVAR_score[0] dbNSFP_MutationTaster_score[0] \
	spidex_dpsi_zscore[0] dbscSNV_ada_score[0] dbscSNV_rf_score[0] clinvar_20200127_CLNSIG[0] hgmd_20193_variant_type[0] \
	GEN[$sample_id].GT GEN[$sample_id].AD[0] GEN[$sample_id].DP[0] GEN[$mother_id].GT GEN[$father_id].GT | awk -F "\t" "BEGIN {OFS = FS} {gsub("/$sample_id/",\"index\",\$0); gsub("/$mother_id/",\"mother\",\$0); gsub("/$father_id/",\"father\",\$0); print \$0,\"na\"}" | sed 's/GEN\[father\].GT\tna/GEN\[father\]\.GT\tGEN\[other\]\.GT/g' 2>$output_dir/logs/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".log" > $output_dir/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".txt"
elif [[ "$analysis" == "Splice" ]]; then
	zcat $vcf_file | \
	SnpSift -Xmx8g -Xmx8g filter -n " ((ANN[ANY].IMPACT has 'HIGH')) " | \
	SnpSift -Xmx8g filter --set analysis/genelists/$gene_list"_2020-01-11.txt" \
	" (ANN[ANY].GENE in SET[0]) & isVariant(GEN[$sample_id]) & \
	( (AF < 0.02) & ((QUAL >= 250) | (FILTER = 'PASS')) & (GEN[$sample_id].DP[0] >= 4) & \
	(! (clinvar_20200127_CLNSIG = 'Benign') & ! (clinvar_20200127_CLNSIG = 'Likely_benign') & ! (clinvar_20200127_CLNSIG = 'Benign/Likely_benign')) & \
	((gADe211_AF < 0.01) | (na gADe211_AC)) & \
	((gADe211_nhomalt = 0) | (na gADe211_nhomalt)) & \
	((gADg30_AF < 0.01) | (na gADg30_AC)) & \
	((gADg30_nhomalt = 0) | (na gADg30_nhomalt)) & \
	((bravo_AF < 0.01) | (na bravo_AC)) & \
	((bravo_Hom = 0) | (na bravo_Hom)) ) & \
	( (CADDv15_PHRED >= 15) & (( spidex_dpsi_zscore[0] <= -2.0 ) | ( dbscSNV_ada_score[0] >= 0.6 ) | ( dbscSNV_rf_score[0] >= 0.6 )) ) " | \
	SnpSift -Xmx8g extractFields -s "," -e "na" - \
	CHROM POS ID REF ALT dbNSFP_hg19_chr dbNSFP_hg19_pos_1_based_ FILTER QUAL AC AF AN ANN[0].GENE ANN[0].FEATUREID ANN[0].EFFECT ANN[0].IMPACT ANN[0].HGVS_C ANN[0].HGVS_P \
	gADe211_AC[0] gADe211_AN[0] gADe211_AF[0] \
	gADg30_AC[0] gADg30_AN[0] gADg30_AF[0] \
	bravo_AC[0] bravo_AN[0] bravo_AF[0] \
	CADDv15_PHRED[0] dbNSFP_MetaLR_score[0] dbNSFP_MetaSVM_score[0] dbNSFP_M_CAP_score[0] dbNSFP_REVEL_score[0] dbNSFP_MVP_score[0] dbNSFP_MPC_score[0] dbNSFP_PrimateAI_score[0] \
	dbNSFP_SIFT_score[0] dbNSFP_Polyphen2_HVAR_score[0] dbNSFP_MutationTaster_score[0] \
	spidex_dpsi_zscore[0] dbscSNV_ada_score[0] dbscSNV_rf_score[0] clinvar_20200127_CLNSIG[0] hgmd_20193_variant_type[0] \
	GEN[$sample_id].GT GEN[$sample_id].AD[0] GEN[$sample_id].DP[0] GEN[$mother_id].GT GEN[$father_id].GT | awk -F "\t" "BEGIN {OFS = FS} {gsub("/$sample_id/",\"index\",\$0); gsub("/$mother_id/",\"mother\",\$0); gsub("/$father_id/",\"father\",\$0); print \$0,\"na\"}" | sed 's/GEN\[father\].GT\tna/GEN\[father\]\.GT\tGEN\[other\]\.GT/g' 2>$output_dir/logs/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".log" > $output_dir/$family_id"_"$individual_id"_"$sample_id".recessive_"$analysis"_"$gene_list".txt"
fi