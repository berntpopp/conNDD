###########################################################
## load required libraries
# install.packages(c("tidyverse", "fs","fuzzyjoin"))
library(tidyverse)
library(fs)
library(fuzzyjoin)
library(readxl)
library(config)			## used for config loading
###########################################################

###########################################################
## define work directory and excel file locations

##set the work directoty from config file
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir))
###########################################################


###########################################################
## define data dir
data_dir <- "filtering/extracted/"


###########################################################
## load other data for gen annotation

dbNSFP4_0_gene <- read_delim("filtering/data/dbNSFP4.0_gene.gz", "\t", escape_double = FALSE, na = ".", trim_ws = TRUE, col_types = cols(ExAC_pLI = col_double(), ExAC_pRec = col_double(), `P(rec)` = col_double(), SORVA_LOForMissense_MAF0.001_HomOrCompoundHet = col_double()))

gevir_gene_ranking <- read_csv("filtering/data/gene_ranking.csv")

# the file FileS2_conNDD-cohort.xlsx needs to be downloaded in current version from Zenodo
FileS2_location <- "data/FileS2_conNDD-cohort.xlsx"

FileS2_individuals <- read_excel(FileS2_location, sheet = "individuals", na = "NA")
FileS2_BAMs <- read_excel(FileS2_location, sheet = "BAMs", na = "NA")

# read homozygous regions from Excel
HomozygosityRegions_2020_02_06 <- read_excel("data/HomozygosityRegions_2020-02-03.xlsx", sheet = "hg38_grouped", na = "NA")


###########################################################
## transform data
dbNSFP4_and_gevir <- dbNSFP4_0_gene %>%
dplyr::left_join(gevir_gene_ranking, by = c("Ensembl_gene" = "gene_id")) %>%
dplyr::select(Gene_name, MIM_phenotype_id, ExAC_pLI, ExAC_pRec, `P(HI)`, `P(rec)`, Known_rec_info, SORVA_LOF_MAF0.005_HetOrHom, SORVA_LOForMissense_MAF0.005_HomOrCompoundHet, gevir_ad_enrichment, loeuf_ad_enrichment, virlof_ad_enrichment, gevir_ar_enrichment, loeuf_ar_enrichment, virlof_ar_enrichment) %>%
separate(MIM_phenotype_id, into = c("MIM_phenotype_id"), sep = ";", remove = TRUE) %>%
mutate(ExAC_pLI = case_when( !is.na(ExAC_pLI) ~ round(ExAC_pLI,4))) %>%
mutate(ExAC_pRec = case_when( !is.na(ExAC_pRec) ~ round(ExAC_pRec,4))) %>%
mutate(`P(HI)` = case_when( !is.na(`P(HI)`) ~ round(`P(HI)`,4))) %>%
mutate(`P(rec)` = case_when( !is.na(`P(rec)`) ~ round(`P(rec)`,4))) %>%
mutate(SORVA_LOF_MAF0.005_HetOrHom = case_when( !is.na(SORVA_LOF_MAF0.005_HetOrHom) ~ round(SORVA_LOF_MAF0.005_HetOrHom,4))) %>%
mutate(SORVA_LOForMissense_MAF0.005_HomOrCompoundHet = case_when( !is.na(SORVA_LOForMissense_MAF0.005_HomOrCompoundHet) ~ round(SORVA_LOForMissense_MAF0.005_HomOrCompoundHet,4))) %>%
mutate(gevir_ad_enrichment = case_when( !is.na(gevir_ad_enrichment) ~ round(gevir_ad_enrichment,4))) %>%
mutate(loeuf_ad_enrichment = case_when( !is.na(loeuf_ad_enrichment) ~ round(loeuf_ad_enrichment,4))) %>%
mutate(virlof_ad_enrichment = case_when( !is.na(virlof_ad_enrichment) ~ round(virlof_ad_enrichment,4))) %>%
mutate(gevir_ar_enrichment = case_when( !is.na(gevir_ar_enrichment) ~ round(gevir_ar_enrichment,4))) %>%
mutate(loeuf_ar_enrichment = case_when( !is.na(loeuf_ar_enrichment) ~ round(loeuf_ar_enrichment,4))) %>%
mutate(virlof_ar_enrichment = case_when( !is.na(virlof_ar_enrichment) ~ round(virlof_ar_enrichment,4)))

Supplementary_Information_JAMA_Psychiatry_2020_02_06 <- FileS2_individuals %>%
dplyr::filter(grepl("^yes", Is_rjp_individual) & Type == "Index") %>%
dplyr::left_join(FileS2_BAMs, by = c("DNA_Identifier" = "DNAid")) %>%
dplyr::filter(primary == "yes", Is_rjp_BAMs == "yes") %>%
dplyr::select(project_FamilyID, IndividualID, DNA_Identifier, Run, PL_attribute)

###########################################################
## load the filenames
txt_files_homozygous <- fs::dir_ls(data_dir, regexp = "homozygous")
txt_files_homozygous

txt_files_dominant <- fs::dir_ls(data_dir, regexp = "dominant")
txt_files_dominant

txt_files_recessive <- fs::dir_ls(data_dir, regexp = "recessive")
txt_files_recessive

###########################################################
## load and merge the file content
homozygous_all <- txt_files_homozygous %>% 
map_dfr(read_tsv, .id = "source", na = c("na"), col_types = "cdccccccddddccccccdddddddddddddddddddddddcccddccc")

dominant_all <- txt_files_dominant %>% 
map_dfr(read_tsv, .id = "source", na = c("na"), col_types = "cdccccccddddccccccdddddddddddddddddddddddcccddccc")

recessive_all <- txt_files_recessive %>% 
map_dfr(read_tsv, .id = "source", na = c("na"), col_types = "cdccccccddddccccccdddddddddddddddddddddddcccddccc")

###########################################################
## process and arrange data

#############
## homozygous

conNDD_homozygous <- homozygous_all %>%
mutate(source = str_remove(source, "filtering/extracted/")) %>%
mutate(source = str_remove(source, ".txt")) %>%
separate(source, into = c("File", "Analysis"), sep = "\\.", remove = FALSE) %>%
mutate(File = gsub("_s", "_", File, perl = TRUE)) %>%
separate(File, into = c("project_FamilyID", "IndividualID", "DNA_Identifier"), sep = "_") %>%
dplyr::rename(CHROM_hg19 = `dbNSFP_hg19_chr`, POS_hg19 = `dbNSFP_hg19_pos_1_based_`, GENE = `ANN[0].GENE`, FEATUREID = `ANN[0].FEATUREID`, EFFECT = `ANN[0].EFFECT`, IMPACT = `ANN[0].IMPACT`, HGVS_C = `ANN[0].HGVS_C`, HGVS_P = `ANN[0].HGVS_P`, gADe211_AC = `gADe211_AC[0]`, gADe211_AN = `gADe211_AN[0]`, gADe211_AF = `gADe211_AF[0]`, gADg30_AC = `gADg30_AC[0]`, gADg30_AN = `gADg30_AN[0]`, gADg30_AF = `gADg30_AF[0]`, bravo_AC = `bravo_AC[0]`, bravo_AN = `bravo_AN[0]`, bravo_AF = `bravo_AF[0]`, CADDv15_PHRED = `CADDv15_PHRED[0]`, MetaLR = `dbNSFP_MetaLR_score[0]`, MetaSVM = `dbNSFP_MetaSVM_score[0]`, M_CAP = `dbNSFP_M_CAP_score[0]`, REVEL = `dbNSFP_REVEL_score[0]`, MVP = `dbNSFP_MVP_score[0]`, MPC = `dbNSFP_MPC_score[0]`, PrimateAI = `dbNSFP_PrimateAI_score[0]`, SIFT = `dbNSFP_SIFT_score[0]`, Polyphen2_HVAR = `dbNSFP_Polyphen2_HVAR_score[0]`, MutationTaster = `dbNSFP_MutationTaster_score[0]`, spidex_dpsi = `spidex_dpsi_zscore[0]`, dbscSNV_ada = `dbscSNV_ada_score[0]`, dbscSNV_rf = `dbscSNV_rf_score[0]`, clinvar_20200127 = `clinvar_20200127_CLNSIG[0]`, hgmd_20193 = `hgmd_20193_variant_type[0]`, index_GT = `GEN[index].GT`, index_AD = `GEN[index].AD[0]`, index_DP = `GEN[index].DP[0]`, mother_GT = `GEN[mother].GT`, father_GT = `GEN[father].GT`, other_GT = `GEN[other].GT`) %>%
mutate(frac_ADDP = round(1-(index_AD/index_DP),8)) %>%
mutate(p_ADDP = round(pbinom(index_AD,index_DP,0.5),8)) %>%
dplyr::select(project_FamilyID, IndividualID, DNA_Identifier, Analysis, CHROM, POS, ID, REF, ALT, CHROM_hg19, POS_hg19, FILTER, QUAL, AC, AF, AN, index_GT, index_AD, index_DP, frac_ADDP, p_ADDP, mother_GT, father_GT, other_GT, GENE, FEATUREID, EFFECT, IMPACT, HGVS_C, HGVS_P, gADe211_AC, gADe211_AN, gADe211_AF, gADg30_AC, gADg30_AN, gADg30_AF, bravo_AC, bravo_AN, bravo_AF, CADDv15_PHRED, MetaLR, MetaSVM, M_CAP, REVEL, MVP, MPC, PrimateAI, SIFT, Polyphen2_HVAR, MutationTaster, spidex_dpsi, dbscSNV_ada, dbscSNV_rf, clinvar_20200127, hgmd_20193) %>%
dplyr::distinct() %>%
group_by(.dots=c("project_FamilyID", "IndividualID", "DNA_Identifier", "CHROM", "POS", "REF", "ALT")) %>%
mutate(Analysis = paste(Analysis, collapse = ',')) %>%
dplyr::distinct() %>% ungroup %>%
mutate(ScoresCount = ifelse(ifelse(is.na(MetaLR),0,MetaLR) >= 0.5,1,0)+ifelse(ifelse(is.na(MetaSVM),0,MetaSVM) > 0,1,0)+ifelse(ifelse(is.na(M_CAP),0,M_CAP) >= 0.025,1,0)+ifelse(ifelse(is.na(REVEL),0,REVEL) >= 0.5,1,0)+ifelse(ifelse(is.na(MVP),0,MVP) >= 0.15,1,0)+ifelse(ifelse(is.na(MPC),0,MPC) >= 2,1,0)+ifelse(ifelse(is.na(PrimateAI),0,PrimateAI) >= 0.803,1,0)) %>%
separate(POS_hg19, into = c("POS_hg19"), sep = ",", remove = TRUE) %>%
separate(CHROM_hg19, into = c("CHROM_hg19"), sep = ",", remove = TRUE) %>%
mutate(CHROM_hg19 = case_when(CHROM_hg19 == 0 ~ CHROM, CHROM_hg19 != 0 ~ CHROM_hg19)) %>%
mutate(CHROM_hg19 = gsub("chr", "", CHROM_hg19, perl = TRUE)) %>%
mutate(CHROM_hg19 = case_when( !is.na(CHROM_hg19) ~ paste0("chr",CHROM_hg19), is.na(CHROM_hg19) ~ CHROM_hg19)) %>%
mutate(gADe211_AF = case_when( !is.na(gADe211_AF) ~ round(gADe211_AF,8))) %>%
mutate(gADg30_AF = case_when( !is.na(gADg30_AF) ~ round(gADg30_AF,8))) %>%
mutate(bravo_AF = case_when( !is.na(bravo_AF) ~ round(bravo_AF,8))) %>%
mutate(MetaLR = case_when( !is.na(MetaLR) ~ round(MetaLR,3))) %>%
mutate(MetaSVM = case_when( !is.na(MetaSVM) ~ round(MetaSVM,3))) %>%
mutate(M_CAP = case_when( !is.na(M_CAP) ~ round(M_CAP,3))) %>%
mutate(REVEL = case_when( !is.na(REVEL) ~ round(REVEL,3))) %>%
mutate(MVP = case_when( !is.na(MVP) ~ round(MVP,3))) %>%
mutate(MPC = case_when( !is.na(MPC) ~ round(MPC,3))) %>%
mutate(PrimateAI = case_when( !is.na(PrimateAI) ~ round(PrimateAI,3))) %>%
mutate(MutationTaster = case_when( !is.na(MutationTaster) ~ round(MutationTaster,3))) %>%
mutate(dbscSNV_ada = case_when( !is.na(dbscSNV_ada) ~ round(dbscSNV_ada,3))) %>%
mutate(Variant_assessed = "") %>%
mutate(Variant_relevance = "") %>%
mutate(Variant_comment = "")

##
conNDD_homozygous_gene <- conNDD_homozygous %>%
dplyr::left_join(dbNSFP4_and_gevir, by = c("GENE" = "Gene_name"))

##
conNDD_homozygous_gene_regions <- conNDD_homozygous_gene %>%
dplyr::left_join(Supplementary_Information_JAMA_Psychiatry_2020_02_06, by = c("DNA_Identifier" = "DNA_Identifier"), suffix = c("", ".y")) %>%
dplyr::select(-project_FamilyID.y,-IndividualID.y) %>%
fuzzy_left_join(HomozygosityRegions_2020_02_06, by = c("project_FamilyID" = "project_FamilyID", "CHROM" = "chrom", "POS" = "chromStart", "POS" = "chromEnd"),  match_fun = list(`==`, `==`, `>=`, `<=`)) %>%
dplyr::select(-chrom,-chromStart,-chromEnd,-project_FamilyID.y,-FamilyID_AlternativeNomenclature,-Family_analysis) %>%
dplyr::rename(project_FamilyID = project_FamilyID.x)

View(conNDD_homozygous_gene_regions)


#############
## dominant by SysID lists

conNDD_dominant <- dominant_all %>%
mutate(source = str_remove(source, "filtering/extracted/")) %>%
mutate(source = str_remove(source, ".txt")) %>%
separate(source, into = c("File", "Analysis"), sep = "\\.", remove = FALSE) %>%
mutate(File = gsub("_s", "_", File, perl = TRUE)) %>%
separate(File, into = c("project_FamilyID", "IndividualID", "DNA_Identifier"), sep = "_") %>%
dplyr::rename(CHROM_hg19 = `dbNSFP_hg19_chr`, POS_hg19 = `dbNSFP_hg19_pos_1_based_`, GENE = `ANN[0].GENE`, FEATUREID = `ANN[0].FEATUREID`, EFFECT = `ANN[0].EFFECT`, IMPACT = `ANN[0].IMPACT`, HGVS_C = `ANN[0].HGVS_C`, HGVS_P = `ANN[0].HGVS_P`, gADe211_AC = `gADe211_AC[0]`, gADe211_AN = `gADe211_AN[0]`, gADe211_AF = `gADe211_AF[0]`, gADg30_AC = `gADg30_AC[0]`, gADg30_AN = `gADg30_AN[0]`, gADg30_AF = `gADg30_AF[0]`, bravo_AC = `bravo_AC[0]`, bravo_AN = `bravo_AN[0]`, bravo_AF = `bravo_AF[0]`, CADDv15_PHRED = `CADDv15_PHRED[0]`, MetaLR = `dbNSFP_MetaLR_score[0]`, MetaSVM = `dbNSFP_MetaSVM_score[0]`, M_CAP = `dbNSFP_M_CAP_score[0]`, REVEL = `dbNSFP_REVEL_score[0]`, MVP = `dbNSFP_MVP_score[0]`, MPC = `dbNSFP_MPC_score[0]`, PrimateAI = `dbNSFP_PrimateAI_score[0]`, SIFT = `dbNSFP_SIFT_score[0]`, Polyphen2_HVAR = `dbNSFP_Polyphen2_HVAR_score[0]`, MutationTaster = `dbNSFP_MutationTaster_score[0]`, spidex_dpsi = `spidex_dpsi_zscore[0]`, dbscSNV_ada = `dbscSNV_ada_score[0]`, dbscSNV_rf = `dbscSNV_rf_score[0]`, clinvar_20200127 = `clinvar_20200127_CLNSIG[0]`, hgmd_20193 = `hgmd_20193_variant_type[0]`, index_GT = `GEN[index].GT`, index_AD = `GEN[index].AD[0]`, index_DP = `GEN[index].DP[0]`, mother_GT = `GEN[mother].GT`, father_GT = `GEN[father].GT`, other_GT = `GEN[other].GT`) %>%
mutate(frac_ADDP = round(1-(index_AD/index_DP),8)) %>%
mutate(p_ADDP = round(pbinom(index_AD,index_DP,0.5),8)) %>%
dplyr::select(project_FamilyID, IndividualID, DNA_Identifier, Analysis, CHROM, POS, ID, REF, ALT, CHROM_hg19, POS_hg19, FILTER, QUAL, AC, AF, AN, index_GT, index_AD, index_DP, frac_ADDP, p_ADDP, mother_GT, father_GT, other_GT, GENE, FEATUREID, EFFECT, IMPACT, HGVS_C, HGVS_P, gADe211_AC, gADe211_AN, gADe211_AF, gADg30_AC, gADg30_AN, gADg30_AF, bravo_AC, bravo_AN, bravo_AF, CADDv15_PHRED, MetaLR, MetaSVM, M_CAP, REVEL, MVP, MPC, PrimateAI, SIFT, Polyphen2_HVAR, MutationTaster, spidex_dpsi, dbscSNV_ada, dbscSNV_rf, clinvar_20200127, hgmd_20193) %>%
dplyr::distinct() %>%
group_by(.dots=c("project_FamilyID", "IndividualID", "DNA_Identifier", "CHROM", "POS", "REF", "ALT")) %>%
mutate(Analysis = paste(Analysis, collapse = ',')) %>%
dplyr::distinct() %>% ungroup %>%
mutate(ScoresCount = ifelse(ifelse(is.na(MetaLR),0,MetaLR) >= 0.5,1,0)+ifelse(ifelse(is.na(MetaSVM),0,MetaSVM) > 0,1,0)+ifelse(ifelse(is.na(M_CAP),0,M_CAP) >= 0.025,1,0)+ifelse(ifelse(is.na(REVEL),0,REVEL) >= 0.5,1,0)+ifelse(ifelse(is.na(MVP),0,MVP) >= 0.15,1,0)+ifelse(ifelse(is.na(MPC),0,MPC) >= 2,1,0)+ifelse(ifelse(is.na(PrimateAI),0,PrimateAI) >= 0.803,1,0)) %>%
separate(POS_hg19, into = c("POS_hg19"), sep = ",", remove = TRUE) %>%
separate(CHROM_hg19, into = c("CHROM_hg19"), sep = ",", remove = TRUE) %>%
mutate(CHROM_hg19 = case_when(CHROM_hg19 == 0 ~ CHROM, CHROM_hg19 != 0 ~ CHROM_hg19)) %>%
mutate(CHROM_hg19 = gsub("chr", "", CHROM_hg19, perl = TRUE)) %>%
mutate(CHROM_hg19 = case_when( !is.na(CHROM_hg19) ~ paste0("chr",CHROM_hg19), is.na(CHROM_hg19) ~ CHROM_hg19)) %>%
mutate(gADe211_AF = case_when( !is.na(gADe211_AF) ~ round(gADe211_AF,8))) %>%
mutate(gADg30_AF = case_when( !is.na(gADg30_AF) ~ round(gADg30_AF,8))) %>%
mutate(bravo_AF = case_when( !is.na(bravo_AF) ~ round(bravo_AF,8))) %>%
mutate(MetaLR = case_when( !is.na(MetaLR) ~ round(MetaLR,3))) %>%
mutate(MetaSVM = case_when( !is.na(MetaSVM) ~ round(MetaSVM,3))) %>%
mutate(M_CAP = case_when( !is.na(M_CAP) ~ round(M_CAP,3))) %>%
mutate(REVEL = case_when( !is.na(REVEL) ~ round(REVEL,3))) %>%
mutate(MVP = case_when( !is.na(MVP) ~ round(MVP,3))) %>%
mutate(MPC = case_when( !is.na(MPC) ~ round(MPC,3))) %>%
mutate(PrimateAI = case_when( !is.na(PrimateAI) ~ round(PrimateAI,3))) %>%
mutate(MutationTaster = case_when( !is.na(MutationTaster) ~ round(MutationTaster,3))) %>%
mutate(dbscSNV_ada = case_when( !is.na(dbscSNV_ada) ~ round(dbscSNV_ada,3))) %>%
mutate(Variant_assessed = "") %>%
mutate(Variant_relevance = "") %>%
mutate(Variant_comment = "")

##
conNDD_dominant_gene <- conNDD_dominant %>%
dplyr::left_join(dbNSFP4_and_gevir, by = c("GENE" = "Gene_name"))

##
conNDD_dominant_gene_regions <- conNDD_dominant_gene %>%
dplyr::left_join(Supplementary_Information_JAMA_Psychiatry_2020_02_06, by = c("DNA_Identifier" = "DNA_Identifier"), suffix = c("", ".y")) %>%
dplyr::select(-project_FamilyID.y,-IndividualID.y)

View(conNDD_dominant_gene_regions)


#############
## recessive by SysID lists

conNDD_recessive <- recessive_all %>%
mutate(source = str_remove(source, "filtering/extracted/")) %>%
mutate(source = str_remove(source, ".txt")) %>%
separate(source, into = c("File", "Analysis"), sep = "\\.", remove = FALSE) %>%
mutate(File = gsub("_s", "_", File, perl = TRUE)) %>%
separate(File, into = c("project_FamilyID", "IndividualID", "DNA_Identifier"), sep = "_") %>%
dplyr::rename(CHROM_hg19 = `dbNSFP_hg19_chr`, POS_hg19 = `dbNSFP_hg19_pos_1_based_`, GENE = `ANN[0].GENE`, FEATUREID = `ANN[0].FEATUREID`, EFFECT = `ANN[0].EFFECT`, IMPACT = `ANN[0].IMPACT`, HGVS_C = `ANN[0].HGVS_C`, HGVS_P = `ANN[0].HGVS_P`, gADe211_AC = `gADe211_AC[0]`, gADe211_AN = `gADe211_AN[0]`, gADe211_AF = `gADe211_AF[0]`, gADg30_AC = `gADg30_AC[0]`, gADg30_AN = `gADg30_AN[0]`, gADg30_AF = `gADg30_AF[0]`, bravo_AC = `bravo_AC[0]`, bravo_AN = `bravo_AN[0]`, bravo_AF = `bravo_AF[0]`, CADDv15_PHRED = `CADDv15_PHRED[0]`, MetaLR = `dbNSFP_MetaLR_score[0]`, MetaSVM = `dbNSFP_MetaSVM_score[0]`, M_CAP = `dbNSFP_M_CAP_score[0]`, REVEL = `dbNSFP_REVEL_score[0]`, MVP = `dbNSFP_MVP_score[0]`, MPC = `dbNSFP_MPC_score[0]`, PrimateAI = `dbNSFP_PrimateAI_score[0]`, SIFT = `dbNSFP_SIFT_score[0]`, Polyphen2_HVAR = `dbNSFP_Polyphen2_HVAR_score[0]`, MutationTaster = `dbNSFP_MutationTaster_score[0]`, spidex_dpsi = `spidex_dpsi_zscore[0]`, dbscSNV_ada = `dbscSNV_ada_score[0]`, dbscSNV_rf = `dbscSNV_rf_score[0]`, clinvar_20200127 = `clinvar_20200127_CLNSIG[0]`, hgmd_20193 = `hgmd_20193_variant_type[0]`, index_GT = `GEN[index].GT`, index_AD = `GEN[index].AD[0]`, index_DP = `GEN[index].DP[0]`, mother_GT = `GEN[mother].GT`, father_GT = `GEN[father].GT`, other_GT = `GEN[other].GT`) %>%
mutate(frac_ADDP = round(1-(index_AD/index_DP),8)) %>%
mutate(p_ADDP = round(pbinom(index_AD,index_DP,0.5),8)) %>%
dplyr::select(project_FamilyID, IndividualID, DNA_Identifier, Analysis, CHROM, POS, ID, REF, ALT, CHROM_hg19, POS_hg19, FILTER, QUAL, AC, AF, AN, index_GT, index_AD, index_DP, frac_ADDP, p_ADDP, mother_GT, father_GT, other_GT, GENE, FEATUREID, EFFECT, IMPACT, HGVS_C, HGVS_P, gADe211_AC, gADe211_AN, gADe211_AF, gADg30_AC, gADg30_AN, gADg30_AF, bravo_AC, bravo_AN, bravo_AF, CADDv15_PHRED, MetaLR, MetaSVM, M_CAP, REVEL, MVP, MPC, PrimateAI, SIFT, Polyphen2_HVAR, MutationTaster, spidex_dpsi, dbscSNV_ada, dbscSNV_rf, clinvar_20200127, hgmd_20193) %>%
dplyr::distinct() %>%
group_by(.dots=c("project_FamilyID", "IndividualID", "DNA_Identifier", "CHROM", "POS", "REF", "ALT")) %>%
mutate(Analysis = paste(Analysis, collapse = ',')) %>%
dplyr::distinct() %>% ungroup %>%
mutate(ScoresCount = ifelse(ifelse(is.na(MetaLR),0,MetaLR) >= 0.5,1,0)+ifelse(ifelse(is.na(MetaSVM),0,MetaSVM) > 0,1,0)+ifelse(ifelse(is.na(M_CAP),0,M_CAP) >= 0.025,1,0)+ifelse(ifelse(is.na(REVEL),0,REVEL) >= 0.5,1,0)+ifelse(ifelse(is.na(MVP),0,MVP) >= 0.15,1,0)+ifelse(ifelse(is.na(MPC),0,MPC) >= 2,1,0)+ifelse(ifelse(is.na(PrimateAI),0,PrimateAI) >= 0.803,1,0)) %>%
separate(POS_hg19, into = c("POS_hg19"), sep = ",", remove = TRUE) %>%
separate(CHROM_hg19, into = c("CHROM_hg19"), sep = ",", remove = TRUE) %>%
mutate(CHROM_hg19 = case_when(CHROM_hg19 == 0 ~ CHROM, CHROM_hg19 != 0 ~ CHROM_hg19)) %>%
mutate(CHROM_hg19 = gsub("chr", "", CHROM_hg19, perl = TRUE)) %>%
mutate(CHROM_hg19 = case_when( !is.na(CHROM_hg19) ~ paste0("chr",CHROM_hg19), is.na(CHROM_hg19) ~ CHROM_hg19)) %>%
mutate(gADe211_AF = case_when( !is.na(gADe211_AF) ~ round(gADe211_AF,8))) %>%
mutate(gADg30_AF = case_when( !is.na(gADg30_AF) ~ round(gADg30_AF,8))) %>%
mutate(bravo_AF = case_when( !is.na(bravo_AF) ~ round(bravo_AF,8))) %>%
mutate(MetaLR = case_when( !is.na(MetaLR) ~ round(MetaLR,3))) %>%
mutate(MetaSVM = case_when( !is.na(MetaSVM) ~ round(MetaSVM,3))) %>%
mutate(M_CAP = case_when( !is.na(M_CAP) ~ round(M_CAP,3))) %>%
mutate(REVEL = case_when( !is.na(REVEL) ~ round(REVEL,3))) %>%
mutate(MVP = case_when( !is.na(MVP) ~ round(MVP,3))) %>%
mutate(MPC = case_when( !is.na(MPC) ~ round(MPC,3))) %>%
mutate(PrimateAI = case_when( !is.na(PrimateAI) ~ round(PrimateAI,3))) %>%
mutate(MutationTaster = case_when( !is.na(MutationTaster) ~ round(MutationTaster,3))) %>%
mutate(dbscSNV_ada = case_when( !is.na(dbscSNV_ada) ~ round(dbscSNV_ada,3))) %>%
mutate(Variant_assessed = "") %>%
mutate(Variant_relevance = "") %>%
mutate(Variant_comment = "") %>%
add_count(GENE, IndividualID, name = "GeneCountPerIndividual") %>%
dplyr::filter(GeneCountPerIndividual > 1) %>%
dplyr::select(-GeneCountPerIndividual)

##
conNDD_recessive_gene <- conNDD_recessive %>%
dplyr::left_join(dbNSFP4_and_gevir, by = c("GENE" = "Gene_name"))

##
conNDD_recessive_gene_regions <- conNDD_recessive_gene %>%
dplyr::left_join(Supplementary_Information_JAMA_Psychiatry_2020_02_06, by = c("DNA_Identifier" = "DNA_Identifier"), suffix = c("", ".y")) %>%
dplyr::select(-project_FamilyID.y,-IndividualID.y)

View(conNDD_recessive_gene_regions)


###########################################################
## export for excel
write_excel_csv2(conNDD_homozygous_gene_regions, paste0("results/conNDD_homozygous.",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "", append = FALSE, col_names = T, quote_escape = "double")
write_excel_csv2(conNDD_dominant_gene_regions, paste0("results/conNDD_dominant.",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "", append = FALSE, col_names = T, quote_escape = "double")
write_excel_csv2(conNDD_recessive_gene_regions, paste0("results/conNDD_recessive.",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "", append = FALSE, col_names = T, quote_escape = "double")


