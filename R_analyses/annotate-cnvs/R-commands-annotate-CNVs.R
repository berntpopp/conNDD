####################################################################
## load libraries
library(readxl)			## for Excel related data import
library(tidyverse)		## for data import and transformation
library(jsonlite)		## needed for HGNC requests
library(fuzzyjoin)		## for gene region joins
####################################################################


####################################################################
## define work directory and excel file locations

##set the work directoty from config file
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir))
####################################################################


####################################################################
## load CNVs from Excel file

#2021-02-16
variants_file <- "data/FileS3_conNDD-variants_2021-01-16.xlsx"

CNVSyndromes_file <- "data/Decipher_CNVSyndomesList_2021-02-02.xlsx"

# the file FileS3_conNDD-variants.xlsx needs to be downloaded in current version from Zenodo
FileS2_location <- "data/FileS2_conNDD-cohort.xlsx"

variants_data <- read_excel(FileS2_location, sheet = "CNVkit_one-vs-all", col_types = c("text", "text", "text", "text", "text", "text", "text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",	"numeric", "text", "text", "text"))

# read CNV syndromes file
CNVSyndromes_data <- read_excel(CNVSyndromes_file, sheet = "CNVSyndomesList")

## import sysid data
sysid_db_humangene_collected <- read_csv("data/sysid_3_human_gene_2021-04-18.csv")
sysid_db_disease_collected <- read_csv("data/sysid_3_disease_2021-04-18.csv")
####################################################################


############################################
## define functions
hgnc_id_from_prevsymbol <- function(symbol_input)  {
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/prev_symbol/", symbol_input))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) symbol else "") %>%
	mutate(score = if (exists('score', where = hgnc_id_from_symbol)) score else 0) %>%
	arrange(desc(score)) %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))

	return(as.integer(hgnc_id_from_symbol$hgnc_id[1]))
}

hgnc_id_from_aliassymbol <- function(symbol_input)  {
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/alias_symbol/", symbol_input))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) symbol else "") %>%
	mutate(score = if (exists('score', where = hgnc_id_from_symbol)) score else 0) %>%
	arrange(desc(score)) %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))

	return(as.integer(hgnc_id_from_symbol$hgnc_id[1]))
}

hgnc_id_from_symbol <- function(symbol_tibble) {
	symbol_list_tibble <- as_tibble(symbol_tibble) %>% select(symbol = value) %>% mutate(symbol = toupper(symbol))
	
	symbol_request <- fromJSON(paste0("http://rest.genenames.org/search/symbol/", str_c(symbol_list_tibble$symbol, collapse = "+OR+")))

	hgnc_id_from_symbol <- as_tibble(symbol_request$response$docs)
	
	hgnc_id_from_symbol <- hgnc_id_from_symbol %>%
	mutate(hgnc_id = if (exists('hgnc_id', where = hgnc_id_from_symbol)) hgnc_id else NA) %>%
	mutate(symbol = if (exists('symbol', where = hgnc_id_from_symbol)) toupper(symbol) else "") %>%
	mutate(hgnc_id = as.integer(str_split_fixed(hgnc_id, ":", 2)[, 2]))
		
	return_tibble <- symbol_list_tibble %>% 
	left_join(hgnc_id_from_symbol, by = "symbol") %>%
	select(hgnc_id)

	return(return_tibble)
}	

hgnc_id_from_symbol_grouped <- function(input_tibble, request_max = 150) {
	input_tibble <- as_tibble(input_tibble)
	
	row_number <- nrow(input_tibble)
	groups_number <- ceiling(row_number/request_max)
	
	input_tibble_request <- input_tibble %>%
	mutate(group = sample(1:groups_number, row_number, replace=T)) %>%
	group_by(group) %>%
	mutate(response = hgnc_id_from_symbol(value)$hgnc_id) %>%
	ungroup()
	
	input_tibble_request_repair <- input_tibble_request %>%
	filter(is.na(response)) %>%
	select(value) %>%
	unique() %>%
	rowwise() %>%
	mutate(response = hgnc_id_from_prevsymbol(value)) %>%
	mutate(response = case_when(!is.na(response) ~ response, is.na(response) ~ hgnc_id_from_aliassymbol(value)))
	
	input_tibble_request <- input_tibble_request %>%
	left_join(input_tibble_request_repair, by = "value") %>%
	mutate(response = case_when(!is.na(response.x) ~ response.x, is.na(response.x) ~ response.y))
	
	return(input_tibble_request$response)
}


gene_info_from_hgnc_id <- function(id_tibble) 
{
	id_list_tibble <- as_tibble(id_tibble) %>% select(id = value)
	
	hgnc_id_request <- fromJSON(paste0("http://rest.genenames.org/fetch/hgnc_id/", str_c(id_list_tibble$id, collapse = "+OR+")))
	
	gene_info_from_hgnc_id <- as_tibble(hgnc_id_request$response$docs) 
	
	gene_info_from_hgnc_id_output <- gene_info_from_hgnc_id %>%
	mutate(prev_symbol = if (exists('prev_symbol', where = gene_info_from_hgnc_id)) prev_symbol else NA) %>%
	mutate(alias_symbol = if (exists('alias_symbol', where = gene_info_from_hgnc_id)) alias_symbol else NA) %>%
	mutate(gene_synonyms = str_c(unlist(c(prev_symbol, alias_symbol)), collapse = ", ")) %>%
	mutate(gene_synonyms = replace_na(gene_synonyms, "-")) %>%
	select(entrez_id, location, gene_symbol = symbol, gene_description = name, gene_synonyms, ensembl_id = ensembl_gene_id) 

	return(gene_info_from_hgnc_id_output)
}
############################################


############################################
## download and reformat of OMIM files and reformat
omim_file_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

# read the text file with the links to the OMIM files
# CAVE: this file has toi be adapted with your own download links
omim_links <- as.tibble(read_lines("data/omim_links.txt")) %>%
	mutate(file_name = str_remove(value, "https.+\\/")) %>%
	mutate(file_name = str_remove(file_name, "\\.txt")) %>%
	mutate(file_name = paste0("data/", file_name, ".", omim_file_date, ".txt"))

downlaod OMIM files
for (row in 1:nrow(omim_links)) {
	download.file(omim_links$value[row], omim_links$file_name[row], mode = "wb")
}

## reformat
mim2gene <- read_delim(omim_links$file_name[1], "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE) %>%
	select(MIM_Number = X1,	MIM_Entry_Type = X2, Entrez_Gene_ID = X3, Approved_Gene_Symbol_HGNC = X4, Ensembl_Gene_ID = X5) %>%
	mutate(omim_id = paste0("OMIM:",MIM_Number)) %>%
	select(-MIM_Number, -Entrez_Gene_ID, -Ensembl_Gene_ID)

mim2gene_hgnc <- mim2gene %>%
	filter(!is.na(Approved_Gene_Symbol_HGNC)) %>%
	unique() %>%
	mutate(hgnc_id = paste0("HGNC:",hgnc_id_from_symbol_grouped(Approved_Gene_Symbol_HGNC)))
	
mimTitles <- read_delim(omim_links$file_name[2], "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE) %>%
	select(Prefix = X1,	MIM_Number = X2, Preferred_Title_symbol = X3, Alternative_Titles_symbols = X4, Included_Titles_symbols = X5) %>%
	mutate(omim_id = paste0("OMIM:",MIM_Number)) %>%
	select(-MIM_Number)

morbidmap <- read_delim(omim_links$file_name[4], "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE) %>%
	select(Phenotype = X1,	Gene_Symbols = X2, Gene_MIM_Number = X3, Cyto_Location = X4) %>%
	select(-Gene_Symbols, -Cyto_Location) %>%
	mutate(Phenotype_AssociationKey = str_remove(Phenotype, ".+ \\(")) %>%
	mutate(Phenotype_AssociationKey = str_remove(Phenotype_AssociationKey, "\\)")) %>%
	mutate(Phenotype_MIM_Number = str_extract(Phenotype, " [0-9][0-9][0-9][0-9][0-9][0-9]")) %>%
	mutate(Phenotype_MIM_Number = as.integer(str_remove(Phenotype_MIM_Number, " "))) %>%
	mutate(Phenotype_HasSpecialCharacter = str_detect(Phenotype, "[\\?\\[\\{]")) %>%
	mutate(Phenotype = str_remove(Phenotype, ", [0-9][0-9][0-9][0-9][0-9][0-9].+")) %>%
	mutate(omim_p_id = case_when(!is.na(Phenotype_MIM_Number) ~ paste0("OMIM:",Phenotype_MIM_Number), is.na(Phenotype_MIM_Number) ~ "")) %>%
	mutate(omim_g_id = case_when(!is.na(Gene_MIM_Number) ~ paste0("OMIM:",Gene_MIM_Number), is.na(Gene_MIM_Number) ~ "")) %>%
	mutate(omim_p_id = na_if(omim_p_id, "")) %>%
	mutate(omim_g_id = na_if(omim_g_id, "")) %>%
	select(-Phenotype_MIM_Number, -Gene_MIM_Number) %>%
	filter(!is.na(omim_p_id))

morbidmap_mim2gene_hgnc <- morbidmap %>%
	left_join(mim2gene_hgnc, by = c("omim_g_id" = "omim_id")) %>%
	select(-Phenotype_AssociationKey, -MIM_Entry_Type, -Phenotype_HasSpecialCharacter, -omim_g_id)

############################################


## on 2021-04-18
############################################
## SysID: load the human_gene and disease tables exported from the SysID database

## reformat the table

repair_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

sysid_db_humangene_collected_repaired <- sysid_db_humangene_collected %>%
	mutate(sysid_id2 = as.integer(str_remove(sysid_id, "SysID_"))) %>%
	arrange(sysid_id2, human_gene_id) %>% 
	mutate(sysid_id3 = as.integer(row_number()+11)) %>%
	mutate(sysid_id4 = case_when(!is.na(sysid_id2) ~ sysid_id2, is.na(sysid_id2) ~ sysid_id3)) %>%
	mutate(sysid_id = paste0("SysID_" ,str_pad(sysid_id4, 4, pad = "0"))) %>%
	select(-sysid_id2, -sysid_id3, -sysid_id4) %>%
	arrange(human_gene_id) %>%
	mutate(human_gene_remark = na_if(human_gene_remark, "")) %>%
	mutate(human_gene_remark = case_when(!is.na(human_gene_remark) ~ paste0(human_gene_remark, "; gene information repaired on ", repair_date), is.na(human_gene_remark) ~ paste0("gene information repaired on ", repair_date)))

sysid_db_humangene_collected_repaired_hgnc <- sysid_db_humangene_collected_repaired %>%
	mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene_symbol))

sysid_db_humangene_collected_repaired_hgnc_new <- sysid_db_humangene_collected_repaired_hgnc %>% 
	rowwise() %>% 
	mutate(new = gene_info_from_hgnc_id(hgnc_id))

sysid_db_humangene_collected_repaired_hgnc_new_select <- sysid_db_humangene_collected_repaired_hgnc_new %>%
	mutate(gene_symbol = new$gene_symbol) %>%
	left_join(sysid_db_disease_collected, by = c("human_gene_id" = "human_gene_id")) %>%
	select(gene_symbol = gene_symbol.x, gene_group)
############################################


############################################
## join the tables
variants_omim_p_ids <- variants_data %>%
	filter(gene != "-") %>%
	select(VarID, gene) %>%
	separate_rows(gene, sep = ",") %>%
	left_join(morbidmap_mim2gene_hgnc, by = c("gene" = "Approved_Gene_Symbol_HGNC")) %>%
	select(VarID, omim_p_id) %>%
	filter(!is.na(omim_p_id)) %>%
	group_by(VarID) %>%
	summarize(omim_p_ids = str_c(omim_p_id, collapse = ", "))

variants_sysid_db_genes <- variants_data %>%
	filter(gene != "-") %>%
	select(VarID, gene) %>%
	separate_rows(gene, sep = ",") %>%
	left_join(sysid_db_humangene_collected_repaired_hgnc_new_select, by = c("gene" = "gene_symbol")) %>%
	select(VarID, gene_group) %>%
	filter(!is.na(gene_group)) %>%
	group_by(VarID) %>%
	summarize(gene_groups = str_c(gene_group, collapse = ", "))

variants_CNVSyndromesAndSysID <- variants_data %>%
	fuzzy_left_join(CNVSyndromes_data, by = c("chromosome" = "chromosome", "start" = "start", "end" = "end"), match_fun = list(`==`, `>=`, `<=`)) %>%
	left_join(variants_omim_p_ids, by = c("VarID" = "VarID")) %>%
	left_join(variants_sysid_db_genes, by = c("VarID" = "VarID")) %>%
	select(VarID, project_FamilyID, IndividualID, DNA_Identifier, Type, Run, chromosome = chromosome.x, start = start.x, end = end.x, gene, log2, cn, depth, probes, weight, CNVPerSample, AffectedGenesCountInCohort, Variant_assessed, Variant_relevance, Variant_comment, SyndromeName, GenotypeAndClass, omim_p_ids, gene_groups)
############################################


############################################
## export for excel

#2021-02-16
write_excel_csv2(variants_CNVSyndromesAndSysID, paste0("results/variants_CNVSyndromesAndSysID.",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "", append = FALSE, col_names = T, quote_escape = "double")
