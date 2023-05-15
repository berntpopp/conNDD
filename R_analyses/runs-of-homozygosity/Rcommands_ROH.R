###########################################################
## load required libraries
library(readr)
library(readxl)
library(tidyverse)
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
## load required files
conNDDcohort <- read_delim("results/ROH/conNDDcohort.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# the file FileS2_conNDD-cohort.xlsx needs to be downloaded in current version from Zenodo
FileS2_conNDDcohort_BAMs <- read_excel("data/FileS2_conNDD-cohort.xlsx", sheet = "BAMs")

FileS2_conNDDcohort_samples <- read_excel("data/FileS2_conNDD-cohort.xlsx", sheet = "samples")

FileS2_conNDDcohort_variants <- read_excel("data/FileS2_conNDD-cohort.xlsx", sheet = "variants")
###########################################################


###########################################################
## reformat the conNDDcohort table and add annotations from other tables
conNDDcohort_samples <- conNDDcohort %>%
	select(sample = `[2]Sample`) %>%
	unique()

conNDDcohort_reformat <- conNDDcohort %>%
	select(sample = `[2]Sample`, chromosome = `[3]Chromosome`, start = `[4]Start`, end = `[5]End`, length = `[6]Length (bp)`)

conNDDcohort_filter <- conNDDcohort_reformat  %>%
	select(sample, chromosome, start, end, length) %>%
filter(sample %in% as.list(slice(conNDDcohort_samples, 1:10))$sample)

conNDDcohort_regionsum <- conNDDcohort_reformat %>% 
	select(sample, length) %>%
	group_by(sample) %>%
	summarise(sum(length))

conNDDcohort_regionsum_samples <- conNDDcohort_regionsum %>%
	mutate(sample = str_sub(sample, 2, str_length(sample))) %>%
	left_join(FileS2_conNDDcohort_samples, by = c("sample" = "DNAid"))

conNDDcohort_regionsum_baminfo <- conNDDcohort_regionsum %>%
	mutate(sample = str_sub(sample, 2, str_length(sample))) %>%
	left_join(FileS2_conNDDcohort_samples, by = c("sample" = "DNAid")) %>%
	left_join(FileS2_conNDDcohort_variants, by = c("IndividualID" = "Index_IndividualID")) %>%
	mutate(FileSize = as.numeric(str_sub(FileSize, 1, str_length(FileSize)-1))) %>%
	filter(primary == "yes") %>%
	filter(NewGene_Classification == "established-gene") %>%
	filter(NewVariant_Classification %in% c("P","LP","VUS"))
###########################################################


###########################################################
## Plotting
set.seed(23)

ggplot(data = conNDDcohort_filter) + 
	geom_rect(aes(xmin = start, xmax = end, ymax = 1, ymin = 0), colour="grey60", fill = "grey60") +
	facet_grid(cols = vars(chromosome), rows = vars(sample)) +
	theme_classic()

ggplot(data = conNDDcohort_regionsum_baminfo, aes(x=PL_attribute, y=`sum(length)`)) + 
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), lwd=0.5, colour="grey20", alpha=0.5) + 
	geom_point(aes(shape=Zygosity, color=ExplainsPhenotype), size=5, position=position_jitter(width=0.2, height=0.0), alpha=0.9) +
	scale_fill_manual(values = c("#f1a340", "#998ec3")) +
	scale_color_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a")) +
	theme_classic()

ggplot(data = conNDDcohort_regionsum_baminfo, aes(x=FileSize, y=`sum(length)`, fill=PL_attribute)) + 
	geom_point(shape=21, size=3, colour="grey50", position=position_jitter(width=0.2, height=0.0), alpha=0.5) +
	stat_smooth(method = "lm", col = "red") +
	theme_classic()
###########################################################


###########################################################
## export for excel
write_excel_csv2(conNDDcohort_regionsum_samples, paste0("conNDDcohort_regionsum_samples.", format(Sys.time(), "%Y-%m-%d"),".csv"), na = "", append = FALSE, col_names = T, quote_escape = "double")