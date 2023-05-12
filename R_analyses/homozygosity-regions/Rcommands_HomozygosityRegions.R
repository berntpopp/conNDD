###########################################################
## load required libraries
library(tidyverse)		## used for data laoding and manipulation
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

# load the file
HomozygosityRegions_2020_02_03_hg38 <- read_delim("results/HomozygosityRegions_2020-02-03.hg38.bed", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

# transform
HomozygosityRegions_2020_02_03_hg38_grouped <- HomozygosityRegions_2020_02_03_hg38 %>%
	dplyr::filter(!grepl("_",X1)) %>%
	group_by(.dots=c("X4", "X1")) %>%
	mutate(X2 = min(X2), X3 = min(X3)) %>%
	dplyr::distinct() %>%
	ungroup()

# write new file
write_excel_csv2(HomozygosityRegions_2020_02_03_hg38_grouped, paste0("results/HomozygosityRegions_2020_02_03_hg38_grouped.",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "NA", append = FALSE, col_names = T, delim = ";", quote_escape = "double")
