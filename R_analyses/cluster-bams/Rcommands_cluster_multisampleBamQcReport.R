## based on https://cmdlinetips.com/2019/05/how-to-do-pca-in-tidyverse-framework/
## based on https://broom.tidyverse.org/reference/tidy.prcomp.html

##################################
## load library
library(tidyverse)
library(readxl)
library(broom)
library(scales)
library(config)			## used for config loading


##################################
## set working directory
## define work directory and excel file locations

##set the work directory from config file
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir))


##################################
## set random seed for session
set.seed(23)

##################################
# read data files

# the file FileS2_conNDD-cohort.xlsx needs to be downloaded in current version from Zenodo
BAMs <- read_excel("data/FileS2_conNDD-cohort.xlsx", sheet = "BAMs") %>%
select(Run, DNAid, PL_attribute)

## load qualimap multisampleBamQcReportdata
multisampleBamQcReport <- read_excel("data/multisampleBamQcReport.xlsx")

multisampleBamQcReport_Samples_groups <- read_excel("data/multisampleBamQcReport.xlsx", sheet = "Samples_groups")

####################################################################
## do k-means clustering on raw multidemensional data
## https://uc-r.github.io/kmeans_clustering

library(cluster)
library(factoextra)

multisampleBamQcReport_sample <- multisampleBamQcReport %>% 
  unite("Sample_Group", c(Sample,Group)) %>%
  column_to_rownames("Sample_Group")

multisampleBamQcReport_sample_scaled <- scale(multisampleBamQcReport_sample)

k4 <- kmeans(multisampleBamQcReport_sample, centers = 4, nstart = 25)

fviz_cluster(k4, data = multisampleBamQcReport_sample, ggtheme = theme_classic(), stand = FALSE)

multisampleBamQcReport_cluster <- add_column(multisampleBamQcReport, k4$cluster)

###################################################
## export table
write_csv2(multisampleBamQcReport_cluster, paste0("multisampleBamQcReport_cluster.3","_",format(Sys.time(), "%Y-%m-%d"),".csv"), na = "")

###################################################
## join with BAMs table

multisampleBamQcReport_cluster_file <- left_join(multisampleBamQcReport_cluster, multisampleBamQcReport_Samples_groups, by = c("Sample", "Group")) %>%
  mutate(File = paste0(File,".MarkDups.bam")) %>%
  mutate(File = str_replace(File, "results/QC/qualimap/mergedbed/", "results/BAMs/hg38/MarkDups/"))

# generate worklist for CNVkit one-vs-all
OneVsAll_worklist <- multisampleBamQcReport_cluster_file %>%
  select(File, `k4$cluster`) %>%
  group_by(`k4$cluster`) %>%
  summarize(Control_Files = stringr::str_c(File, collapse = " ")) %>%
  ungroup() %>%
  right_join(multisampleBamQcReport_cluster_file, by = "k4$cluster") %>%
  select(Sample, Group, File, Control_Files, Cluster = `k4$cluster`) %>%
  mutate(Control_Files = str_replace(Control_Files, paste0(File,"\\s?"), "")) %>%
  mutate(Cluster = paste0("cluster",Cluster)) %>%
  left_join(BAMs, by = c("Sample" = "DNAid", "Group" = "Run"))

View(OneVsAll_worklist)

OneVsAll_worklist_SOLID <- OneVsAll_worklist %>%
filter(PL_attribute == "SOLID")

OneVsAll_worklist_ILLUMINA <- OneVsAll_worklist %>%
filter(PL_attribute == "ILLUMINA")

write_csv2(OneVsAll_worklist_SOLID, paste0("CNVkit.one-vs-all.worklist_SOLID", ".", format(Sys.time(), "%Y-%m-%d"), ".csv"), na = "")
write_csv2(OneVsAll_worklist_ILLUMINA, paste0("CNVkit.one-vs-all.worklist_ILLUMINA", ".", format(Sys.time(), "%Y-%m-%d"),".csv"), na = "")