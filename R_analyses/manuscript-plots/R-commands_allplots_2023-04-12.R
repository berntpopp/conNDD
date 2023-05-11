###########################################################
## load required libraries
library(tidyverse)		## used for data laoding and manipulation
library(readxl)			## used to load excel data
library(showtext)
library(ggalluvial)
library(RColorBrewer)
library(ggplot2)
library(shadowtext)
library(cowplot)
library(scales)
library(config)			## used for config loading
###########################################################



###########################################################
## define work directory and excel file locations

##set the work directoty from config file
## read config
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"))
## set working directory
setwd(paste0(config_vars$projectsdir, project_name, script_path))

# the file FileS2_conNDD-cohort.xlsx needs to be downloaded in current version from Zenodo 
FileS2_location <- "data/FileS2_conNDD-cohort.xlsx"

###########################################################



###########################################################
## define functions

age_in_years_double_from_ym <- function(ym_text) {
	ym_text_tibble <- as_tibble(ym_text) %>% mutate(age = tolower(value))
	## one should implement basic error handling and logik checks here
	## handling of input that is NA or does not match the pattern or obviously is not in months format
	
	return_tibble <- ym_text_tibble %>%
		mutate(value = str_replace(value, "m", "")) %>%
		mutate(value = str_split(value, pattern = "y", n = 2)) %>%
		mutate(value = purrr::map(value, setNames, c("y","m"))) %>%
		unnest_wider(value) %>%
		mutate(age_y = as.double(y), age_m = as.double(m)) %>%
		mutate(age_years = (age_y*12 + age_m)/12) %>%
		select(age_years)
	
	return(return_tibble)
}

###########################################################


###########################################################
## load the data from the excel files

FileS2_individuals_data <- read_excel(FileS2_location, sheet = "individuals", na = "NA")

FileS2_families_data <- read_excel(FileS2_location, sheet = "families", na = "NA")

FileS2_variants_data <- read_excel(FileS2_location, sheet = "variants", na = "NA")

FileS2_BAMs_data <- read_excel(FileS2_location, sheet = "BAMs", na = "NA")
 colnames(FileS2_BAMs_data)[21] <- "Coverage_mean"

FileS2_samples_data <- read_excel(FileS2_location, sheet = "samples", na = "NA")
###########################################################



###########################################################
##  A) plot of cohort regarding family structure
## get total number of simplex and multiplex cases, resp.
n_simplex <- FileS2_families_data %>%
  select(FamilyStructure, SexAffected) %>%
  group_by(FamilyStructure) %>%
  tally() %>%
  filter(FamilyStructure == "simplex") %>%
  pull(n)

n_multiplex <- FileS2_families_data %>%
  select(FamilyStructure, SexAffected) %>%
  group_by(FamilyStructure) %>%
  tally() %>%
  filter(FamilyStructure == "multiplex") %>%
  pull(n)
  
## get the information about the family structure from the family tab and calculate percentage
FileS2_families_data_family_structure <- FileS2_families_data %>%
  select(SexAffected, FamilyStructure) %>%
  group_by(FamilyStructure, SexAffected) %>%
  tally() %>%
  mutate(total_n = case_when(
    FamilyStructure == "simplex" ~ n_simplex,
    FamilyStructure == "multiplex" ~ n_multiplex)
  ) %>%
  mutate(percentage = (n/total_n*100)) %>%
  mutate(n_sex = case_when(
    SexAffected == "f" ~ paste0(intToUtf8(9792), ": ", n),
    SexAffected == "m" ~ paste0(intToUtf8(9794), ": ", n),
    SexAffected == "fm" ~ paste0(intToUtf8(9792), "/ ", intToUtf8(9794), ": ", n))
  )
view(FileS2_families_data_family_structure)


## generate stacket bar plot with family structure information
plot_family_structure <- ggplot(FileS2_families_data_family_structure, aes(x = FamilyStructure, y = percentage, fill = SexAffected)) +
geom_col( width=0.7) +
  geom_bar(position="stack", stat = "identity", width = 0.3, show.legend = FALSE) +
  scale_fill_manual(values=c("#66c2a5","#8da0cb","#fc8d62")) + 
  geom_text(aes(label= n_sex ), position = position_stack(vjust= 0.5), colour = "black", size = 5) +
  geom_text(data= FileS2_families_data_family_structure, aes(x = FamilyStructure, y = 110, label= paste0("n:", total_n ), fontface = "bold")) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), labels = function(x) paste0(x, "%")) +
  ggtitle(label= "Family structure") +
  theme_light() +
  theme(plot.title = element_text(size = 16), legend.position="none", text = element_text(size=16, family = "sans"), panel.grid.major = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.minor = element_blank())
plot_family_structure
###########################################################



###########################################################
## B) Pyramid Plot of age
## inspired by:
## https://community.rstudio.com/t/ggplot2-alter-scale-on-pyramid-plot/14934/2
## https://community.rstudio.com/t/ggplot2-alter-scale-on-pyramid-plot/14934/2

## reformat the "AgeLastInvestigation" column to be a double number in years
FileS2_individuals_data_reformat <- FileS2_individuals_data %>%
	filter(!is.na(AgeLastInvestigation)) %>%
	filter(Sex != "undetermined") %>%
	mutate(Age_years = age_in_years_double_from_ym(AgeLastInvestigation))

## 
FileS2_individuals_data_reformat_Sex_and_AgeGroups <- FileS2_individuals_data_reformat %>%
select(IndividualID, Age_years, Sex) %>%
mutate(AgeGroup  = case_when(Age_years >= 0 & Age_years < 2 ~ "0-2",
Age_years >= 2 & Age_years < 4 ~ "2-4",
Age_years >= 4 & Age_years < 6 ~ "4-6",
Age_years >= 6 & Age_years < 8 ~ "6-8",
Age_years >= 8 & Age_years < 10 ~ "8-10",
Age_years >= 10 & Age_years < 12 ~ "10-12",
Age_years >= 12 & Age_years < 14 ~ "12-14",
Age_years >= 14 & Age_years < 16 ~ "14-16",
Age_years >= 16 & Age_years < 18 ~ "16-18",
Age_years >= 18 & Age_years < 50 ~ "18-50")) %>%
group_by(AgeGroup, Sex) %>%
mutate(Sex = factor(Sex, levels = c("female", "male"))) %>%
tally() %>%
mutate(n = case_when(Sex == "male" ~ as.double(n)*(-1), Sex == "female" ~ as.double(n)))

## integrate correct sorting for AgeGroups (currently lexigraphic), and then remove the leading zeros from group name
FileS2_individuals_data_reformat_Sex_and_AgeGroups$AgeGroup <- factor(FileS2_individuals_data_reformat_Sex_and_AgeGroups$AgeGroup, levels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10-12", "12-14", "14-16", "16-18", "18-50"), ordered=TRUE)

##
FileS2_individuals_data_reformat_Sex_Count <- FileS2_individuals_data_reformat %>%
select(Sex) %>%
group_by(Sex) %>%
tally()

female_count <- paste0(intToUtf8(9792), ": ", FileS2_individuals_data_reformat_Sex_Count$n[1])
male_count <- paste0(intToUtf8(9794), ": ", FileS2_individuals_data_reformat_Sex_Count$n[2])

#############
##generate plot
##to do:
## what does "data = ." mean and can we replace it?

ctr_width <- 6
brks <- seq(0, 20, by = 5)

plot_age_pyramid <- FileS2_individuals_data_reformat_Sex_and_AgeGroups %>%
	mutate(start = if_else(Sex == "male", 0, ctr_width),
			 end = n + if_else(Sex == "male", 0, ctr_width),
			 mid = (start + end)/2,
			 wid = abs(end - start)) %>%
	ggplot() +
	geom_tile(aes(mid, AgeGroup, fill = Sex, width = wid), height = 0.8, color="black") +
	geom_text(data = . %>% distinct(AgeGroup), aes(ctr_width/2, AgeGroup, label = AgeGroup), size = 4.0) +
	scale_x_continuous(breaks = c(-brks, ctr_width + brks), labels = c(brks, brks)) +
	scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb")) + 
	geom_hline(yintercept=9.5, linetype="dashed", color = "red") +
	annotate("text", x = 7, y = 11.5, label = female_count, size = 4.0, hjust = 0, colour = "black", fontface = "bold", family = "sans") +
	annotate("text", x = -8, y = 11.5, label = male_count, size = 4.0, hjust = 0, colour = "black", fontface = "bold", family = "sans") +
	expand_limits(y = c(0, 12.0)) +
    ggtitle(label= "Age groups by sex") +
	theme_light() +
	theme(axis.text.y = element_blank()) + 
	theme(plot.title = element_text(size = 16), legend.position="none", text = element_text(size=16, family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
plot_age_pyramid

###########################################################


###########################################################
## 3) plot of cohort regarding platform
## get the information about the platform and DNA identified from the BAMs tab
FileS2_BAMs_data_platform <- FileS2_BAMs_data %>%
select(DNAid, PL_attribute, ReAnalyze, primary, Coverage_mean) 

## get the information about the platform and DNA identified from the samples tab
FileS2_samples_data_platform <- FileS2_samples_data %>%
select(DNAid, Type)

## join the information about the platform and DNA identifier from BAMs and samples tab
Files2_platform <- left_join(x = FileS2_samples_data_platform, y = FileS2_BAMs_data_platform, by = "DNAid") %>%
filter(Type == "Index", ReAnalyze == "yes", primary == "yes")
view(Files2_platform)

## count total Illumina and Solid runs
Files2_platform_count <- Files2_platform %>%
  group_by(PL_attribute)%>%
  tally()

## generate violin bar plot with platform and coverage information

plot_platform_coverage <- ggplot(Files2_platform, aes(x = PL_attribute, y =  Coverage_mean, fill = PL_attribute)) +
  geom_violin() +
  scale_fill_manual(values=c("#e9a3c9","#a1d76a")) +
  geom_jitter(width=0.15, alpha=0.5) +
  geom_text(data= Files2_platform_count, aes(x = PL_attribute, y = 320, label=paste("n:", n, sep=""), fontface = "bold")) +
  ggtitle(label= "Mean coverage by platform") +
  theme_light() +
  theme(plot.title = element_text(size = 16), legend.position="none", text = element_text(size=16, family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
plot_platform_coverage
###########################################################


###########################################################
## 4) plot of cohort regarding runs of homozygosity
## get the information about the platform identified from the BAMs tab
FileS2_BAMs_data_RoH <- FileS2_BAMs_data %>%
select(DNAid, PL_attribute, ReAnalyze, primary) 

## get the information about the runs of homozygosity and DNA identifier from the samples tab
FileS2_individuals_data_RoH <- FileS2_individuals_data %>%
select(DNAid, SumROHExome_primaryBAM)

## get the information about the platform and DNA identified from the samples tab
FileS2_samples_data_platform <- FileS2_samples_data %>%
select(DNAid, Type)

## join the information about the platform and DNA identifier from BAMs and samples tab
Files2_RoH <- left_join(x = FileS2_BAMs_data_RoH, y = FileS2_individuals_data_RoH, by = "DNAid") %>%
left_join(FileS2_samples_data_platform, by = "DNAid") %>%
filter(Type == "Index", ReAnalyze == "yes", primary == "yes") %>%
drop_na(SumROHExome_primaryBAM)
view(Files2_RoH)

## count total Illumina and Solid runs
Files2_RoH_count <- Files2_RoH %>%
  group_by(PL_attribute)%>%
  tally()

## generate violin bar plot with platform and RoH information

plot_platform_RoH <- ggplot(Files2_RoH, aes(x = PL_attribute, y =  SumROHExome_primaryBAM, fill = PL_attribute)) +
  geom_violin() +
  scale_fill_manual(values=c("#e9a3c9","#a1d76a")) +
  geom_jitter(width=0.15, alpha=0.5) +
  geom_text(data= Files2_RoH_count, aes(x = PL_attribute, y = 7.5e+08, label=paste("n:", n, sep=""), fontface = "bold")) +
  ggtitle(label= "RoH by platform") +
  theme_light() +
  theme(plot.title = element_text(size = 16), legend.position="none", text = element_text(size=16, family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_y_continuous(labels = label_bytes("MB"))

plot_platform_RoH

###########################################################



###########################################################
## 5) Alluvial Plot of variant data

FileS2_families_data_selectForAlluvial <- FileS2_families_data %>%
	select(project_FamilyID)

FileS2_variants_data_selectForAlluvial <- FileS2_variants_data %>%
	select(project_FamilyID, Index_IndividualID, OldGene_Classification, OldVariant_Classification, NewGene_Classification, NewVariant_Classification) %>%
	unite("Conclusion_Old", OldGene_Classification:OldVariant_Classification, sep = "_", remove = TRUE, na.rm = TRUE) %>%
	unite("Conclusion_New", NewGene_Classification:NewVariant_Classification, sep = "_", remove = TRUE, na.rm = TRUE) %>%
	mutate_all(na_if,"")

FileS2_families_and_variants_ForAlluvial <- FileS2_families_data_selectForAlluvial %>%
	left_join(FileS2_variants_data_selectForAlluvial, by = "project_FamilyID") %>%
	mutate(across(everything(), ~replace_na(.x, "unsolved"))) %>%
	select(-Index_IndividualID) %>%
	pivot_longer(!project_FamilyID, names_to = "Group", values_to = "Conclusion") %>%
	mutate(Conclusion = case_when(Conclusion == "established-gene_B" ~ "unsolved", Conclusion != "established-gene_B" ~ Conclusion)) %>%
	mutate(Conclusion = case_when(Conclusion == "weak candidate-gene" ~ "candidate-gene", 
		Conclusion == "published candidate-gene" ~ "candidate-gene", 
		Conclusion != "published candidate-gene" & Conclusion != "weak candidate-gene" ~ Conclusion)) %>%
	mutate(Conclusion = str_remove_all(Conclusion, "established-gene_"))

FileS2_families_and_variants_ForAlluvial$Conclusion <- factor(FileS2_families_and_variants_ForAlluvial$Conclusion, levels = c("P", "LP", "VUS", "candidate-gene", "unsolved"), ordered=TRUE)

FileS2_families_and_variants_ForAlluvial_grouped <- FileS2_families_and_variants_ForAlluvial %>%
	group_by(project_FamilyID, Group) %>% 
	mutate(Conclusion = min(Conclusion)) %>%
	ungroup() %>%
	unique() %>%
	pivot_wider(names_from = Group, values_from = Conclusion) %>%
	group_by(Conclusion_Old, Conclusion_New) %>%
	tally() %>%
	ungroup() %>%
	mutate(Families = n) %>%
	mutate(Conclusion_Old_integer = as.double(Conclusion_Old), Conclusion_New_integer = as.double(Conclusion_New)) %>%
	mutate(Conclusion_Old_integer = case_when(Conclusion_Old_integer == 2 ~ 1, Conclusion_Old_integer != 2 ~ Conclusion_Old_integer)) %>%
	mutate(Conclusion_New_integer = case_when(Conclusion_New_integer == 2 ~ 1, Conclusion_New_integer != 2 ~ Conclusion_New_integer)) %>%
	mutate(Conclusion_Old_integer = case_when(Conclusion_Old_integer == 3 ~ 4, Conclusion_Old_integer != 3 ~ Conclusion_Old_integer)) %>%
	mutate(Conclusion_New_integer = case_when(Conclusion_New_integer == 3 ~ 4, Conclusion_New_integer != 3 ~ Conclusion_New_integer)) %>%
	mutate(Clinical_change = case_when(abs(Conclusion_Old_integer - Conclusion_New_integer) > 1 ~ "yes", abs(Conclusion_Old_integer - Conclusion_New_integer) <= 1 ~ "no"))

plot_alluvial_families_conclusion_and_Clinical_change <- ggplot(FileS2_families_and_variants_ForAlluvial_grouped,
       aes(y = Families, axis1 = Conclusion_New, axis2 = Conclusion_Old)) +
  geom_alluvium(aes(fill = Clinical_change), knot.pos = 0, width = 2/12, size=0.5, color = "black") +
  geom_stratum(width = 2/12, color = "black", size=0.5, fill = c(brewer.pal(n = 5, name = 'RdBu'), brewer.pal(n = 5, name = 'RdBu'))) +
  geom_shadowtext(stat = "stratum", aes(label = after_stat(stratum)), size=4.5, color = "black", bg.colour="white", bg.r = 0.05) +
  scale_x_discrete(limits = c("2022", "2017"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  scale_fill_manual(values=c("#D1E5F0","#C95D6A")) +
  scale_y_continuous(breaks= c(seq(from = 0, to = 160, by = 20),152), sec.axis = sec_axis(~./152,labels = function(b) { paste0(round(b * 100, 0), "%")})) +
  coord_flip() +
  ggtitle(label= "Changes of assessment per family") +
  theme_light() +
  theme(plot.title = element_text(size = 16),legend.position="none", text = element_text(size=16, family = "sans"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), axis.line.y = element_line(colour = "black", arrow = grid::arrow(angle = 25, length = unit(0.5, "cm"), ends = "first", type = "closed"))) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
plot_alluvial_families_conclusion_and_Clinical_change

###########################################################
#list of families with clinical change

families_clinical_change <- FileS2_families_and_variants_ForAlluvial %>%
	group_by(project_FamilyID, Group) %>% 
	mutate(Conclusion = min(Conclusion)) %>%
	ungroup() %>%
	unique() %>%
	pivot_wider(names_from = Group, values_from = Conclusion) %>%
	group_by(Conclusion_Old, Conclusion_New) %>%
	mutate(Conclusion_Old_integer = as.double(Conclusion_Old), Conclusion_New_integer = as.double(Conclusion_New)) %>%
	mutate(Conclusion_Old_integer = case_when(Conclusion_Old_integer == 2 ~ 1, Conclusion_Old_integer != 2 ~ Conclusion_Old_integer)) %>%
	mutate(Conclusion_New_integer = case_when(Conclusion_New_integer == 2 ~ 1, Conclusion_New_integer != 2 ~ Conclusion_New_integer)) %>%
	mutate(Conclusion_Old_integer = case_when(Conclusion_Old_integer == 3 ~ 4, Conclusion_Old_integer != 3 ~ Conclusion_Old_integer)) %>%
	mutate(Conclusion_New_integer = case_when(Conclusion_New_integer == 3 ~ 4, Conclusion_New_integer != 3 ~ Conclusion_New_integer)) %>%
	mutate(Clinical_change = case_when(abs(Conclusion_Old_integer - Conclusion_New_integer) > 1 ~ "yes", abs(Conclusion_Old_integer -       Conclusion_New_integer) <= 1 ~ "no")) %>%
	filter(Clinical_change == "yes")

##########################################################
##6) plot regarding clinical change

# get percentage of clinical change
FileS2_families_clinical_change <- FileS2_families_and_variants_ForAlluvial_grouped %>%
  group_by(Clinical_change) %>%
  summarise(n = sum(Families)) %>%
  mutate(percentage = (n/sum(n)*100)) %>%
  mutate(help_col = "clinChange")
FileS2_families_clinical_change

#convert 'Clinical_change' to factor and specify level order
FileS2_families_clinical_change$Clinical_change <- factor(FileS2_families_clinical_change$Clinical_change, levels = c("yes", "no"))

## generate stacket bar plot regarding clinical change
plot_clinical_change <- ggplot(FileS2_families_clinical_change, aes(x = help_col, y = percentage, fill = forcats::fct_rev(Clinical_change))) +
  geom_bar(position="stack", stat = "identity", width = 0.5, show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#D1E5F0", "#C95D6A"),.4)) +
  geom_text(aes(label= paste0(Clinical_change, "\n","(", n, ")" )), position = position_stack(vjust= 0.5), colour = "black", size = 4) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), labels = function(x) paste0(x, "%"), position = "right", sec.axis =     sec_axis(~., labels = NULL)) +
  ggtitle(label= "Clinical change") +
  theme_light() +
  theme(plot.title = element_text(size = 16), legend.position="none", text = element_text(size=16, family = "sans"),      panel.grid.major = element_blank(), axis.title.x = element_blank(),axis.ticks.x=element_blank(), axis.text.x            =element_blank(), axis.title.y = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
plot_clinical_change

#########################
spaceholder <- ggplot() + geom_blank() + theme_void()
##########
Figure1a <- plot_grid(plot_family_structure, plot_age_pyramid, plot_platform_coverage,plot_platform_RoH , ncol = 4, align = "h", rel_widths = c(1, 1, 1, 1))
Figure1b <- plot_grid(plot_alluvial_families_conclusion_and_Clinical_change, plot_clinical_change, ncol = 2, align = "h", rel_widths = c(8, 1))
Figure1 <- plot_grid(Figure1a, spaceholder, Figure1b, nrow = 3, align = "v", labels = c("A", "B", ""), label_size = 26, rel_heights = c(4,0.4, 4))
Figure1


##########################################################
## save the plot
showtext_auto()
ggplot2::ggsave(paste0( toString(Sys.Date()),".pdf"), width = 36, height = 20, units = "cm", dpi = 300, Figure1, limitsize = FALSE)
showtext_auto(FALSE)


