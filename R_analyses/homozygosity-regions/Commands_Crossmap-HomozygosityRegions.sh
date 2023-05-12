#!/bin/sh

#########################################
# initialize anaconda
export PATH="$(config_get conda_path)"
#########################################


#########################################
# go to workdirektory
cd "$(config_get work_directory)"
#########################################


## lifting HomozygosityRegions to hg38
CrossMap.py bed ChainFiles/hg19ToHg38.over.chain.gz data/HomozygosityRegions_2020-02-03.bed results/HomozygosityRegions_2020-02-03.hg38.bed
