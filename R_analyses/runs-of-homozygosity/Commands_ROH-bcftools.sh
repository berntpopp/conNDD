#!/bin/sh

#########################################
# initialize anaconda
```
export PATH="$(config_get conda_path)"
```
#########################################


#########################################
# go to workdirektory
cd "$(config_get work_directory)"
#########################################


#########################################
# make project folders
mkdir results/ROH/


#########################################
# call ROH with bcftools
bcftools roh results/VCFs/annotation/conNDDcohort.annoated.vcf.gz > results/ROH/conNDDcohort.txt