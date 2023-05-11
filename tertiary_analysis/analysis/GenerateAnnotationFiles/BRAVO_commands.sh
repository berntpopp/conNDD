#!/bin/sh

#############
## 2020-01-27

# go to workdirektory
cd "$(config_get work_directory)"

# go to annotation folder
$ cd analysis/annotation/hg38/bravo/

# download files
curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _ga=GA1.2.736518053.1577369830; remember_token="bernt.popp@gmail.com|d94d053fe12c4913950374f7647a2ec1530332dc0c194c0d3cc8261dc71946be1508e71b8da6828edda289092332ca6060bc39fd23c8ad92c41822c7d52d32c8"; _gid=GA1.2.1964649665.1580115153; _gat_gtag_UA_73910830_2=1' --compressed > bravo-dbsnp-all.vcf.gz

tabix -f -p vcf bravo-dbsnp-all.vcf.gz

#########################################
# subset to target for faster annotation:
sort -k 1,1 -k2,2n -k3,2n  "$(config_get work_directory)"/analysis/targets/SureSelectBED/hg19/S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.bed > S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed

bcftools view --threads 24 -O z -l 1 -T S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.bed bravo-dbsnp-all.vcf.gz -o bravo-dbsnp-all.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz

tabix -f -p vcf bravo-dbsnp-all.S03723314-and-S04380110-and-S07604514_Covered.hg38.lifted.plusChrM.addedlines.plusminus100.sorted.vcf.gz