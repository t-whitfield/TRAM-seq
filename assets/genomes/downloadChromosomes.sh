#!/bin/bash

wd=`pwd`
ad=${wd}"/assets/genomes"
cd ${ad}
# Download each chromosome individually from UCSC
for chr in {1..22} X Y M; do
    echo "Downloading chr${chr}..."
    wget "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr${chr}.fa.gz"
    gunzip "chr${chr}.fa.gz"
done

cd ${wd}

exit
