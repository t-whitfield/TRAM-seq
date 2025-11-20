#!/bin/bash

wd=`pwd`
ad=${wd}"/assets/genomes"
cd ${ad}
# Download each chromosome individually from UCSC.
for chr in {1..22} X Y M; do
    echo "Downloading chr${chr}..."
    wget "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr${chr}.fa.gz"
    gunzip "chr${chr}.fa.gz"
done

# Concatenate these fasta files into a single file for the
# assembly, including canonical chromosomes only, no unplaced contigs.
cat chr{1..22}.fa chrX.fa chrY.fa chrM.fa > hg38.fa

# Index the whole-genome assembly.
samtools faidx hg38.fa

cd ${wd}

exit
