#!/bin/bash

wd=`pwd`
ad=${wd}"/assets/STARindex"
cd ${ad}

echo "================================================"
echo "Preparing homo sapiens STAR index (optimized for"
echo "150 bp read libraries)"
echo "================================================"

tcores=`nproc`
mcores="32"
cores=$(( ${tcores} < ${mcores} ? ${tcores} : ${mcores} ))

mkdir GRCh38.104.canonical_overhang_150

# Create the STAR index using the available number of threads.
STAR --runMode genomeGenerate --runThreadN ${cores} --genomeDir GRCh38.104.canonical_overhang_150 --genomeFastaFiles ../genomes/hg38.fa --limitGenomeGenerateRAM 24000000000 --sjdbGTFfile ../gtf/Homo_sapiens.GRCh38.104.canonical.gtf --sjdbOverhang 150

echo ""
echo "================================================"
echo "STAR index ready."
echo "================================================"

cd ${wd}

exit
