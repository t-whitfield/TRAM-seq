#!/bin/bash

wd=`pwd`
ad=${wd}"/assets/gtf"
cd ${ad}

echo "================================================"
echo "Preparing homo sapiens ENSEMBL release 104 GTF"
echo "================================================"

# Decompress annotations.
zmore Homo_sapiens.GRCh38.104.chr.gtf.gz | grep -v ^\# | awk -F "\t" '{ if($1 ~ /^[1-9XYM]/) print "chr"$0 }' | perl -pe 's/^chrMT/chrM/' > Homo_sapiens.GRCh38.104.canonical.gtf

echo ""
echo "================================================"
echo "GTF ready."
echo "================================================"

cd ${wd}

exit
