# Reference Genome Index for STAR

This pipeline uses the STAR aligner (https://github.com/alexdobin/STAR) and
accordingly uses and indexed genome.

## Quick setup (for 150 bp read libraries):
```bash
# Run the provided script to download individual chromosomes
bash assets/STARindex/buildIndex.sh
