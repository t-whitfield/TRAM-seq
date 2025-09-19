# TRAM-seq analysis pipeline

A nextflow pipeline for analyzing TRAM-seq (transcriptome-wide accessibility mapping by sequencing) data to quantify RNA structure and changes across different experimental conditions.

## Overview

This pipeline processes TRAM-seq sequencing data to:
- Generate and process pileup files for mutation detection
- Filter coverage and normalize DMS modification rates
- Compute element-level statistics for genomic regions and transcripts

## Pipeline workflow

### Main processing flow

```mermaid
flowchart TD
    subgraph "Input Data"
        A1[Control BAMs<br/>C1, C2, C3]
        A2[DMS BAMs<br/>D1, D2, D3]
        A3[Stress Control BAMs<br/>A1, A2, A3]
        A4[Stress DMS BAMs<br/>AD1, AD2, AD3]
    end
    
    subgraph "Chromosome Processing"
        B[EXTRACT_CHROMOSOME_BAM<br/>Split by chromosome]
        C[GENERATE_PILEUP<br/>Mutation detection]
        D[PROCESS_PILEUP<br/>Strand separation]
    end
    
    subgraph "Quality Control & Normalization"
        E[FILTER_COVERAGE<br/>Remove low coverage sites]
        F[NORMALIZE_RATES<br/>DMS vs Control normalization]
    end
    
    subgraph "Genome Assembly"
        G[COMBINE_CHROMOSOMES<br/>Create genome-wide tracks]
        H[FIND_COMMON_SITES<br/>Identify shared sites]
        I[SET_COMMON_SITES<br/>Apply common filter]
    end
    
    subgraph "Statistical Analysis"
        J[COMPUTE_ELEMENT_STATS<br/>Basic regions: 3'UTR, 5'UTR, CDS]
        K[COMPUTE_ELEMENT_STATS_SPLICED<br/>Transcript-level analysis]
    end
    
    A1 & A2 & A3 & A4 --> B
    B --> C
    C --> D
    D --> E
    E --> F
    F --> G
    G --> H
    H --> I
    I --> J
    I --> K
    
    style A1 fill:#ffcdd2
    style A2 fill:#f8bbd9
    style A3 fill:#e1bee7
    style A4 fill:#d1c4e9
    style J fill:#c8e6c9
    style K fill:#c8e6c9
```

## Quick start
### Clone the repository
```
git clone https://github.com/whitehead/TRAM-seq.git
cd TRAM-seq
```

### Set up reference genome files
```
bash assets/genomes/download_chromosomes.sh
```

### Run the pipeline
```
nextflow run main.nf -params-file params.yaml
```

## Requirements
### Software dependencies

nextflow (>=22.10.0)\
Java (>=11)\
samtools\
bcftools\
bedtools\
python (>=3.7)

## Setup

### 1. Reference genome files

The pipeline requires individual chromosome FASTA files:
```
bash assets/genomes/download_chromosomes.sh
```

### 2. Input data structure

Organize your BAM files in the data directory.
