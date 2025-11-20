process MAP_READS {
    label 'high_memory'
    tag "${experiment}"
    publishDir "${params.outdir}/mapping_stats", mode: 'copy', pattern: "*_stats.tsv"
    publishDir "${params.outdir}/trim_reports", mode: 'copy', pattern: "*_trimreport.txt"
    publishDir "${params.outdir}/dedup_reports", mode: 'copy', pattern: "*_dedupreport.txt"
    
    input:
    val experiment
    
    output:
    tuple val(experiment), path("${experiment}.bam"), path("${experiment}.bai"), emit: bam_files
    path "${experiment}_full_stats.tsv", emit: full_stats
    path "${experiment}_primary_stats.tsv", emit: primary_stats
    path "${experiment}_trimreport.txt", emit: trim_report
    path "${experiment}_dedupreport.txt", emit: dedup_report
    
    script:
    """
    set -euo pipefail
    echo "Processing experiment: ${experiment}"
    
    # Check if input FASTQ files exist
    if [[ ! -f "${params.fastqpath}/${experiment}_R1.fastq.gz" ]]; then
        echo "ERROR: R1 FASTQ file not found: ${params.fastqpath}/${experiment}_R1.fastq.gz"
        exit 1
    fi
    
    if [[ ! -f "${params.fastqpath}/${experiment}_R2.fastq.gz" ]]; then
        echo "ERROR: R2 FASTQ file not found: ${params.fastqpath}/${experiment}_R2.fastq.gz"
        exit 1
    fi
    
    echo "=== Step 1: Trimming with cutadapt ==="
    # Trim adapters, quality filter at Q20, trim 10bp from 5' of R2,
    # minimum length 10bp
    cutadapt \\
        -a ${params.adapter_r1} \\
        -A ${params.adapter_r2} \\
        -u ${params.trim_5p_r1} \\
        -U ${params.trim_5p_r2} \\
        -q ${params.quality_cutoff} \\
        -j ${task.cpus} \\
        -m ${params.min_length} \\
        -o ${experiment}_R1_trimmed.fq.gz \\
        -p ${experiment}_R2_trimmed.fq.gz \\
        ${params.fastqpath}/${experiment}_R1.fastq.gz \\
        ${params.fastqpath}/${experiment}_R2.fastq.gz \\
        > ${experiment}_trimreport.txt
    
    if [[ ! -s ${experiment}_R1_trimmed.fq.gz ]]; then
        echo "ERROR: Trimming failed - no output generated"
        exit 1
    fi
    
    echo "=== Step 2: Deduplication with clumpify.sh ==="
    # Check if BBMap path exists
    if [[ ! -d "${params.bbmappath}" ]]; then
        echo "ERROR: BBMap directory not found: ${params.bbmappath}"
        exit 1
    fi
    
    if [[ ! -f "${params.bbmappath}/clumpify.sh" ]]; then
        echo "ERROR: clumpify.sh not found in: ${params.bbmappath}"
        exit 1
    fi
    
    # Deduplicate using clumpify from bbmap package
    ${params.bbmappath}/clumpify.sh \\
        in1=${experiment}_R1_trimmed.fq.gz \\
        in2=${experiment}_R2_trimmed.fq.gz \\
        out1=${experiment}_R1_dedup.fq.gz \\
        out2=${experiment}_R2_dedup.fq.gz \\
        dedupe=t \\
        optical=f \\
        subs=${params.dedup_subs} \\
        -Xmx${params.dedup_memory} \\
        2> ${experiment}_dedupreport.txt
    
    if [[ ! -s ${experiment}_R1_dedup.fq.gz ]]; then
        echo "ERROR: Deduplication failed - no output generated"
        exit 1
    fi
    
    # Clean up intermediate trimmed files
    rm ${experiment}_R1_trimmed.fq.gz ${experiment}_R2_trimmed.fq.gz
    
    echo "=== Step 3: Mapping with STAR ==="
    # Check if STAR genome index exists
    if [[ ! -d "${params.star_index}" ]]; then
        echo "ERROR: STAR genome index not found: ${params.star_index}"
        exit 1
    fi
    
    # Map using STAR with basic settings
    STAR \\
        --genomeLoad NoSharedMemory \\
        --genomeDir ${params.star_index} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --readFilesIn ${experiment}_R1_dedup.fq.gz ${experiment}_R2_dedup.fq.gz \\
        --outFileNamePrefix ${experiment}_ \\
        --outSAMtype BAM Unsorted
    
    if [[ ! -f ${experiment}_Aligned.out.bam ]]; then
        echo "ERROR: STAR mapping failed - no BAM file generated"
        exit 1
    fi
    
    # Clean up intermediate dedup files
    rm ${experiment}_R1_dedup.fq.gz ${experiment}_R2_dedup.fq.gz
    
    echo "=== Step 4: Sorting BAM file ==="
    # Sort the BAM file
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${experiment}_full.bam \\
        ${experiment}_Aligned.out.bam
    
    # Clean up unsorted BAM
    rm ${experiment}_Aligned.out.bam
    
    echo "=== Step 5: Filtering for primary alignments ==="
    # Remove secondary (flag 256) and unmapped (flag 4) reads
    # Keep only primary and unique alignments,
    #  filter for MAPQ>${params.mapping_quality}
    samtools view \\
        -b \\
        -F 260 \\
        -@ ${task.cpus} \\
        -q ${params.mapping_quality} \\
        ${experiment}_full.bam \\
        > ${experiment}.bam
    
    if [[ ! -s ${experiment}.bam ]]; then
        echo "ERROR: Filtering failed - no output generated"
        exit 1
    fi
    
    echo "=== Step 6: Generating alignment statistics ==="
    # Generate flagstat for full BAM
    samtools flagstat \\
        -O tsv \\
        -@ ${task.cpus} \\
        ${experiment}_full.bam \\
        > ${experiment}_full_stats.tsv
    
    # Generate flagstat for primary BAM
    samtools flagstat \\
        -O tsv \\
        -@ ${task.cpus} \\
        ${experiment}.bam \\
        > ${experiment}_primary_stats.tsv
    
    echo "=== Step 7: Indexing BAM files ==="
    # Index primary BAM file (this is what will be used downstream)
    samtools index \\
        -@ ${task.cpus} \\
        ${experiment}.bam \\
        ${experiment}.bai
    
    # Clean up full sorted BAM (we only need the primary filtered version)
    rm ${experiment}_full.bam
    
    # Check final outputs
    if [[ ! -s ${experiment}.bam ]] || [[ ! -s ${experiment}.bai ]]; then
        echo "ERROR: Final BAM or index file missing or empty"
        exit 1
    fi
    
    echo "=== Mapping Summary ==="
    echo "Experiment: ${experiment}"
    echo "Primary BAM size: `du -h ${experiment}.bam | cut -f1`"
    echo "Trim report: ${experiment}_trimreport.txt"
    echo "Dedup report: ${experiment}_dedupreport.txt"
    echo "Stats: ${experiment}_primary_stats.tsv"
    echo "Mapping completed successfully"
    """
    
    stub:
    """
    touch ${experiment}.bam
    touch ${experiment}.bai
    touch ${experiment}_full_stats.tsv
    touch ${experiment}_primary_stats.tsv
    touch ${experiment}_trimreport.txt
    touch ${experiment}_dedupreport.txt
    """
}
