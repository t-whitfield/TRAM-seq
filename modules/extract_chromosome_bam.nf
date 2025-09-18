process EXTRACT_CHROMOSOME_BAM {
    label 'standard'
    tag "${experiment}_${chromosome}"
    publishDir "${params.outdir}/bam_files", mode: 'copy', pattern: "*.{bam,bai}"
    
    input:
    tuple val(experiment), val(chromosome)
    
    output:
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}.bam"), path("${experiment}_${chromosome}.bai"), emit: bam_files
    
    script:
    """
    set -euo pipefail
    echo "Processing ${experiment} chromosome ${chromosome}"
    echo "Input BAM: ${params.datapath}/${experiment}.bam"
    
    # Check if input BAM file exists
    if [[ ! -f "${params.datapath}/${experiment}.bam" ]]; then
        echo "ERROR: Input BAM file not found: ${params.datapath}/${experiment}.bam"
        exit 1
    fi
    
    # Extract chromosome-specific bam file (output to work directory)
    samtools view "${params.datapath}/${experiment}.bam" "${chromosome}" -b > "${experiment}_${chromosome}.bam"
    
    # Check if extraction was successful
    if [[ ! -s "${experiment}_${chromosome}.bam" ]]; then
        echo "WARNING: No reads found for ${chromosome} in ${experiment}"
        # Create empty BAM file with header
        samtools view -H "${params.datapath}/${experiment}.bam" | samtools view -bS - > "${experiment}_${chromosome}.bam"
    fi
    
    # Index chromosome-specific bam file (output to work directory)
    samtools index "${experiment}_${chromosome}.bam" "${experiment}_${chromosome}.bai"
    
    echo "Completed BAM extraction for ${experiment}_${chromosome}"
    """
    
    stub:
    """
    touch ${experiment}_${chromosome}.bam
    touch ${experiment}_${chromosome}.bai
    """
}