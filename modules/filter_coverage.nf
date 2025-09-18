process FILTER_COVERAGE {
    label 'standard'
    tag "${experiment}_${chromosome}"
    publishDir "${params.outdir}/filtered_coverage", mode: 'copy', pattern: "*.bedGraph"
    publishDir "${params.outdir}/compressed_pileup", mode: 'copy', pattern: "*.pileup.gz"
    
    input:
    tuple val(experiment), val(chromosome), path(pos_pileup), path(neg_pileup)
    
    output:
    tuple val(experiment), val(chromosome), path("${experiment}_coverage_${chromosome}_pos.bedGraph"), emit: coverage_pos_bedgraph
    tuple val(experiment), val(chromosome), path("${experiment}_coverage_${chromosome}_neg.bedGraph"), emit: coverage_neg_bedgraph
    tuple val(experiment), val(chromosome), path("${experiment}_mmRate_${chromosome}_pos.bedGraph"), emit: mmrate_pos_bedgraph
    tuple val(experiment), val(chromosome), path("${experiment}_mmRate_${chromosome}_neg.bedGraph"), emit: mmrate_neg_bedgraph
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_pos.pileup.gz"), emit: compressed_pos_pileup
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_neg.pileup.gz"), emit: compressed_neg_pileup
    
    script:
    """
    set -euo pipefail
    echo "Filtering coverage for ${experiment}_${chromosome} with threshold ${params.coverage_threshold}"
    
    # Check if input pileup files exist and have content
    if [[ ! -s "${pos_pileup}" ]]; then
        echo "ERROR: Positive strand pileup file is empty or missing: ${pos_pileup}"
        exit 1
    fi
    
    if [[ ! -s "${neg_pileup}" ]]; then
        echo "ERROR: Negative strand pileup file is empty or missing: ${neg_pileup}"
        exit 1
    fi
    
    echo "Processing positive strand coverage..."
    # Filter coverage for positive strand (column 4 is coverage)
    awk -v cut=${params.coverage_threshold} '{ if (\$4 >= cut) {OFS="\\t"; print \$1,\$2-1,\$2,\$4}}' "${pos_pileup}" | \\
        grep -v "pos" > "${experiment}_coverage_${chromosome}_pos.bedGraph"
    
    echo "Processing negative strand coverage..."
    # Filter coverage for negative strand (column 4 is coverage)
    awk -v cut=${params.coverage_threshold} '{ if (\$4 >= cut) {OFS="\\t"; print \$1,\$2-1,\$2,\$4}}' "${neg_pileup}" | \\
        grep -v "pos" > "${experiment}_coverage_${chromosome}_neg.bedGraph"
    
    echo "Processing positive strand mismatch rates..."
    # Filter mismatch rates for positive strand (column 7 is mismatch rate)
    awk -v cut=${params.coverage_threshold} '{ if (\$4 >= cut) {OFS="\\t"; print \$1,\$2-1,\$2,\$7}}' "${pos_pileup}" | \\
        grep -v "pos" > "${experiment}_mmRate_${chromosome}_pos.bedGraph"
    
    echo "Processing negative strand mismatch rates..."
    # Filter mismatch rates for negative strand (column 7 is mismatch rate)
    awk -v cut=${params.coverage_threshold} '{ if (\$4 >= cut) {OFS="\\t"; print \$1,\$2-1,\$2,\$7}}' "${neg_pileup}" | \\
        grep -v "pos" > "${experiment}_mmRate_${chromosome}_neg.bedGraph"
    
    # Compress original pileup files
    echo "Compressing pileup files..."
    gzip -c "${pos_pileup}" > "${experiment}_${chromosome}_pos.pileup.gz"
    gzip -c "${neg_pileup}" > "${experiment}_${chromosome}_neg.pileup.gz"
    
    # Validate output files were created
    for file in "${experiment}_coverage_${chromosome}_pos.bedGraph" \\
                "${experiment}_coverage_${chromosome}_neg.bedGraph" \\
                "${experiment}_mmRate_${chromosome}_pos.bedGraph" \\
                "${experiment}_mmRate_${chromosome}_neg.bedGraph"; do
        if [[ ! -f \$file ]]; then
            echo "ERROR: Failed to create output file: \$file"
            exit 1
        fi
        echo "Created: \$file (`wc -l < \$file` lines)"
    done
    
    echo "Completed coverage filtering for ${experiment}_${chromosome}"
    """
    
    stub:
    """
    touch ${experiment}_coverage_${chromosome}_pos.bedGraph
    touch ${experiment}_coverage_${chromosome}_neg.bedGraph
    touch ${experiment}_mmRate_${chromosome}_pos.bedGraph
    touch ${experiment}_mmRate_${chromosome}_neg.bedGraph
    touch ${experiment}_${chromosome}_pos.pileup.gz
    touch ${experiment}_${chromosome}_neg.pileup.gz
    """
}