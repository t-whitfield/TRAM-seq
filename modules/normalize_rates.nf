process NORMALIZE_RATES {
    label 'standard'
    tag "${dms_experiment}_vs_${control_experiment}_${chromosome}_${strand}"
    publishDir "${params.outdir}/normalized_rates", mode: 'copy', pattern: "*_norm_*.bedGraph"
    
    input:
    tuple val(dms_experiment), val(control_experiment), val(chromosome), val(strand)
    
    output:
    tuple val(dms_experiment), val(chromosome), val(strand), 
          path("${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph"), emit: normalized_coverage
    tuple val(dms_experiment), val(chromosome), val(strand), 
          path("${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph"), emit: normalized_mmrate
    
    script:
    """
    set -euo pipefail
    echo "Normalizing ${dms_experiment} vs ${control_experiment} on ${chromosome} ${strand} strand"
    
    # Define file paths based on the expected output from FILTER_COVERAGE
    DMS_COVERAGE="${params.outdir}/filtered_coverage/${dms_experiment}_coverage_${chromosome}_${strand}.bedGraph"
    CTRL_COVERAGE="${params.outdir}/filtered_coverage/${control_experiment}_coverage_${chromosome}_${strand}.bedGraph"
    DMS_MMRATE="${params.outdir}/filtered_coverage/${dms_experiment}_mmRate_${chromosome}_${strand}.bedGraph"
    CTRL_MMRATE="${params.outdir}/filtered_coverage/${control_experiment}_mmRate_${chromosome}_${strand}.bedGraph"
    
    # Check if all required input files exist
    echo "Checking input files..."
    for file in "\$DMS_COVERAGE" "\$CTRL_COVERAGE" "\$DMS_MMRATE" "\$CTRL_MMRATE"; do
        if [[ ! -f "\$file" ]]; then
            echo "ERROR: Required input file not found: \$file"
            exit 1
        fi
        if [[ ! -s "\$file" ]]; then
            echo "WARNING: Input file is empty: \$file"
        else
            echo "Input file validated: \$file (`wc -l < "\$file"` lines)"
        fi
    done
    
    # Check if bedtools is available
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: bedtools is not available in PATH"
        exit 1
    fi
    
    echo "Processing coverage normalization..."
    # Intersect coverage files and keep DMS coverage values
    bedtools intersect \\
        -wb -wa \\
        -a "\$DMS_COVERAGE" \\
        -b "\$CTRL_COVERAGE" | \\
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' > "${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph"
    
    # Validate coverage normalization output
    if [[ ! -s "${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph" ]]; then
        echo "WARNING: No overlapping coverage regions found between ${dms_experiment} and ${control_experiment}"
        touch "${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph"
    else
        echo "Coverage normalization completed: `wc -l < ${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph` regions"
    fi
    
    echo "Processing mismatch rate normalization..."
    # Intersect mismatch rate files and subtract control from DMS (max with 0)
    bedtools intersect \\
        -wb -wa \\
        -a "\$DMS_MMRATE" \\
        -b "\$CTRL_MMRATE" | \\
        awk 'BEGIN{OFS="\\t"} {
            max = 0
            if (\$4 > \$8) max = \$4 - \$8
            print \$1, \$2, \$3, max
        }' > "${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph"
    
    # Validate mismatch rate normalization output
    if [[ ! -s "${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph" ]]; then
        echo "WARNING: No overlapping mismatch rate regions found between ${dms_experiment} and ${control_experiment}"
        touch "${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph"
    else
        echo "Mismatch rate normalization completed: `wc -l < ${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph` regions"
    fi
    
    # Report statistics
    echo "=== Normalization Summary ==="
    echo "DMS experiment: ${dms_experiment}"
    echo "Control experiment: ${control_experiment}"
    echo "Chromosome: ${chromosome}"
    echo "Strand: ${strand}"
    echo "Normalized coverage regions: `wc -l < ${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph`"
    echo "Normalized mismatch rate regions: `wc -l < ${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph`"
    
    echo "Completed normalization for ${dms_experiment}_${chromosome}_${strand}"
    """
    
    stub:
    """
    touch ${dms_experiment}_coverage_norm_${chromosome}_${strand}.bedGraph
    touch ${dms_experiment}_mmRate_norm_${chromosome}_${strand}.bedGraph
    """
}
