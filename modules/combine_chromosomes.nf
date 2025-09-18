process COMBINE_CHROMOSOMES {
    label 'standard'
    tag "${experiment}_${assay}_${strand}"
    publishDir "${params.outdir}/genome_wide_bedgraphs", mode: 'copy', pattern: "*.bedGraph"
    
    input:
    tuple val(experiment), val(assay), val(strand)
    
    output:
    tuple val(experiment), val(assay), val(strand), 
          path("${experiment}_${assay}_${strand}.bedGraph"), emit: genome_wide_bedgraph
    
    when:
    // Only process DMS experiments (not controls)
    experiment in (params.dms_experiments + params.stress_dms_experiments)
    
    script:
    """
    set -euo pipefail
    echo "Combining chromosomes for ${experiment}_${assay}_${strand}"
    
    # Check if normalized rates directory exists
    NORM_DIR="${params.outdir}/normalized_rates"
    if [[ ! -d "\$NORM_DIR" ]]; then
        echo "ERROR: Normalized rates directory not found: \$NORM_DIR"
        exit 1
    fi
    
    # Initialize output file
    > ${experiment}_${assay}_${strand}.bedGraph
    
    # Track files found and missing
    files_found=0
    files_missing=0
    
    echo "Processing chromosomes in order..."
    
    # Process each chromosome in order
    for chr in ${params.chromosomes.join(' ')}; do
        input_file="\$NORM_DIR/${experiment}_${assay}_norm_\${chr}_${strand}.bedGraph"
        
        if [[ -f "\$input_file" ]]; then
            if [[ -s "\$input_file" ]]; then
                echo "  Adding \$chr: `wc -l < "\$input_file"` lines"
                cat "\$input_file" >> ${experiment}_${assay}_${strand}.bedGraph
                files_found=\$((files_found + 1))
            else
                echo "  WARNING: \$chr file is empty: \$input_file"
                files_missing=\$((files_missing + 1))
            fi
        else
            echo "  WARNING: \$chr file not found: \$input_file"
            files_missing=\$((files_missing + 1))
        fi
    done
    
    # Validate final output
    if [[ -s ${experiment}_${assay}_${strand}.bedGraph ]]; then
        total_lines=\$(wc -l < ${experiment}_${assay}_${strand}.bedGraph)
        echo "Successfully combined \$files_found chromosomes into ${experiment}_${assay}_${strand}.bedGraph"
        echo "Total lines in combined file: \$total_lines"
        
        # Sort the combined bedGraph file by chromosome and position
        echo "Sorting combined bedGraph file..."
        sort -k1,1 -k2,2n ${experiment}_${assay}_${strand}.bedGraph > ${experiment}_${assay}_${strand}.sorted.bedGraph
        mv ${experiment}_${assay}_${strand}.sorted.bedGraph ${experiment}_${assay}_${strand}.bedGraph
        
        echo "Final sorted file: `wc -l < ${experiment}_${assay}_${strand}.bedGraph` lines"
    else
        echo "WARNING: No data found for ${experiment}_${assay}_${strand}"
        # Create empty file to avoid pipeline failure
        touch ${experiment}_${assay}_${strand}.bedGraph
    fi
    
    # Report summary
    echo "=== Combination Summary ==="
    echo "Experiment: ${experiment}"
    echo "Assay: ${assay}"
    echo "Strand: ${strand}"
    echo "Files found: \$files_found"
    echo "Files missing: \$files_missing"
    echo "Output file: ${experiment}_${assay}_${strand}.bedGraph"
    
    echo "Completed chromosome combination for ${experiment}_${assay}_${strand}"
    """
    
    stub:
    """
    touch ${experiment}_${assay}_${strand}.bedGraph
    """
}
