process FIND_COMMON_SITES {
    label 'standard'
    tag "common_sites_${strand}"
    publishDir "${params.outdir}/common_sites", mode: 'copy', pattern: "*.bedGraph"
    
    input:
    val strand
    
    output:
    tuple val(strand), path("commonSites_coverage_${strand}.bedGraph"), emit: common_sites
    
    script:
    """
    set -euo pipefail
    echo "Finding common sites for ${strand} strand across all DMS experiments"
    
    # Check if genome-wide bedGraphs directory exists
    BEDGRAPH_DIR="${params.outdir}/genome_wide_bedgraphs"
    if [[ ! -d "\$BEDGRAPH_DIR" ]]; then
        echo "ERROR: Genome-wide bedGraphs directory not found: \$BEDGRAPH_DIR"
        exit 1
    fi
    
    # Define all DMS experiments
    DMS_EXPERIMENTS="${params.dms_experiments.join(' ')} ${params.stress_dms_experiments.join(' ')}"
    echo "Processing DMS experiments: \$DMS_EXPERIMENTS"
    
    # Convert to array for easier handling
    dms_array=(${params.dms_experiments.join(' ')} ${params.stress_dms_experiments.join(' ')})
    
    # Check if all required input files exist
    echo "Validating input files..."
    for exp in \${dms_array[@]}; do
        input_file="\$BEDGRAPH_DIR/\${exp}_coverage_${strand}.bedGraph"
        if [[ ! -f "\$input_file" ]]; then
            echo "ERROR: Required input file not found: \$input_file"
            exit 1
        fi
        if [[ ! -s "\$input_file" ]]; then
            echo "WARNING: Input file is empty: \$input_file"
        else
            line_count=`wc -l < "\$input_file"`
            echo "  Validated: \$exp (\$line_count lines)"
        fi
    done
    
    # Check if bedtools is available
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: bedtools is not available in PATH"
        exit 1
    fi
    
    # Start with the first two experiments
    first_exp=\${dms_array[0]}
    second_exp=\${dms_array[1]}
    
    echo "Starting intersection with \$first_exp and \$second_exp..."
    bedtools intersect \\
        -wb -wa \\
        -a "\$BEDGRAPH_DIR/\${first_exp}_coverage_${strand}.bedGraph" \\
        -b "\$BEDGRAPH_DIR/\${second_exp}_coverage_${strand}.bedGraph" | \\
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' > temp.bedGraph
    
    # Check if initial intersection produced results
    if [[ ! -s temp.bedGraph ]]; then
        echo "WARNING: No common sites found between \$first_exp and \$second_exp"
        touch commonSites_coverage_${strand}.bedGraph
        exit 0
    fi
    
    initial_sites=`wc -l < temp.bedGraph`
    echo "Initial intersection: \$initial_sites sites"
    
    # Intersect with remaining experiments
    remaining_experiments=(\${dms_array[@]:2})  # Skip first two elements
    
    for exp in \${remaining_experiments[@]}; do
        echo "Intersecting with \$exp..."
        
        bedtools intersect \\
            -wb -wa \\
            -a temp.bedGraph \\
            -b "\$BEDGRAPH_DIR/\${exp}_coverage_${strand}.bedGraph" | \\
            awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' > temp2.bedGraph
        
        # Check if intersection produced results
        if [[ ! -s temp2.bedGraph ]]; then
            echo "WARNING: No common sites after intersecting with \$exp"
            touch commonSites_coverage_${strand}.bedGraph
            rm -f temp.bedGraph temp2.bedGraph
            exit 0
        fi
        
        mv temp2.bedGraph temp.bedGraph
        current_sites=`wc -l < temp.bedGraph`
        echo "  After \$exp: \$current_sites sites remaining"
    done
    
    # Finalize the output
    mv temp.bedGraph commonSites_coverage_${strand}.bedGraph
    
    # Validate final output
    if [[ -s commonSites_coverage_${strand}.bedGraph ]]; then
        final_sites=`wc -l < commonSites_coverage_${strand}.bedGraph`
        echo "Successfully identified \$final_sites common sites across all DMS experiments"
        
        # Sort the final file
        echo "Sorting common sites file..."
        sort -k1,1 -k2,2n commonSites_coverage_${strand}.bedGraph > commonSites_coverage_${strand}.sorted.bedGraph
        mv commonSites_coverage_${strand}.sorted.bedGraph commonSites_coverage_${strand}.bedGraph
        
        echo "Final sorted common sites: `wc -l < commonSites_coverage_${strand}.bedGraph` sites"
    else
        echo "WARNING: No common sites found across all experiments"
        touch commonSites_coverage_${strand}.bedGraph
    fi
    
    # Report summary
    echo "=== Common Sites Summary ==="
    echo "Strand: ${strand}"
    echo "DMS experiments processed: \${#dms_array[@]}"
    echo "Final common sites: `wc -l < commonSites_coverage_${strand}.bedGraph`"
    echo "Output file: commonSites_coverage_${strand}.bedGraph"
    
    echo "Completed common sites identification for ${strand} strand"
    """
    
    stub:
    """
    touch commonSites_coverage_${strand}.bedGraph
    """
}
