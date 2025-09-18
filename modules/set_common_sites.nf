process SET_COMMON_SITES {
    label 'standard'
    tag "${experiment}_${assay}_${strand}"
    publishDir "${params.outdir}/common_filtered", mode: 'copy', pattern: "*.bedGraph"
    
    input:
    tuple val(experiment), val(assay), val(strand)
    
    output:
    tuple val(experiment), val(assay), val(strand), 
          path("${experiment}_${assay}_${strand}.bedGraph"), emit: common_filtered_bedgraph
    
    when:
    // Only process DMS experiments (not controls)
    experiment in (params.dms_experiments + params.stress_dms_experiments)
    
    script:
    """
    set -euo pipefail
    echo "Setting common sites for ${experiment}_${assay}_${strand}"
    
    # Check if genome-wide bedGraphs directory exists
    BEDGRAPH_DIR="${params.outdir}/genome_wide_bedgraphs"
    if [[ ! -d "\$BEDGRAPH_DIR" ]]; then
        echo "ERROR: Genome-wide bedGraphs directory not found: \$BEDGRAPH_DIR"
        exit 1
    fi
    
    # Check if common sites directory exists
    COMMON_SITES_DIR="${params.outdir}/common_sites"
    if [[ ! -d "\$COMMON_SITES_DIR" ]]; then
        echo "ERROR: Common sites directory not found: \$COMMON_SITES_DIR"
        exit 1
    fi
    
    # Define file paths
    input_file="\$BEDGRAPH_DIR/${experiment}_${assay}_${strand}.bedGraph"
    common_sites_file="\$COMMON_SITES_DIR/commonSites_coverage_${strand}.bedGraph"
    
    # Check if input file exists
    if [[ ! -f "\$input_file" ]]; then
        echo "ERROR: Input bedGraph file not found: \$input_file"
        exit 1
    fi
    
    if [[ ! -s "\$input_file" ]]; then
        echo "WARNING: Input bedGraph file is empty: \$input_file"
        touch ${experiment}_${assay}_${strand}.bedGraph
        exit 0
    fi
    
    # Check if common sites file exists and has content
    if [[ ! -f "\$common_sites_file" ]]; then
        echo "ERROR: Common sites file not found: \$common_sites_file"
        exit 1
    fi
    
    if [[ ! -s "\$common_sites_file" ]]; then
        echo "WARNING: Common sites file is empty: \$common_sites_file"
        touch ${experiment}_${assay}_${strand}.bedGraph
        exit 0
    fi
    
    # Check if bedtools is available
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: bedtools is not available in PATH"
        exit 1
    fi
    
    # Report input file statistics
    input_lines=`wc -l < "\$input_file"`
    common_lines=`wc -l < "\$common_sites_file"`
    echo "Input file: \$input_lines lines"
    echo "Common sites: \$common_lines lines"
    
    # Create temporary file (mimicking the original script behavior)
    cp "\$input_file" temp_${experiment}_${assay}_${strand}.bedGraph
    
    echo "Intersecting with common sites..."
    # Intersect with common sites to retain only common positions
    bedtools intersect \\
        -wb -wa \\
        -a temp_${experiment}_${assay}_${strand}.bedGraph \\
        -b "\$common_sites_file" | \\
        awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4}' > ${experiment}_${assay}_${strand}.bedGraph
    
    # Clean up temporary file
    rm temp_${experiment}_${assay}_${strand}.bedGraph
    
    # Validate output
    if [[ -s ${experiment}_${assay}_${strand}.bedGraph ]]; then
        output_lines=`wc -l < ${experiment}_${assay}_${strand}.bedGraph`
        echo "Successfully filtered to common sites: \$output_lines lines"
        
        # Sort the output file
        echo "Sorting filtered bedGraph file..."
        sort -k1,1 -k2,2n ${experiment}_${assay}_${strand}.bedGraph > ${experiment}_${assay}_${strand}.sorted.bedGraph
        mv ${experiment}_${assay}_${strand}.sorted.bedGraph ${experiment}_${assay}_${strand}.bedGraph
        
        final_lines=`wc -l < ${experiment}_${assay}_${strand}.bedGraph`
        echo "Final sorted file: \$final_lines lines"
    else
        echo "WARNING: No common sites found for ${experiment}_${assay}_${strand}"
        touch ${experiment}_${assay}_${strand}.bedGraph
    fi
    
    # Report summary
    echo "=== Common Sites Filtering Summary ==="
    echo "Experiment: ${experiment}"
    echo "Assay: ${assay}"
    echo "Strand: ${strand}"
    echo "Input lines: \$input_lines"
    echo "Common sites: \$common_lines"
    echo "Output lines: `wc -l < ${experiment}_${assay}_${strand}.bedGraph`"
    echo "Output file: ${experiment}_${assay}_${strand}.bedGraph"
    
    echo "Completed common sites filtering for ${experiment}_${assay}_${strand}"
    """
    
    stub:
    """
    touch ${experiment}_${assay}_${strand}.bedGraph
    """
}