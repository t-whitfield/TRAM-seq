process COMPUTE_ELEMENT_STATS_SPLICED {
    label 'standard'
    tag "element_stats_${reference}"
    publishDir "${params.outdir}/element_statistics", mode: 'copy', pattern: "*.{txt,csv,bed}"
    
    input:
    val reference
    
    output:
    path "coverage*AllGenes.txt", emit: coverage_stats, optional: true
    path "mmRate*AllGenes.txt", emit: mmrate_stats, optional: true
    path "processedExperiments*.csv", emit: processed_stats, optional: true
    path "*_collapseGenes.bed", emit: collapsed_genes, optional: true
    
    script:
    def reference_tag = params.reference_tags[reference]
    def dms_count=params.dms_experiments.size()
    def stress_dms_count=params.stress_dms_experiments.size()
    def dms_list=params.dms_experiments
    def stress_dms_list=params.stress_dms_experiments
    """
    set -euo pipefail
    echo "Computing element statistics for reference: ${reference}.bed"
    echo "Using reference tag: ${reference_tag}"
    
    # Check if required directories exist
    COMMON_FILTERED_DIR="${params.outdir}/common_filtered"
    if [[ ! -d "\$COMMON_FILTERED_DIR" ]]; then
        echo "ERROR: Common filtered directory not found: \$COMMON_FILTERED_DIR"
        exit 1
    fi
    
    ANNOT_DIR="${params.annotpath}"
    if [[ ! -d "\$ANNOT_DIR" ]]; then
        echo "ERROR: Annotations directory not found: \$ANNOT_DIR"
        exit 1
    fi

    HEADER_DIR="${params.headpath}"
    if [[ ! -d "\$HEADER_DIR" ]]; then
        echo "ERROR: Header directory not found: \$HEADER_DIR"
        exit 1
    fi
        
    # Check if required annotation files exist
    for strand in pos neg; do
        annot_file="\$ANNOT_DIR/${reference}_\${strand}.bed"
        if [[ ! -f "\$annot_file" ]]; then
            echo "ERROR: Annotation file not found: \$annot_file"
            exit 1
        fi
        echo "Found annotation file: \$annot_file (`wc -l < "\$annot_file"` lines)"
    done
    
    # Check if required python scripts exist in bin directory
    for script in updateBed.py process_replicate_signal.py; do
    	script_path="${projectDir}/bin/\$script"
    	if [[ ! -f "\$script_path" ]]; then
           echo "ERROR: Required script not found: \$script_path"
           exit 1
    	fi
    	echo "Found script: \$script_path"
    done
    
    # Check if header files exist
    for assay in coverage mmRate; do
        header_file="\$HEADER_DIR/\${assay}Head6.txt"
        if [[ ! -f "\$header_file" ]]; then
            echo "ERROR: Header file not found: \$header_file"
            exit 1
        fi
        echo "Found header file: \$header_file"
    done
    
    # Define DMS experiments
    dms_experiments="${params.dms_experiments.join(' ')} ${params.stress_dms_experiments.join(' ')}"
    echo "Processing DMS experiments: \$dms_experiments"
    
    # Check if bedtools is available
    if ! command -v bedtools &> /dev/null; then
        echo "ERROR: bedtools is not available in PATH"
        exit 1
    fi
    
    # Check if python is available
    if ! command -v python &> /dev/null; then
        echo "ERROR: python is not available in PATH"
        exit 1
    fi
    
    echo "=== Step 1: Intersecting with annotations and grouping by regions ==="
    
    # First loop: collect per-nucleotide coverage and mismatch rates for each region
    for treat in \$dms_experiments; do
        echo "Processing experiment: \$treat"
        
        for assay in coverage mmRate; do
            echo "  Processing assay: \$assay"
            
            for strand in pos neg; do
                echo "    Processing strand: \$strand"
                
                # Check if input file exists
                input_file="\$COMMON_FILTERED_DIR/\${treat}_\${assay}_\${strand}.bedGraph"
                if [[ ! -f "\$input_file" ]]; then
                    echo "    ERROR: Input file not found: \$input_file"
                    exit 1
                fi
                
                if [[ ! -s "\$input_file" ]]; then
                    echo "    WARNING: Input file is empty: \$input_file"
                    touch \${treat}_\${assay}_\${strand}.bed
                    touch \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.bed
                    continue
                fi
                
                input_lines=`wc -l < "\$input_file"`
                echo "    Input file: \$input_lines lines"
                
                # Intersect with annotations
                bedtools intersect \\
                    -wb -wa \\
                    -a "\$ANNOT_DIR/${reference}_\${strand}.bed" \\
                    -b "\$input_file" > \${treat}_\${assay}_\${strand}.bed
                
                if [[ ! -s \${treat}_\${assay}_\${strand}.bed ]]; then
                    echo "    WARNING: No intersections found for \$treat \$assay \$strand"
                    touch \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.bed
                    continue
                fi
                
                intersect_lines=`wc -l < \${treat}_\${assay}_\${strand}.bed`
                echo "    Intersections: \$intersect_lines lines"
                
                # Group by regions
                bedtools groupby \\
                    -i \${treat}_\${assay}_\${strand}.bed \\
                    -g 1,2,3,4,5,6 \\
                    -c 10 \\
                    -o collapse > \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.txt
                
                # Reformat output
                awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4,\$7,\$6}' \\
                    \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.txt > \\
                    \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.bed
                
                # Clean up intermediate files
                rm \${treat}_\${assay}_\${strand}.bed \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.txt
                
                collapse_lines=`wc -l < \${treat}_\${assay}_\${strand}_${reference_tag}_collapse.bed`
                echo "    Collapsed regions: \$collapse_lines lines"
            done
            
            # Combine positive and negative strands
            echo "  Combining strands for \$treat \$assay"
            cat \${treat}_\${assay}_pos_${reference_tag}_collapse.bed \\
                \${treat}_\${assay}_neg_${reference_tag}_collapse.bed > \\
                \${treat}_\${assay}_${reference_tag}_collapsed.bed
            
            # Sort combined file
            sort -k 1,1 -k4,4 -k2,2n \\
                \${treat}_\${assay}_${reference_tag}_collapsed.bed > \\
                \${treat}_\${assay}_${reference_tag}_collapse.bed
            
            # Clean up intermediate files
            rm \${treat}_\${assay}_pos_${reference_tag}_collapse.bed \\
               \${treat}_\${assay}_neg_${reference_tag}_collapse.bed \\
               \${treat}_\${assay}_${reference_tag}_collapsed.bed
            
            final_lines=`wc -l < \${treat}_\${assay}_${reference_tag}_collapse.bed`
            echo "  Final collapsed file: \$final_lines lines"
        done
    done
    
    echo "=== Step 2: Aggregating regions into transcripts ==="
    
    # Second loop: aggregate regions into transcripts
    for treat in \$dms_experiments; do
        echo "Aggregating transcripts for: \$treat"
        
        for assay in coverage mmRate; do
            echo "  Processing assay: \$assay"
            
            # Create temporary files for each column
            bedtools groupby -i \${treat}_\${assay}_${reference_tag}_collapse.bed -g 4 -c 1 -o first | awk '{print \$2}' > \${treat}_\${assay}_${reference_tag}_collapseGenes_chr.tmp
            bedtools groupby -i \${treat}_\${assay}_${reference_tag}_collapse.bed -g 4 -c 2 -o min | awk '{print \$2}' > \${treat}_\${assay}_${reference_tag}_collapseGenes_start.tmp
            bedtools groupby -i \${treat}_\${assay}_${reference_tag}_collapse.bed -g 4 -c 3 -o max | awk '{print \$2}' > \${treat}_\${assay}_${reference_tag}_collapseGenes_end.tmp
            bedtools groupby -i \${treat}_\${assay}_${reference_tag}_collapse.bed -g 4 -c 6 -o first | awk '{print \$2}' > \${treat}_\${assay}_${reference_tag}_collapseGenes_strand.tmp
            bedtools groupby -i \${treat}_\${assay}_${reference_tag}_collapse.bed -g 4 -c 5 -o collapse > \${treat}_\${assay}_${reference_tag}_collapseGenes_sig.tmp
            
            # Combine columns using pr
            pr -m -t -J \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes_chr.tmp \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes_start.tmp \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes_end.tmp \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes_sig.tmp \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes_strand.tmp | \\
                awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4,\$5,\$6}' > \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes.tmp
            
            # Sort the result
            sort -k 1,1 -k 2,2n \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes.tmp > \\
                \${treat}_\${assay}_${reference_tag}_collapseGenes.bed
            
            # Clean up temporary files
            rm \${treat}_\${assay}_${reference_tag}_collapseGenes*.tmp
            
            genes_lines=`wc -l < \${treat}_\${assay}_${reference_tag}_collapseGenes.bed`
            echo "  Collapsed genes: \$genes_lines lines"
        done
    done
    
    echo "=== Step 3: Updating coordinates to reflect full transcript lengths ==="
    
    # Third loop: update coordinates using Python script
    for treat in \$dms_experiments; do
        echo "Updating coordinates for: \$treat"
        
        for assay in coverage mmRate; do
            echo "  Processing assay: \$assay"
            
            mv \${treat}_\${assay}_${reference_tag}_collapseGenes.bed \${treat}_\${assay}_${reference_tag}_collapseGenes.tmp
            
            python ${projectDir}/bin/updateBed.py \\
                --ref "\$ANNOT_DIR/${reference}.bed" \\
                --query \${treat}_\${assay}_${reference_tag}_collapseGenes.tmp \\
                --out \${treat}_\${assay}_${reference_tag}_collapseGenes.bed
            
            rm \${treat}_\${assay}_${reference_tag}_collapseGenes.tmp
            
            updated_lines=`wc -l < \${treat}_\${assay}_${reference_tag}_collapseGenes.bed`
            echo "  Updated coordinates: \$updated_lines lines"
        done
    done
    
    echo "=== Checking replicate availability for Steps 4 and 5 ==="

    # Check if we have 3 replicates for two conditions

    echo "Normal condition DMS experiments: ${dms_count}"
    echo "Stress condition DMS experiments: ${stress_dms_count}"

    if [[ ${dms_count} -eq 3 && ${stress_dms_count} -eq 3 ]]; then
        echo "Found 3 replicates for both conditions - proceeding with Steps 4 and 5"
    	# Define shell arrays for experiments
    	dms_array=(${params.dms_experiments.join(' ')})
    	stress_dms_array=(${params.stress_dms_experiments.join(' ')})

        echo "=== Step 4: Creating final output files ==="
    
    # Fourth loop: create final output files for each assay
        for assay in coverage mmRate; do
            echo "Creating final output for assay: \$assay"
        
    # Extract data from each experiment
            awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4,\$6,\$5}' \${dms_array[0]}_\${assay}_${reference_tag}_collapseGenes.bed > \${dms_array[0]}.tmp
            awk '{print \$5}' \${dms_array[1]}_\${assay}_${reference_tag}_collapseGenes.bed > \${dms_array[1]}.tmp
            awk '{print \$5}' \${dms_array[2]}_\${assay}_${reference_tag}_collapseGenes.bed > \${dms_array[2]}.tmp
            awk '{print \$5}' \${stress_dms_array[0]}_\${assay}_${reference_tag}_collapseGenes.bed > \${stress_dms_array[0]}.tmp
            awk '{print \$5}' \${stress_dms_array[1]}_\${assay}_${reference_tag}_collapseGenes.bed > \${stress_dms_array[1]}.tmp
            awk '{print \$5}' \${stress_dms_array[2]}_\${assay}_${reference_tag}_collapseGenes.bed > \${stress_dms_array[2]}.tmp
        
    # Combine all experiments
            pr -m -t -J \${dms_array[0]}.tmp \${dms_array[1]}.tmp \${dms_array[2]}.tmp \${stress_dms_array[0]}.tmp \${stress_dms_array[1]}.tmp \${stress_dms_array[2]}.tmp | \\
                awk 'BEGIN{OFS="\\t"} {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11}' > \\
                \${assay}_${reference_tag}_AllGenes.tmp
        
    # Add header
            cat "\$HEADER_DIR/\${assay}Head6.txt" \${assay}_${reference_tag}_AllGenes.tmp > \${assay}_${reference_tag}_AllGenes.txt
        
    # Clean up temporary files
            rm -f \${dms_array[0]}.tmp \${dms_array[1]}.tmp \${dms_array[2]}.tmp \${stress_dms_array[0]}.tmp \${stress_dms_array[1]}.tmp \${stress_dms_array[2]}.tmp
            rm -f *_\${assay}_${reference_tag}_collapse*.bed
            rm -f \${assay}_${reference_tag}_AllGenes.tmp
        
	    final_lines=`wc -l < \${assay}_${reference_tag}_AllGenes.txt`
            echo "Final \$assay file: \$final_lines lines"
        done
    
        echo "=== Step 5: Processing replicate signals ==="
    
    # Final step: process replicate signals
        python ${projectDir}/bin/process_replicate_signal.py \\
            --covin coverage_${reference_tag}_AllGenes.txt \\
            --mmin mmRate_${reference_tag}_AllGenes.txt \\
            --outfile processedExperiments_${reference_tag}.csv
     
        if [[ -s processedExperiments_${reference_tag}.csv ]]; then
            processed_lines=`wc -l < processedExperiments_${reference_tag}.csv`
            echo "Final processed file: \$processed_lines lines"
        else
	    echo "WARNING: Final processed file is empty or missing"
        fi
    
        echo "=== Element Statistics Summary ==="
        echo "Reference: ${reference}.bed"
        echo "Coverage stats: `wc -l < coverage_${reference_tag}_AllGenes.txt` lines"
        echo "mmRate stats: `wc -l < mmRate_${reference_tag}_AllGenes.txt` lines"
        echo "Processed stats: `wc -l < processedExperiments_${reference_tag}.csv` lines"

    else
        echo "Insufficient replicates for full analysis (need 3 replicates for both conditions)"
        echo "Skipping Steps 4 and 5 - only individual gene files will be available"
        
        # Create empty output files to satisfy the process outputs
        touch coverage_${reference_tag}_AllGenes.txt
        touch mmRate_${reference_tag}_AllGenes.txt  
        touch processedExperiments_${reference_tag}.csv
        
        echo "=== Element Statistics Summary ==="
        echo "Reference: ${reference}.bed"
        echo "Individual gene files created, but no combined analysis performed"
    fi
    
    echo "Completed element statistics computation"
    """
    
    stub:
    """
    def reference_tag = params.reference_tags[reference]
    touch coverage_${reference_tag}_AllGenes.txt
    touch mmRate_${reference_tag}_AllGenes.txt
    touch processedExperiments_${reference_tag}.csv
    """
}
