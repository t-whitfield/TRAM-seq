process PROCESS_PILEUP {
    label 'standard'
    tag "${experiment}_${chromosome}"
    publishDir "${params.outdir}/processed_pileup", mode: 'copy'
    
    input:
    tuple val(experiment), val(chromosome), path(pileup_file)
    
    output:
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_pos.pileup"), emit: pos_pileup
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_neg.pileup"), emit: neg_pileup
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_vars.pileup.gz"), emit: vars_pileup
    
    script:
    """
    set -euo pipefail
    echo "Processing pileup for ${experiment}_${chromosome}"

    # Check if required python scripts exist in bin directory
    script_path="${projectDir}/bin/process_pileup_strand_select.py"
    if [[ ! -f "\$script_path" ]]; then
       echo "ERROR: Required script not found: \$script_path"
       exit 1
    fi
    echo "Found script: \$script_path"
 
    # Check if Python script is available in PATH
    if ! command -v process_pileup_strand_select.py &> /dev/null; then
        echo "ERROR: Python script not found in PATH: process_pileup_strand_select.py"
        exit 1
    fi
        
    # Check if input pileup file has content (skip header if present)
    content_lines=\$(tail -n +2 "${pileup_file}" | wc -l)
    if [[ \$content_lines -eq 0 ]]; then
        echo "WARNING: Empty input pileup file for ${experiment}_${chromosome}"
        # Create empty output files with headers
        echo -e "ID\\tCHROM\\tPOS\\tREF\\tTYPE\\tALT\\tIMF\\tDP\\tAD" > "${experiment}_${chromosome}_pos.pileup"
        echo -e "ID\\tCHROM\\tPOS\\tREF\\tTYPE\\tALT\\tIMF\\tDP\\tAD" > "${experiment}_${chromosome}_neg.pileup"
        echo -e "ID\\tCHROM\\tPOS\\tREF\\tTYPE\\tALT\\tIMF\\tDP\\tAD" > "${experiment}_${chromosome}_vars.pileup"
    else
        # Process raw pileup file to get ratios of interest
        echo "Running Python processing script..."
        python ${projectDir}/bin/process_pileup_strand_select.py \\
            --pileup "${pileup_file}" \\
            --outfile1 "${experiment}_${chromosome}_pos.pileup" \\
            --outfile2 "${experiment}_${chromosome}_neg.pileup" \\
            --outfile3 "${experiment}_${chromosome}_vars.pileup"
        
        # Check if Python script ran successfully
        if [[ \$? -ne 0 ]]; then
            echo "ERROR: Python processing script failed"
            exit 1
        fi
    fi
    
    # Verify output files were created
    for file in "${experiment}_${chromosome}_pos.pileup" "${experiment}_${chromosome}_neg.pileup" "${experiment}_${chromosome}_vars.pileup"; do
        if [[ ! -f "\$file" ]]; then
            echo "ERROR: Expected output file not created: \$file"
            exit 1
        fi
    done
    
    # Compress variants file
    gzip "${experiment}_${chromosome}_vars.pileup"
    
    echo "Completed pileup processing for ${experiment}_${chromosome}"
    """
    
    stub:
    """
    touch "${experiment}_${chromosome}_pos.pileup"
    touch "${experiment}_${chromosome}_neg.pileup"
    touch "${experiment}_${chromosome}_vars.pileup.gz"
    """
}