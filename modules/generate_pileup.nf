process GENERATE_PILEUP {
    label 'standard'
    tag "${experiment}_${chromosome}"
    publishDir "${params.outdir}/pileup_files", mode: 'copy', pattern: "*_norm.pileupt"
    
    input:
    tuple val(experiment), val(chromosome), path(bam_file), path(bai_file)
    
    output:
    tuple val(experiment), val(chromosome), path("${experiment}_${chromosome}_norm.pileupt"), emit: pileup_files
    
    script:
    """
    set -euo pipefail
    echo "Generating pileup for ${experiment}_${chromosome}"
    
    # Check if reference file exists
    if [[ ! -f "${params.refpath}/${chromosome}.fa" ]]; then
        echo "ERROR: Reference file not found: ${params.refpath}/${chromosome}.fa"
        exit 1
    fi
    
    # Get raw pileup file from bcftools
    echo "Running bcftools mpileup..."
    bcftools mpileup \\
        -d ${params.max_depth} \\
        -q ${params.mapping_quality} \\
        -Q ${params.base_quality} \\
        -a FORMAT/AD,INFO/AD \\
        -Ov \\
        -f ${params.refpath}/${chromosome}.fa \\
        "${bam_file}" > "${experiment}_${chromosome}.vcf"
    
    # Check if VCF was created successfully
    if [[ ! -s "${experiment}_${chromosome}.vcf" ]]; then
        echo "ERROR: Failed to generate VCF file"
        exit 1
    fi
    
    # Normalize the VCF file
    echo "Normalizing VCF..."
    bcftools norm \\
        -f "${params.refpath}/${chromosome}.fa" \\
        "${experiment}_${chromosome}.vcf" > "${experiment}_${chromosome}_norm.vcf"
    
    # Clean up intermediate file
    rm -f "${experiment}_${chromosome}.vcf"
    
    # Convert to pileup format
    echo "Converting to pileup format..."
    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%TYPE\\t%ALT{0}\\t%IMF\\t%DP\\t%AD\\n' \\
        "${experiment}_${chromosome}_norm.vcf" | \\
        awk 'BEGIN{OFS="\\t"} {print \$1"_"\$2,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' | \\
        sed 's/<\\*>/R/g' > "${experiment}_${chromosome}_norm.pileupt"
    
    # Clean up VCF file
    rm -f "${experiment}_${chromosome}_norm.vcf"
    
    # Check if final output was created
    if [[ ! -s ${experiment}_${chromosome}_norm.pileupt ]]; then
        echo "WARNING: Empty pileup file generated for ${experiment}_${chromosome}"
        # Create minimal pileup file to avoid pipeline failure
        echo -e "ID\\tCHROM\\tPOS\\tREF\\tTYPE\\tALT\\tIMF\\tDP\\tAD" > "${experiment}_${chromosome}_norm.pileupt"
    fi
    
    echo "Completed pileup generation for ${experiment}_${chromosome}"
    """
    
    stub:
    """
    touch ${experiment}_${chromosome}_norm.pileupt
    """
}
