include { EXTRACT_CHROMOSOME_BAM } from '../modules/extract_chromosome_bam'
include { GENERATE_PILEUP } from '../modules/generate_pileup'
include { PROCESS_PILEUP } from '../modules/process_pileup'
include { FILTER_COVERAGE } from '../modules/filter_coverage'
include { NORMALIZE_RATES } from '../modules/normalize_rates'
include { COMBINE_CHROMOSOMES } from '../modules/combine_chromosomes'
include { FIND_COMMON_SITES } from '../modules/find_common_sites'
include { SET_COMMON_SITES } from '../modules/set_common_sites'
include { COMPUTE_ELEMENT_STATS } from '../modules/compute_element_stats'
include { COMPUTE_ELEMENT_STATS_SPLICED } from '../modules/compute_element_stats_spliced'

workflow TRAM_SEQ {
    
    // Create channels for all experiments
    all_experiments = Channel.from(
        params.control_experiments + 
        params.dms_experiments + 
        params.stress_control_experiments + 
        params.stress_dms_experiments
    )
    
    // Create channel for chromosomes
    chromosomes_ch = Channel.from(params.chromosomes)
    
    // Create experiment-chromosome combinations
    experiment_chromosome_ch = all_experiments
        .combine(chromosomes_ch)
    
    // Extract chromosome-specific BAM files
    EXTRACT_CHROMOSOME_BAM(experiment_chromosome_ch)
    
    // Generate pileup files
    GENERATE_PILEUP(EXTRACT_CHROMOSOME_BAM.out.bam_files)

    // Process pileup files
    PROCESS_PILEUP(GENERATE_PILEUP.out.pileup_files)

    // Combine pos and neg pileup files for each experiment-chromosome combination
    pileup_combined_ch = PROCESS_PILEUP.out.pos_pileup
        .join(PROCESS_PILEUP.out.neg_pileup, by: [0, 1])
        .map { experiment, chromosome, pos_file, neg_file ->
            [experiment, chromosome, pos_file, neg_file]
        }

    // Filter out sites with coverage that is below threshold
    FILTER_COVERAGE(pileup_combined_ch)

    // Normalize rates from DMS-treated samples using control samples
    dms_experiments_ch = Channel.from(params.dms_experiments + params.stress_dms_experiments)
    chromosomes_ch = Channel.from(params.chromosomes)
    strands_ch = Channel.from(['pos', 'neg'])

    // Create experiment pairs channel and combine with chromosome and strand
    normalization_input_ch = dms_experiments_ch
    	.map { dms_exp -> 
             def control_exp = params.experiment_pairs[dms_exp]
             if (control_exp == null) {
             	error "No control experiment found for DMS experiment: ${dms_exp}. Check your experiment_pairs in parameters."
             }
             [dms_exp, control_exp]
    	}
    	.combine(chromosomes_ch)
    	.combine(strands_ch)

    // Wait for FILTER_COVERAGE to complete
   filter_coverage_complete = FILTER_COVERAGE.out.coverage_pos_bedgraph
    	.mix(FILTER_COVERAGE.out.coverage_neg_bedgraph)
    	.mix(FILTER_COVERAGE.out.mmrate_pos_bedgraph)
    	.mix(FILTER_COVERAGE.out.mmrate_neg_bedgraph)
    	.collect()
    	.map { "ready" }

   // Combine normalization input with completion signal
   filtered = normalization_input_ch
    	.combine(filter_coverage_complete)
    	.map { tuple ->
             def dms_exp = tuple[0]
             def control_exp = tuple[1] 
             def chromosome = tuple[2]
             def strand = tuple[3]
             def ready_signal = tuple[4]
             [dms_exp, control_exp, chromosome, strand]
    	}
    
    // Wait for FILTER_COVERAGE to complete before starting normalization
    NORMALIZE_RATES(filtered)        

    // Combine chromosomes into genome-wide bedGraph files
    dms_experiments_only = Channel.from(params.dms_experiments + params.stress_dms_experiments)
    assays_ch = Channel.from(['coverage', 'mmRate'])
    
    combination_input_ch = dms_experiments_only
        .combine(assays_ch)
        .combine(strands_ch)
    
    // Gather all NORMALIZE_RATES outputs
    normalize_rates_complete = NORMALIZE_RATES.out.normalized_coverage
        .mix(NORMALIZE_RATES.out.normalized_mmrate)
        .collect()
        .map { "ready" }
    
    // Combine combination input with completion signal
    normalized = combination_input_ch
        .combine(normalize_rates_complete)
        .map { experiment, assay, strand, ready_signal ->
            [experiment, assay, strand]
        }
    
    COMBINE_CHROMOSOMES(normalized)

    // Find common sites to retain across all DMS experiments
    common_sites_input_ch = Channel.from(['pos', 'neg'])
    
    // Wait for all COMBINE_CHROMOSOMES coverage outputs
    combine_chromosomes_complete = COMBINE_CHROMOSOMES.out.genome_wide_bedgraph
        .filter { experiment, assay, strand, file -> assay == 'coverage' }
        .collect()
        .map { "ready" }
    
    // Combine common sites input with completion signal
    combined = common_sites_input_ch
        .combine(combine_chromosomes_complete)
        .map { strand, ready_signal ->
            strand
        }
    
    FIND_COMMON_SITES(combined)

    // Retain common sites for all DMS experiments
    set_common_input_ch = dms_experiments_only
        .combine(assays_ch)
        .combine(strands_ch)

    // Wait for all FIND_COMMON_SITES outputs
    find_common_sites_complete = FIND_COMMON_SITES.out.common_sites
        .collect()
        .map { "ready" }

    // Combine set common input with completion signal
    set_combined = set_common_input_ch
        .combine(find_common_sites_complete)
        .map { experiment, assay, strand, ready_signal ->
            [experiment, assay, strand]
        }
    
    SET_COMMON_SITES(set_combined)

    // Wait for all SET_COMMON_SITES outputs
    set_common_sites_complete = SET_COMMON_SITES.out.common_filtered_bedgraph
        .collect()
        .map { "ready" }

    if (params.references) {
    // Compute element statistics
        reference_ch = Channel.from(params.references)
        
    // Combine reference input with completion signal
        elements = reference_ch
            .combine(set_common_sites_complete)
            .map { reference, ready_signal ->
                reference
            }
    
        COMPUTE_ELEMENT_STATS(elements)
    }

    if (params.references_spliced) {
    // Compute element statistics
        spliced_reference_ch = Channel.from(params.references_spliced)
        
    // Combine reference input with completion signal
        elements = spliced_reference_ch
            .combine(set_common_sites_complete)
            .map { reference, ready_signal ->
                reference
            }
    
        COMPUTE_ELEMENT_STATS_SPLICED(elements)
    }

    emit:
    pos_pileup = PROCESS_PILEUP.out.pos_pileup
    neg_pileup = PROCESS_PILEUP.out.neg_pileup
    vars_pileup = PROCESS_PILEUP.out.vars_pileup
    
    coverage_pos_bedgraph = FILTER_COVERAGE.out.coverage_pos_bedgraph
    coverage_neg_bedgraph = FILTER_COVERAGE.out.coverage_neg_bedgraph
    mmrate_pos_bedgraph = FILTER_COVERAGE.out.mmrate_pos_bedgraph
    mmrate_neg_bedgraph = FILTER_COVERAGE.out.mmrate_neg_bedgraph
    compressed_pos_pileup = FILTER_COVERAGE.out.compressed_pos_pileup
    compressed_neg_pileup = FILTER_COVERAGE.out.compressed_neg_pileup

    normalized_coverage = NORMALIZE_RATES.out.normalized_coverage
    normalized_mmrate = NORMALIZE_RATES.out.normalized_mmrate

    genome_wide_bedgraphs = COMBINE_CHROMOSOMES.out.genome_wide_bedgraph

    common_sites = FIND_COMMON_SITES.out.common_sites

    common_filtered_bedgraphs = SET_COMMON_SITES.out.common_filtered_bedgraph

    coverage_stats = COMPUTE_ELEMENT_STATS.out.coverage_stats
    mmrate_stats = COMPUTE_ELEMENT_STATS.out.mmrate_stats
    processed_stats = COMPUTE_ELEMENT_STATS.out.processed_stats
    spliced_coverage_stats = COMPUTE_ELEMENT_STATS_SPLICED.out.coverage_stats
    spliced_mmrate_stats = COMPUTE_ELEMENT_STATS_SPLICED.out.mmrate_stats
    spliced_processed_stats = COMPUTE_ELEMENT_STATS_SPLICED.out.processed_stats
}
