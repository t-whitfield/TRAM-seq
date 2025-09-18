#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRAM_SEQ } from './workflows/tram_seq'

workflow {
    TRAM_SEQ()
}
