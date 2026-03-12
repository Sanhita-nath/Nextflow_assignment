#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq.gz'
params.outdir = 'outputs/Sanhita/'
params.adapters = 'adapters.fa'
params.reference = 'LG12.fasta*'

log.info """
      LIST OF PARAMETERS
================================
Reads            : ${params.reads}
Output-folder    : ${params.outdir}
Adapters         : ${params.adapters}
"""

// Create read channel
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).map { sample, reads -> tuple(sample, reads.collect { it.toAbsolutePath() }) }
adapter_ch = Channel.fromPath(params.adapters)
reference_ch = Channel.fromPath(params.reference, checkIfExists: true).collect()

include { fastqc } from './modules/fastqc'
include { trimmomatic } from './modules/trimmomatic'
include { bwa_mem2 } from './modules/bwa_mem2'
 
//Added second argument for trimmomatic
// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trim_out = trimmomatic(read_pairs_ch, adapter_ch)
    bwa_mem2(trim_out.trimmed_fq, reference_ch)
}


