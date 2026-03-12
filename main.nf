#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = 'data/*_{1,2}.fq.gz'
params.outdir = 'outputs/Sanhita/'
params.adapters = 'adapters.fa'
params.reference = 'LG12.fasta'

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
reference_ch = Channel.fromPath(params.reference, checkIfExists: true)


// Define fastqc process
process fastqc {
    label "fastqc"
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    script:
    """
    fastqc ${reads}
    """
}

// Process trimmomatic
process trimmomatic {
    label "trimmomatic"
    publishDir "${params.outdir}/trimmed-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path adapters_file

    output:
    tuple val("${sample}"), path("${sample}*.trimmed.fq.gz"), emit: trimmed_fq
    tuple val("${sample}"), path("${sample}*.discarded.fq.gz"), emit: discarded_fq

    script:
    """
    trimmomatic PE -phred33 ${reads[0]} ${reads[1]} ${sample}_1.trimmed.fq.gz ${sample}_1.discarded.fq.gz ${sample}_2.trimmed.fq.gz ${sample}_2.discarded.fq.gz ILLUMINACLIP:${adapters_file}:2:30:10
    """
}

// define alignment_process
process bwa-mem2 {
    label "bwa-mem2"
    publishDir "${params.outdir}/aligned-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path reference

    output:
    tuple val(sample), path("${sample}.sam")
    
    script:
    """
    bwa-mem2 mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} > ${sample}.sam
    """
}
    
//Added second argument for trimmomatic
// Run the workflow
workflow {
    read_pairs_ch.view()
    fastqc(read_pairs_ch)
    trimmomatic(read_pairs_ch, adapter_ch)
}


