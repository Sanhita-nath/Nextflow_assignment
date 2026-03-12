process bwa_mem2 {
    label "bwa_mem2"
    publishDir "${params.outdir}/aligned-reads-${sample}/", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    path reference_files

    output:
    tuple val(sample), path("${sample}.sam")

    script:
    """
    bwa-mem2 mem -t ${task.cpus} LG12.fasta ${reads[0]} ${reads[1]} > ${sample}.sam
    """
}
