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
