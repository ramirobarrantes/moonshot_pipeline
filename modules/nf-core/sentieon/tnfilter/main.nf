process SENTIEON_TNFILTER {
    tag "${meta.id}"
    label 'process_medium'

    // Sentieon is not available via conda — container only.
    container "community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(orientation), path(contamination), path(contamination_segments)
    path fasta
    path fai
    val  sentieon_license

    output:
    tuple val(meta), path("${prefix}_tnhap2.vcf.gz"),     emit: vcf
    tuple val(meta), path("${prefix}_tnhap2.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('sentieon'), eval("sentieon driver --version 2>&1 | head -1 | sed 's/sentieon-genomics-//'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def tumor_sm  = meta.tumor_sample  ?: meta.id
    def normal_sm = meta.normal_sample ?: meta.id
    """
    export SENTIEON_LICENSE=${sentieon_license}

    sentieon driver \\
        -r ${fasta} \\
        --algo TNfilter \\
            -v              ${vcf} \\
            --tumor_sample  ${tumor_sm} \\
            --normal_sample ${normal_sm} \\
            --contamination ${contamination} \\
            --tumor_segments ${contamination_segments} \\
            --orientation_priors ${orientation} \\
            ${args} \\
            ${prefix}_tnhap2.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_tnhap2.vcf.gz
    touch ${prefix}_tnhap2.vcf.gz.tbi
    """
}
