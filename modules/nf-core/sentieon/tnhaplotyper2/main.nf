process SENTIEON_TNHAPLOTYPER2 {
    tag "${meta.id}"
    label 'process_high'

    // Sentieon is not available via conda — container only.
    container "community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5"

    input:
    tuple val(meta), path(tumor_cram), path(tumor_crai), path(normal_cram), path(normal_crai)
    path fasta
    path fai
    path germline_vcf
    path germline_vcf_tbi
    path contamination_vcf
    path contamination_vcf_tbi

    output:
    tuple val(meta), path("${prefix}_tnhap2-tmp.vcf.gz"),          emit: vcf_tmp
    tuple val(meta), path("${prefix}_tnhap2-tmp.vcf.gz.tbi"),      emit: vcf_tmp_tbi
    tuple val(meta), path("${prefix}_orientation"),                 emit: orientation
    tuple val(meta), path("${prefix}_contamination"),               emit: contamination
    tuple val(meta), path("${prefix}_contamination-segments"),      emit: contamination_segments
    tuple val("${task.process}"), val('sentieon'), eval("sentieon driver --version 2>&1 | head -1 | sed 's/sentieon-genomics-//'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    prefix           = task.ext.prefix ?: "${meta.id}"
    def tumor_sm     = meta.tumor_sample  ?: meta.id
    def normal_sm    = meta.normal_sample ?: meta.id
    def germline_arg = germline_vcf ? "--germline_vcf ${germline_vcf}" : ''
    """
    sentieon driver \\
        -r ${fasta} \\
        -t ${task.cpus} \\
        -i ${tumor_cram} \\
        -i ${normal_cram} \\
        --algo TNhaplotyper2 \\
            --tumor_sample  ${tumor_sm} \\
            --normal_sample ${normal_sm} \\
            ${germline_arg} \\
            ${prefix}_tnhap2-tmp.vcf.gz \\
        --algo OrientationBias \\
            --tumor_sample  ${tumor_sm} \\
            ${prefix}_orientation \\
        --algo ContaminationModel \\
            --tumor_sample  ${tumor_sm} \\
            --normal_sample ${normal_sm} \\
            --vcf           ${contamination_vcf} \\
            --tumor_segments ${prefix}_contamination-segments \\
            ${prefix}_contamination \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_tnhap2-tmp.vcf.gz
    touch ${prefix}_tnhap2-tmp.vcf.gz.tbi
    touch ${prefix}_orientation
    touch ${prefix}_contamination
    touch ${prefix}_contamination-segments
    """
}
