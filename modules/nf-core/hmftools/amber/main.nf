process HMFTOOLS_AMBER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9e16ad88e25f7eb04cbf7fce03fe3f8d4e4c1b6e7349148d63c3cf0dac4e7c14/data'
        : 'community.wave.seqera.io/library/hmftools-amber:4.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fai
    path het_pon
    path het_pon_tbi
    val  ref_genome_version

    output:
    tuple val(meta), path("amber_out"), emit: amber_dir
    tuple val("${task.process}"), val('hmftools-amber'), eval("amber -version 2>&1 | head -1 | sed 's/.*AMBER version: //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def tumor_sm    = meta.tumor_sample  ?: meta.id
    def normal_sm   = meta.normal_sample ?: meta.id
    def het_pon_arg = het_pon            ? "-loci ${het_pon}" : ''
    """
    mkdir -p amber_out

    amber \\
        -tumor            ${tumor_sm} \\
        -tumor_bam        ${tumor_bam} \\
        -reference        ${normal_sm} \\
        -reference_bam    ${normal_bam} \\
        -output_dir       amber_out \\
        ${het_pon_arg} \\
        -ref_genome       ${fasta} \\
        -ref_genome_version ${ref_genome_version} \\
        -threads          ${task.cpus} \\
        ${args}
    """

    stub:
    def tumor_sm = meta.tumor_sample ?: meta.id
    """
    mkdir -p amber_out
    touch amber_out/${tumor_sm}.amber.baf.tsv.gz
    touch amber_out/${tumor_sm}.amber.baf.tsv.gz.tbi
    touch amber_out/${tumor_sm}.amber.baf.pcf
    touch amber_out/${tumor_sm}.amber.qc
    """
}
