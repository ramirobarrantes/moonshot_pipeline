process HMFTOOLS_PURPLE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4f/4f2e9c5c8d1a8e6e7a9b3c0d1f2e3a4b5c6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f/data'
        : 'community.wave.seqera.io/library/hmftools-purple:4.0.2--hdfd78af_0'}"

    input:
    tuple val(meta), path(amber_dir), path(cobalt_dir), path(somatic_vcf), path(somatic_vcf_tbi)
    path fasta
    path fai
    path gc_profile
    val  ref_genome_version
    path ensembl_data_dir

    output:
    tuple val(meta), path("${tumor_sm}.purple.purity.tsv"),          emit: purity
    tuple val(meta), path("${tumor_sm}.purple.purity.range.tsv"),    emit: purity_range
    tuple val(meta), path("${tumor_sm}.purple.qc"),                  emit: qc
    tuple val(meta), path("${tumor_sm}.purple.cnv.somatic.tsv"),     emit: cnv_somatic
    tuple val(meta), path("${tumor_sm}.purple.cnv.gene.tsv"),        emit: cnv_gene
    tuple val(meta), path("${tumor_sm}.purple.segment.tsv"),         emit: segments
    tuple val(meta), path("${tumor_sm}.purple.somatic.vcf.gz"),      emit: somatic_vcf,    optional: true
    tuple val(meta), path("${tumor_sm}.purple.somatic.vcf.gz.tbi"),  emit: somatic_vcf_tbi, optional: true
    tuple val(meta), path("plot"),                                   emit: plots,          optional: true
    tuple val("${task.process}"), val('hmftools-purple'), eval("purple -version 2>&1 | head -1 | sed 's/.*PURPLE version: //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    tumor_sm         = meta.tumor_sample  ?: meta.id
    def normal_sm    = meta.normal_sample ?: meta.id
    def somatic_arg  = somatic_vcf        ? "-somatic_vcf ${somatic_vcf}" : ''
    """
    purple \\
        -tumor               ${tumor_sm} \\
        -reference           ${normal_sm} \\
        -amber               ${amber_dir} \\
        -cobalt              ${cobalt_dir} \\
        -gc_profile          ${gc_profile} \\
        -ref_genome          ${fasta} \\
        -ref_genome_version  ${ref_genome_version} \\
        -ensembl_data_dir    ${ensembl_data_dir} \\
        ${somatic_arg} \\
        -output_dir          . \\
        -threads             ${task.cpus} \\
        ${args}
    """

    stub:
    tumor_sm = meta.tumor_sample ?: meta.id
    """
    touch ${tumor_sm}.purple.purity.tsv
    touch ${tumor_sm}.purple.purity.range.tsv
    touch ${tumor_sm}.purple.qc
    touch ${tumor_sm}.purple.cnv.somatic.tsv
    touch ${tumor_sm}.purple.cnv.gene.tsv
    touch ${tumor_sm}.purple.segment.tsv
    mkdir -p plot
    """
}
