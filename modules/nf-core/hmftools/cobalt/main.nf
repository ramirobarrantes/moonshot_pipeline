process HMFTOOLS_COBALT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hmftools-cobalt:2.2--hdfd78af_0'
        : 'quay.io/biocontainers/hmftools-cobalt:2.2--hdfd78af_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fai
    path gc_profile
    val ref_genome_version

    output:
    tuple val(meta), path("cobalt_out"), emit: cobalt_dir
    tuple val("${task.process}"), val('hmftools-cobalt'), eval("cobalt -version 2>&1 | head -1 | sed 's/.*COBALT version: //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def tumor_sm  = meta.tumor_sample  ?: meta.id
    def normal_sm = meta.normal_sample ?: meta.id
    """
    mkdir -p cobalt_out

    export _JAVA_OPTIONS="-Xmx220g"
    cobalt \\
        -tumor            ${tumor_sm} \\
        -tumor_bam        ${tumor_bam} \\
        -reference        ${normal_sm} \\
        -reference_bam    ${normal_bam} \\
        -output_dir       cobalt_out \\
        -gc_profile           ${gc_profile} \\
        -ref_genome           ${fasta} \\
        -ref_genome_version   ${ref_genome_version} \\
        -threads              ${task.cpus} \\
        ${args}
    """

    stub:
    def tumor_sm = meta.tumor_sample ?: meta.id
    """
    mkdir -p cobalt_out
    touch cobalt_out/${tumor_sm}.cobalt.ratio.tsv.gz
    touch cobalt_out/${tumor_sm}.cobalt.ratio.tsv.gz.tbi
    touch cobalt_out/${tumor_sm}.cobalt.ratio.pcf
    touch cobalt_out/${tumor_sm}.cobalt.gc.median.tsv
    """
}
