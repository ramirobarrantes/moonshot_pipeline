process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/03/0330fc16eba3c55ee8cb3e3c5b2fe8d64b6a49e77e50c20b6a0d6c1d2e2e3a7a/data'
        : 'docker.io/ensemblorg/ensembl-vep:release_114.0'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  cache_dir
    val   genome
    val   cache_version

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.summary.html"), emit: report, optional: true
    tuple val("${task.process}"), val('ensembl-vep'), eval("vep --help 2>&1 | grep 'Versions:' -A1 | tail -1 | sed 's/[^0-9.]//g'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vep \\
        --input_file   ${vcf} \\
        --output_file  ${prefix}.vep.vcf \\
        --vcf \\
        --format       vcf \\
        --assembly     ${genome} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache    ${cache_dir} \\
        --offline \\
        --fork         ${task.cpus} \\
        ${args}

    bgzip --threads ${task.cpus} ${prefix}.vep.vcf
    tabix -p vcf ${prefix}.vep.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vep.vcf.gz
    touch ${prefix}.vep.vcf.gz.tbi
    """
}
