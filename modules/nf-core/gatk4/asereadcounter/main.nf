process GATK4_ASEREADCOUNTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b8c5d6e7f8a9b0c1d2e3f4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d0e1f2a3b4/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(vcf), path(vcf_tbi)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.ase.csv"), emit: ase_counts
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version 2>&1 | head -1 | sed 's/The Genome Analysis Toolkit (GATK) //'"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? "-Xmx${(task.memory.mega * 0.8).intValue()}m" : '-Xmx4g'
    """
    # Filter to passing biallelic SNPs
    gatk --java-options "${avail_mem} -XX:-UsePerfData" \\
        SelectVariants \\
        -R ${fasta} \\
        -V ${vcf} \\
        --select-type-to-include SNP \\
        --restrict-alleles-to BIALLELIC \\
        --exclude-filtered \\
        -O biallelic_snp.vcf.gz

    # Count allele-specific reads at each passing SNP position
    gatk --java-options "${avail_mem} -XX:-UsePerfData" \\
        ASEReadCounter \\
        -R ${fasta} \\
        -I ${bam} \\
        -V biallelic_snp.vcf.gz \\
        --tmp-dir . \\
        ${args} \\
        -O ${prefix}.ase.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ase.csv
    """
}
