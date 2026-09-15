process BCFTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_low'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0':
        'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.sorted.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.sorted.vcf.gz.tbi") , emit: tbi
    path("versions.yml")                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/bin/bash -eo pipefail

    # truvari collapse's own output isn't guaranteed coordinate-sorted (observed directly:
    # tabix fails on it with "Chromosome blocks not continuous" otherwise) -- this step exists
    # specifically to re-sort before indexing rather than assuming a VCF is safe to index as-is.
    bcftools sort -O z -o ${prefix}.sorted.vcf.gz ${args} ${vcf}
    tabix -p vcf ${prefix}.sorted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.sorted.vcf.gz
    touch ${prefix}.sorted.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
