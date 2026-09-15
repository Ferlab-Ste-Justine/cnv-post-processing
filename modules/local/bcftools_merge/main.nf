process BCFTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0':
        'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.merged.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.merged.vcf.gz.tbi") , emit: tbi
    path("versions.yml")                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/bin/bash -eo pipefail

    bcftools merge -O u ${args} ${vcfs} \\
        | bcftools sort -O u - \\
        | bcftools norm -m -any -O z -o ${prefix}.merged.vcf.gz -

    tabix -p vcf ${prefix}.merged.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
