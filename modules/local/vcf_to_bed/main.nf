process VCF_TO_BED {
    tag "$meta.id"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0':
        'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.bed") , emit: bed
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail

    bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%ID\\n' ${vcf} \\
        | awk 'BEGIN{OFS="\\t"} {print \$1, \$2-1, \$3, \$4}' > ${prefix}.sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    touch ${prefix}.sites.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
