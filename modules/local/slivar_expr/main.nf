process SLIVAR_EXPR {
    tag "$meta.id"
    label 'process_single'

    // No official biocontainer for slivar -- brentp (its author) publishes a Docker Hub image
    // directly, same convention as this repo's EXOMISER module for the same reason. Explicitly
    // qualified with docker.io/ so nextflow.config's global `docker.registry = 'quay.io'` doesn't
    // rewrite this into a quay.io reference that doesn't exist (confirmed directly -- an
    // unqualified 'brentp/slivar:v0.3.1' resolved to quay.io/brentp/slivar:v0.3.1 and failed to pull).
    container 'docker.io/brentp/slivar:v0.3.1'

    input:
    tuple val(meta), path(vcf), path(tbi), path(ped)
    path(js)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path("versions.yml")            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail

    slivar expr \\
        --vcf ${vcf} \\
        --ped ${ped} \\
        --js ${js} \\
        ${args} \\
        -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar 2>&1 | grep -o 'version: [^ ]*' | sed 's/version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slivar: \$(slivar 2>&1 | grep -o 'version: [^ ]*' | sed 's/version: //')
    END_VERSIONS
    """
}
