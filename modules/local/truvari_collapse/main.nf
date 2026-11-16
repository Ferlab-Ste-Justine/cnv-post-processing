process TRUVARI_COLLAPSE {
    tag "$meta.id"
    label 'process_low'

    // Official biocontainer -- same truvari version already validated against the CNV-tuned
    // thresholds used here (see
    // /cephfs/jtrembla/projects/bioinfo-214-CNV-post-processing/preprocessing.sh). Doesn't bundle
    // bcftools/tabix, so the post-collapse sort+index is a separate step (BCFTOOLS_SORT_INDEX)
    // rather than folded into this module -- one tool per container.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:5.4.0--pyhdfd78af_0':
        'quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.collapsed.vcf.gz") , emit: vcf
    tuple val(meta), path("*.removed.vcf.gz")   , emit: removed
    path("versions.yml")                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/bin/bash -eo pipefail

    truvari collapse \\
        -i ${vcf} \\
        -o ${prefix}.collapsed.vcf.gz \\
        -c ${prefix}.removed.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(truvari version | sed 's/^Truvari v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.collapsed.vcf.gz
    echo "" | gzip > ${prefix}.removed.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(truvari version | sed 's/^Truvari v//')
    END_VERSIONS
    """
}
