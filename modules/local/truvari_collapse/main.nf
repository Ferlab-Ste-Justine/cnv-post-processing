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
    // Proband priority (TRUVARI_COLLAPSE_COHORT, see conf/modules.config): when the family has
    // probands (meta.probands, from its PED), keep a proband's CNV as the representative of every
    // cluster it's in, instead of the leftmost call. truvari can't prefer a sample, so
    // bin/proband_priority.py rewrites QUAL (proband records 1, others 0, real QUAL saved in
    // INFO) for --keep maxqual, and puts the real QUAL back afterwards. Ties fall back to
    // truvari's leftmost rule, so clusters without a proband call are unchanged.
    def prioritize = task.ext.proband_priority && meta.probands
    def input_vcf = prioritize ? "${prefix}.proband_marked.vcf.gz" : "${vcf}"
    def keep_args = prioritize ? '--keep maxqual' : ''
    def mark = prioritize
        ? "proband_priority.py mark --vcf-in ${vcf} --vcf-out ${input_vcf} ${meta.probands.collect { p -> "--proband '${p}'" }.join(' ')}"
        : ''
    def restore = prioritize
        ? "proband_priority.py restore ${prefix}.collapsed.vcf.gz ${prefix}.removed.vcf.gz"
        : ''

    """
    #!/bin/bash -eo pipefail

    ${mark}

    truvari collapse \\
        -i ${input_vcf} \\
        -o ${prefix}.collapsed.vcf.gz \\
        -c ${prefix}.removed.vcf.gz \\
        ${args} ${keep_args}

    ${restore}

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
