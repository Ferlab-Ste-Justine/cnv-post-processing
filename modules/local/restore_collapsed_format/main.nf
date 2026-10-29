process RESTORE_COLLAPSED_FORMAT {
    tag "$meta.id"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    tuple val(meta), path(collapsed_vcf), path(removed_vcf)

    output:
    tuple val(meta), path("*.restored.vcf.gz") , emit: vcf
    path("versions.yml")                       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // truvari collapse only copies an absorbed call's FORMAT into the kept record when it reads
    // that call's GT as set, which DRAGEN's `./1` (gain of unknown zygosity) and haploid `1`
    // aren't. This copies them back from truvari's removed VCF (MatchId -> CollapseId), so no
    // sample loses its call. Output unsorted, like truvari's: BCFTOOLS_SORT_INDEX follows.

    """
    #!/bin/bash -eo pipefail

    restore_collapsed_format.py \\
        --collapsed-vcf ${collapsed_vcf} \\
        --removed-vcf ${removed_vcf} \\
        --vcf-out ${prefix}.restored.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.restored.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
