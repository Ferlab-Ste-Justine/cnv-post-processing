process FILL_CN_FROM_SEGMENTS {
    tag "$meta.id"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    // segment_vcfs: each family member's original per-sample VCF, REF blocks included, matched to
    // the family VCF's columns by their own sample name. Staged one per folder: two samplesheet
    // rows can point at files with the same name. ped: [] when the family has none.
    tuple val(meta), path(vcf), path(tbi), path(segment_vcfs, stageAs: 'segments??/*'), path(ped)

    output:
    tuple val(meta), path("*.cn_filled.vcf.gz") , emit: vcf
    path("versions.yml")                        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Records under 5 kb are left missing: the per-sample caller misses many small losses the
    // joint caller finds, and a parent filled as normal there turns the child's inherited CNV into
    // `denovo` (see docs/journal.md).
    def min_size = task.ext.min_size != null ? task.ext.min_size : 5000
    // The share of a record one CN must cover in the sample's own segments (seg_c50, 2026-11-10).
    def min_cover = task.ext.min_cover != null ? task.ext.min_cover : 0.5
    // The PED's sex column gives each sample's expected CN on non-PAR chrX/chrY: without one, a
    // CN found there is still written, with GT left missing.
    def ped_arg = ped ? "--ped ${ped}" : ''

    """
    #!/bin/bash -eo pipefail

    fill_cn_from_segments.py \\
        --vcf-in ${vcf} \\
        --vcf-out ${prefix}.cn_filled.vcf.gz \\
        --min-size ${min_size} \\
        --min-cover ${min_cover} \\
        ${ped_arg} \\
        ${segment_vcfs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.cn_filled.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
