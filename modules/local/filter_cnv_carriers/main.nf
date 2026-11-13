process FILTER_CNV_CARRIERS {
    tag "${meta.proband ? "${meta.id}.${meta.proband}" : meta.id}"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: vcf_tbi
    path("versions.yml")                                                      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // A proband VCF (JOIN_PROBAND_CNV) is named after its family and proband, like its slivar VCF
    def prefix = task.ext.prefix ?: (meta.proband ? "${meta.id}.${meta.proband}" : "${meta.id}")
    // Sets each record's FILTER from its carriers' FORMAT/FT: PASS when at least one carrier's FT
    // is PASS, the filters they name otherwise. ext.args --pass-only (--cnv_pass_only) also drops
    // the records that don't pass. See bin/filter_cnv_carriers.py.
    """
    #!/bin/bash -eo pipefail

    filter_cnv_carriers.py \\
        --vcf-in ${vcf} \\
        --vcf-out ${prefix}.filtered.vcf.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: (meta.proband ? "${meta.id}.${meta.proband}" : "${meta.id}")
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
