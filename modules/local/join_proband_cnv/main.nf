process JOIN_PROBAND_CNV {
    tag "${meta.id}.${meta.proband}"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    tuple val(meta), path(family_vcf), path(family_tbi), path(proband_vcf), path(proband_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi
    path("versions.yml")                                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.proband}"
    // One record per call of the proband (meta.proband, its VCF sample name): the call itself,
    // from its VEP-annotated per-sample VCF (proband_vcf), with the FORMAT of the family record
    // (family_vcf) whose proband FORMAT/OID names it. slivar runs on the output. See
    // bin/join_proband_cnv.py.
    """
    #!/bin/bash -eo pipefail

    join_proband_cnv.py \\
        --family-vcf ${family_vcf} \\
        --proband-vcf ${proband_vcf} \\
        --proband '${meta.proband}' \\
        --vcf-out ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.proband}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
