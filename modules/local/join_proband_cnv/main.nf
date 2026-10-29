process JOIN_PROBAND_CNV {
    tag "${meta.id}.${meta.proband}"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    tuple val(meta), path(slivar_vcf), path(slivar_tbi), path(proband_vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf_tbi
    path("versions.yml")                                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.proband}"
    // One record per slivar record where the proband (meta.proband, its VCF sample name) has a
    // call: the proband's own call from its caller VCF (proband_vcf), joined through its
    // FORMAT/OID, plus slivar's tags and the family's FORMAT. See bin/join_proband_cnv.py.
    """
    #!/bin/bash -eo pipefail

    join_proband_cnv.py \\
        --slivar-vcf ${slivar_vcf} \\
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
