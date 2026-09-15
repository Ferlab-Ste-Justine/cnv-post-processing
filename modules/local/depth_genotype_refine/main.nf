process DEPTH_GENOTYPE_REFINE {
    tag "$meta.id"
    label 'process_single'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1':
        'quay.io/biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(ratio_beds)

    output:
    tuple val(meta), path("*.refined.vcf.gz") , emit: vcf
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Resolves missing (./.) genotypes to hom-ref (0/0) only where a sample's mosdepth-derived
    // depth ratio at that site actually supports "no CNV here" -- schema.txt's "missing ->
    // hom-ref vs missing -> no-data" goal, not new CNV calling. This 0.8-1.2 band is a reasoned
    // default (symmetric around 1.0), NOT a validated clinical threshold -- needs tuning against
    // real data (see the NA12878/91/92 truvari-bench ground truth) before being trusted
    // clinically.
    def low = task.ext.low ?: '0.8'
    def high = task.ext.high ?: '1.2'

    """
    #!/bin/bash -eo pipefail

    refine_genotypes.py \\
        --vcf-in ${vcf} \\
        --vcf-out ${prefix}.refined.vcf.gz \\
        --low ${low} --high ${high} \\
        ${ratio_beds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.refined.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
