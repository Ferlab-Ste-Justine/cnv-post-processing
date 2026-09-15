process MOSDEPTH_RATIO {
    tag "$meta.id"
    label 'process_medium'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.11--h0ec343a_1':
        'quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_1' }"

    input:
    tuple val(meta), path(cram), path(crai), path(sites_bed) // sites_bed varies per family, unlike fasta below
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.ratio.bed") , emit: ratio_bed
    path("versions.yml")                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-x -n'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail

    mosdepth -f ${fasta} --by ${sites_bed} -t ${task.cpus} ${args} ${prefix} ${cram}

    # Ratio of each region's mean depth to this sample's genome-wide mean depth (the "total" row
    # in mosdepth's own summary output) -- a simpler normalization than duphold's GC-matched bins,
    # but transparent and not dependent on an unmaintained tool. See DEPTH_GENOTYPE_REFINE for how
    # this is used.
    genome_mean=\$(awk -F'\\t' '\$1=="total"{print \$4}' ${prefix}.mosdepth.summary.txt)

    # mosdepth's --by preserves the input BED's 4th column (here, the DRAGEN record ID) as the
    # region name and appends its own computed mean depth as a 5th column -- so the depth to use
    # is \$5, not \$4.
    zcat ${prefix}.regions.bed.gz \\
        | awk -F'\\t' -v mean="\${genome_mean}" 'BEGIN{OFS="\\t"} { ratio = (\$5>0 && mean>0) ? \$5/mean : "NA"; print \$1, \$2, \$3, ratio }' \\
        > ${prefix}.ratio.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^mosdepth //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    touch ${prefix}.ratio.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^mosdepth //')
    END_VERSIONS
    """
}
