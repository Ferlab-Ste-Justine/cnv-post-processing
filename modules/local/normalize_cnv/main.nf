process NORMALIZE_CNV {
    tag "$meta.id"
    label 'process_low'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0':
        'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.normalized.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.normalized.vcf.gz.tbi") , emit: tbi
    path("versions.yml")                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // DRAGEN's own jointly-called VCFs declare INFO/SVLEN as Number=. rather than Number=A, so
    // `bcftools norm -m -any` can't cleanly re-split it per-allele on originally-multiallelic
    // sites -- some resulting records keep a comma-separated, multi-value SVLEN. htsjdk (used by
    // Exomiser) parses per the header's declared Number regardless of the actual value's
    // cardinality, so it reads this as a list even on single-value records and crashes trying to
    // cast it to a scalar (ClassCastException, confirmed directly against this exact file).
    // DRAGEN's own per-sample output was directly confirmed to have zero multiallelic sites this
    // same session, so this only needs to apply where task.ext.strip_svlen is set (the
    // DRAGEN_JOINT path). SVLEN isn't needed downstream anyway -- VEP/Exomiser derive size from
    // POS/END/SVTYPE.
    def strip_svlen = task.ext.strip_svlen ? '| bcftools annotate --no-version -x INFO/SVLEN' : ''

    // Only DRAGEN's CNV VCF conventions have been validated against this normalization step so
    // far (see /cephfs/jtrembla/projects/bioinfo-214-CNV-post-processing/preprocessing.sh). Fail
    // loudly for other callers rather than silently applying possibly-wrong assumptions.
    // DRAGEN_JOINT reuses this same module (see workflows/cnv_post_processing.nf) since the
    // filter/sort/split logic below is identical whether the input is per-sample or a family-level
    // jointly-called VCF.
    if (!(meta.caller in ['DRAGEN', 'DRAGEN_JOINT'])) {
        error("NORMALIZE_CNV: caller '${meta.caller}' is not yet supported (only DRAGEN and DRAGEN_JOINT are implemented).")
    }

    """
    #!/bin/bash -eo pipefail

    # Drop non-variant records (ALT=".", e.g. DRAGEN's REF segmentation blocks confirming normal
    # copy number) before anything else touches them -- they're not CNVs, so no downstream stage
    # (collapse, merge, depth refinement, VEP) has any use for them. Filtered on ALT rather than
    # DRAGEN's own ID naming convention so this keeps working for other callers.
    #
    # DRAGEN CNV VCFs can also carry multiallelic <DEL>,<DUP> records at a single site, which
    # breaks downstream truvari comparisons unless split first.
    # --no-version on every view/norm/annotate call below: bcftools appends a
    # "##bcftools_<cmd>Command=...; Date=<now>" header line by default, which made every
    # normalized VCF differ byte-for-byte between otherwise-identical runs (confirmed directly --
    # broke nf-test's pipeline-level snapshot for reasons unrelated to any real content change).
    # bcftools sort has no such flag because it never adds this header in the first place.
    bcftools view --no-version -e 'ALT="."' ${vcf} \\
        | bcftools sort -O u \\
        | bcftools norm --no-version -m -any -O u ${args} - \\
        ${strip_svlen} \\
        | bcftools view --no-version -O z -o ${prefix}.normalized.vcf.gz -

    tabix -p vcf ${prefix}.normalized.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash -eo pipefail
    echo "" | gzip > ${prefix}.normalized.vcf.gz
    touch ${prefix}.normalized.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //')
    END_VERSIONS
    """
}
