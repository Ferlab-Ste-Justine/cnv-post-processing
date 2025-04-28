process EXOMISER {
    tag "$meta.id"
    label 'process_medium'

    // Commenting this out because Conda is not supported at the moment.
    // conda "${moduleDir}/environment.yml"

    // We use a custom docker image for this module (should match exomiser 14.0.0)
    container 'registry.hub.docker.com/ferlabcrsj/exomiser:2.8.1'

    input:
    tuple val(meta), path(vcf_file), path(index_file), path(pheno_file), path(analysis_file)
    path data_dir
    val exomiser_genome
    val exomiser_data_version
    tuple path(local_frequency_path), path(local_frequency_index_path)
    tuple val(remm_version), val(remm_file_name)
    tuple val(cadd_version), val(cadd_snv_file_name),val(cadd_indel_file_name)

    output:
    tuple val(meta), path("*vcf.gz")         , optional:true, emit: vcf
    tuple val(meta), path("*vcf.gz.tbi")     , optional:true, emit: tbi
    tuple val(meta), path("*html")           , optional:true, emit: html
    tuple val(meta), path("*json")           , optional:true, emit: json
    tuple val(meta), path("*genes.tsv")      , optional:true, emit: genetsv
    tuple val(meta), path("*variants.tsv")   , optional:true, emit: variantstsv
    path("versions.yml")            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def application_properties_args = task.ext.application_properties_args ?: ''

    def local_frequency_file_args = ""
    if (local_frequency_path) {
        log.info("Using LOCAL frequency file {}", local_frequency_path)
        local_frequency_file_args = "--exomiser.${exomiser_genome}.local-frequency-path=/`pwd`/${local_frequency_path}"
    }

    def remm_args = ""
    if (remm_version) {
        log.info("Using REMM version {}", remm_version)
        remm_args += "--exomiser.remm.version=\"${remm_version}\""
        remm_args += " --exomiser.${exomiser_genome}.remm-path=/`pwd`/${data_dir}/remm/${remm_file_name}"
    }

    def cadd_args = ""
    if (cadd_version) {
        log.info("Using CADD version {}", cadd_version)
        cadd_args += "--cadd.version=\"${cadd_version}\""
        cadd_args += " --exomiser.${exomiser_genome}.cadd-snv-path=/`pwd`/${data_dir}/cadd/${cadd_version}/${cadd_snv_file_name}"
        cadd_args += " --exomiser.${exomiser_genome}.cadd-indel-path=/`pwd`/${data_dir}/cadd/${cadd_version}/${cadd_indel_file_name}"
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[EXOMISER] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    #!/bin/bash -eo pipefail

    java -Xmx${avail_mem}M -cp \$( cat /app/jib-classpath-file ) \$( cat /app/jib-main-class-file ) \\
        --vcf ${vcf_file} \\
        --assembly "${exomiser_genome}" \\
        --analysis "${analysis_file}" \\
        --sample ${pheno_file} \\
        --output-format=HTML,JSON,TSV_GENE,TSV_VARIANT,VCF \\
        --output-directory=/`pwd` \\
        ${args} \\
        --exomiser.data-directory=/`pwd`/${data_dir} \\
        ${local_frequency_file_args} \\
        ${remm_args} \\
        ${cadd_args} \\
        --exomiser.${exomiser_genome}.data-version="${exomiser_data_version}" \\
        --exomiser.phenotype.data-version="${exomiser_data_version}" \\
        ${application_properties_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "exomiser": "\$(cat /EXOMISER_VERSION.txt)"
    END_VERSIONS
    """


    stub:


    // Exomiser provides multiple ways to configure the output filename, which can
    // differ between Exomiser versions. Implementing a fully realistic stub to
    // account for these variations would be complex and significantly increase
    // the stub's implementation effort. To keep things simple, the stub uses a
    // fixed output filename prefix. The chosen prefix aligns with the default
    // behavior of the pipeline.
    def output_filename_prefix = meta.id + ".exomiser"

    """
    #!/bin/bash -eo pipefail
    touch ${output_filename_prefix}.genes.tsv
    touch ${output_filename_prefix}.html
    touch ${output_filename_prefix}.json
    touch ${output_filename_prefix}.variants.tsv
    touch ${output_filename_prefix}.vcf.gz
    touch ${output_filename_prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        "exomiser": "\$(cat /EXOMISER_VERSION.txt)"
    END_VERSIONS
    """
}
