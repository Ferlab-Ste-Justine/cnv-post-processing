/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { EXOMISER } from '../modules/local/exomiser/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EXOMISER_WORKFLOW {
    take:
    ch_input
    exomiser_genome
    exomiser_data_version
    exomiser_data_dir
    analysis_wes_path
    analysis_wgs_path
    local_frequency_file
    local_frequency_index_file
    remm_version
    remm_filename
    cadd_version
    cadd_snv_filename
    cadd_indel_filename

    main:
    def ch_versions = Channel.empty()
    def ch_input_for_exomiser = ch_input.map{
        meta, vcf ->
            def tbi = file(vcf + ".tbi")
            def pheno = file(meta["pheno"])
            def analysis_file = meta.sequencingType == "WES"? analysis_wes_path : analysis_wgs_path
            [meta, vcf, tbi.exists() ? tbi : [], pheno, analysis_file]
    }
    EXOMISER(ch_input_for_exomiser,
        exomiser_data_dir,
        exomiser_genome,
        exomiser_data_version,
        [local_frequency_file, local_frequency_index_file],
        remm_version? [remm_version, remm_filename] : ["", ""],
        cadd_version? [cadd_version, cadd_snv_filename, cadd_indel_filename] : ["", "", ""]
    )
    ch_versions = ch_versions.mix(EXOMISER.out.versions)

    emit:
    versions=ch_versions
}

workflow CNV_POST_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    def exomiser_local_frequency_file = params.exomiser_local_frequency_path? file(params.exomiser_local_frequency_path) : []
    def exomiser_local_frequency_index_file = params.exomiser_local_frequency_index_path? file(params.exomiser_local_frequency_index_path) : []
    def exomiser_data_dir = params.exomiser_data_dir? file(params.exomiser_data_dir) : []
    def exomiser_analysis_wes_path = params.exomiser_analysis_wes? file(params.exomiser_analysis_wes) : []
    def exomiser_analysis_wgs_path = params.exomiser_analysis_wgs? file(params.exomiser_analysis_wgs) : []

    def ch_versions = Channel.empty()

    EXOMISER_WORKFLOW(
        ch_samplesheet,
        params.exomiser_genome,
        params.exomiser_data_version,
        exomiser_data_dir,
        exomiser_analysis_wes_path,
        exomiser_analysis_wgs_path,
        exomiser_local_frequency_file,
        exomiser_local_frequency_index_file,
        params.exomiser_remm_version,
        params.exomiser_remm_filename,
        params.exomiser_cadd_version,
        params.exomiser_cadd_snv_filename,
        params.exomiser_cadd_indel_filename
    )
    ch_versions = ch_versions.mix(EXOMISER_WORKFLOW.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'cnv-post-processing_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions = ch_versions   // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
