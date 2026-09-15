//
// Subworkflow with functionality specific to the Ferlab-Ste-Justine/cnv-post-processing pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    // Samples are grouped by familyId so that downstream cohort-level steps (bcftools merge,
    // truvari collapse, mosdepth, mendelian2, slivar) know each family's size upfront
    // (needed for groupKey-based grouping) and so that a family's familyPed/familyPheno are
    // validated identical across all of its rows before the main workflow body runs.

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map { meta, vcf, cram ->
            [meta.familyId, [meta, vcf, cram]]
        }
        .tap { ch_sample_simple }
        .groupTuple()
        .map { familyId, items ->
            validatePedFiles(familyId, items)
            validatePhenopacketFiles(familyId, items)
            validateJointCallerExclusivity(familyId, items)
            [familyId, items.size()]
        }
        .combine(ch_sample_simple, by: 0)
        .map { _familyId, sampleSize, item ->
            def meta = item[0]
            def new_meta = meta + [sampleSize: sampleSize, id: "${meta.familyId}.${meta.sample}".toString()]
            [new_meta, item[1], item[2]]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    // These attributes won't be available in the completion workflow because the workflow
    // object will be null. So we need to extract them from the workflow object beforehand.
    def Map summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    def metadata_dict = [:] + summary_params
    metadata_dict["Core Nextflow options"]["Start date"] =  workflow.start.toString()
    metadata_dict["Core Nextflow options"]["Command Line"] = workflow.commandLine
    metadata_dict["Core Nextflow options"]["CommitId"] = workflow.commitId

    def config_files = workflow.configFiles
    def command_line = workflow.commandLine

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)

        addMetadataFile(metadata_dict, "${outdir}/pipeline_info/metadata.txt")
        copyConfigurationFiles(config_files, "${outdir}/pipeline_info/configs/")

        def nextflow_log_file_path = getNextflowLogFilePath(command_line)
        copyNextflowLogFile(nextflow_log_file_path,  "${outdir}/pipeline_info/nextflow.log")
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS FROM NF-CORE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (meta, _vcf, _cram) = input
    // caller='DRAGEN_JOINT' rows skip the per-sample normalize/collapse chain entirely (see
    // workflows/cnv_post_processing.nf) -- exomiser (single-sample mode) reads its input straight
    // from that chain's output, so a DRAGEN_JOINT row's 'pheno' would never actually be used by
    // anything. Fail loudly here rather than silently ignoring it, which would otherwise look
    // like a missing feature rather than a samplesheet mistake.
    if (meta.caller == 'DRAGEN_JOINT' && meta.pheno) {
        error("Sample '${meta.sample}' (family '${meta.familyId}'): caller='DRAGEN_JOINT' rows don't go through the per-sample chain exomiser (single-sample mode) reads from, so 'pheno' would silently never be used. Remove 'pheno' for this row (use 'familyPheno' for exomiser family mode instead).")
    }
    return input
}

// All rows of the same family must reference the same familyPed value (or all omit it).
// familyPed is optional overall (families without one skip the mendelian2/slivar steps),
// but it must not be ambiguous within a family.
def validatePedFiles(family_id, items) {
    def ped_files = items.collect { entry -> entry[0].familyPed }.findAll { ped -> ped }.unique(false)
    if (ped_files.size() > 1) {
        error("Samples in the same family must reference the same familyPed value in the input samplesheet. Found ${ped_files} in family ${family_id}.")
    }
}

// Same rule as validatePedFiles, for familyPheno (required for exomiser family mode).
def validatePhenopacketFiles(family_id, items) {
    def pheno_files = items.collect { entry -> entry[0].familyPheno }.findAll { pheno -> pheno }.unique(false)
    if (pheno_files.size() > 1) {
        error("Samples in the same family must reference the same familyPheno value in the input samplesheet. Found ${pheno_files} in family ${family_id}.")
    }
}

// caller='DRAGEN_JOINT' means 'vcf' is already a family-level, jointly-called, multi-sample VCF
// -- there's no per-sample VCF to normalize/collapse/merge, so a family using it doesn't fit the
// one-row-per-sample shape the rest of the samplesheet assumes. Enforce that exactly one row
// represents the whole family in that case, rather than silently processing only one of several
// rows or mixing the two input shapes together.
def validateJointCallerExclusivity(family_id, items) {
    def callers = items.collect { entry -> entry[0].caller }
    def joint_count = callers.count { caller -> caller == 'DRAGEN_JOINT' }
    if (joint_count > 0 && items.size() > 1) {
        error("Family ${family_id} has caller='DRAGEN_JOINT' but ${items.size()} rows -- a DRAGEN_JOINT family must have exactly one row (the joint VCF already contains every sample), and cannot mix joint and per-sample rows.")
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS FROM FERLAB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def addMetadataFile(metadata_dict, destination){
    // Generate the report
    def report = _printMap(metadata_dict)

    // Write the report to the specified destination
    file(destination).write(report)
}

// Helper function to recursively print nested maps
// We could use JsonOutput.toJson() instead, but it is not stable and can cause stack overflow
// Note that, for now, lists will be printed as is even if they contain complex objects
def _printMap(map, indent =0){
        def report = ""
        map.each { key, value ->
            def prefix = "  " * indent // Indentation for nested levels
            if (value instanceof Map) {
                report += "${prefix}${key}:\n"
                report += _printMap(value, indent + 1) // Recursively process nested maps
                report += "\n"
            } else {
                report += "${prefix}${key}: ${value}\n"
            }
        }
        return report
}

def copyConfigurationFiles(config_files, destination){
    file(destination).mkdirs()
    config_files.each { config_file ->
        file(config_file).copyTo(destination)
    }
}

def copyNextflowLogFile(source, destination){
    // This log is provided to help pipeline operators easily locate the nextflow log file copies.
    // It should be shown at the command line when the pipeline is run.
    log.info "Copying nextflow log file ${source} to ${destination}"
    file(source).copyTo(destination)
}

// Infer the nextflow log file path:
// 1. Check the '-log' command-line option.
// 2. If not set, check the 'NXF_LOG_FILE' environment variable.
// 3. Default to '.nextflow.log' if neither is provided.
def getNextflowLogFilePath(command_line) {

    // Check if the -log option is present in the command line
    def tokens = parseCommandLine(command_line)
    def log_option_index = tokens.lastIndexOf("-log")
    if (log_option_index >= 0) {
        return tokens[log_option_index + 1]
    }

    // Check if the NXF_LOG_FILE environment variable is set.
    // If it is, use that value, otherwise default to '.nextflow.log'
    return System.getenv("NXF_LOG_FILE") ?: '.nextflow.log'
}

// This code takes a command-line string (e.g., "arg1 'arg2 with spaces' arg3") and
// splits it into individual tokens (arguments). It ensures that quoted arguments
// (single or double quotes) are treated as a single token, even if they contain spaces.
def parseCommandLine(command_line) {
    if (!command_line) {
        error "Command line not provided"
    }

    // The regex contains 3 groups, separated by the pipe operator '|'. The first group
    // capture arguments between double quotes, the second group captures arguments between
    // single quotes, and the third group captures unquoted arguments.
    def matcher = command_line =~ /"([^"]+)"|'([^']+)'|(\S+)/

    // We recusively iterate over the matcher to extract the tokens
    def tokens = []
    matcher.each { match ->
        tokens << (match[1] ?: match[2] ?: match[3])
    }
    return tokens
}
