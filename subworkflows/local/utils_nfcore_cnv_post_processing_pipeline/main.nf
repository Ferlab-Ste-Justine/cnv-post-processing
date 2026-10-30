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
    // help/help_full/show_hidden/cli_typecast are new required inputs as of the synced
    // subworkflow version (2026-09) -- this pipeline has never wired up a --help CLI flag, so
    // they're passed as false/empty here to keep that unchanged rather than adding the feature
    // as a side effect of the sync.
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        false,
        false,
        false,
        "",
        "",
        "nextflow run . -profile <docker/singularity/...> --input samplesheet.csv --outdir <OUTDIR>",
        false
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
    // A familyPed means slivar will classify this family, and assets/cnv-slivar-functions.js
    // decides expected copy number on chrX/chrY from hard-coded GRCh38 PAR coordinates
    // (CNV_PAR_X / CNV_PAR_Y). On any other build those coordinates are wrong, so every
    // sex-chromosome call near a PAR would be silently misclassified -- fail before anything
    // runs instead.
    if (meta.familyPed && params.vep_genome != 'GRCh38') {
        error("Family '${meta.familyId}' has a familyPed, so slivar will classify it, but vep_genome is '${params.vep_genome}': slivar's chrX/chrY expected-copy-number model (assets/cnv-slivar-functions.js) only knows the GRCh38 pseudoautosomal regions. Use GRCh38 data, or remove familyPed to skip slivar.")
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

// Checks a family's PED against the samples of the family VCF slivar is about to read (port of
// snv-post-processing BIOINFO-240). slivar silently mishandles a PED that doesn't match its VCF
// and still exits 0: a VCF sample without a PED row is left out of the family (a PED matching no
// VCF sample gives no tags at all), and a child whose named father isn't in the VCF loses the link
// to its mother too. The CNV rules need a parent in the family, so those children get no tag.
// Two rules: every VCF sample has a PED row, and a sequenced sample's named father/mother is in
// the VCF (an unsequenced parent is written as 0). Rows for relatives not in the VCF are fine.
// Stops the run, or only warns when warn_only (--skip_ped_validation) is set.
def validatePedAgainstVcf(family_id, vcf, ped, warn_only = false) {
    def vcf_samples = vcfSampleNames(vcf)
    def rows = ped.readLines()
        .findAll { line -> !line.startsWith('#') }
        .collect { line -> line.trim().split(/\s+/) }
        .findAll { cols -> cols.size() >= 4 }
    def ped_samples = rows.collect { cols -> cols[1] }
    def problems = []

    def without_row = vcf_samples.findAll { sample -> !(sample in ped_samples) }
    if (without_row) {
        problems << "VCF samples without a PED row: ${without_row.join(', ')} (the PED must use the VCF's sample IDs)"
    }
    rows.findAll { cols -> cols[1] in vcf_samples }.each { cols ->
        [father: cols[2], mother: cols[3]].each { role, parent ->
            if (parent != '0' && !(parent in vcf_samples)) {
                problems << "${cols[1]}'s ${role} '${parent}' isn't in the VCF (write an unsequenced parent as 0)"
            }
        }
    }

    if (problems) {
        def message = ["Family '${family_id}': the PED file doesn't match the family VCF (${ped}):"]
        message += problems.collect { problem -> "  - ${problem}" }
        message << "Samples in the VCF: ${vcf_samples.join(', ')}."
        message << "slivar would leave out or misread these family members, so inheritance and parental-origin tags would be missing or wrong."
        if (warn_only) {
            message << "--skip_ped_validation is set, so slivar runs anyway."
            log.warn(message.join('\n'))
        } else {
            message << "To run anyway, with this as a warning, set --skip_ped_validation."
            error(message.join('\n'))
        }
    }
}

// Sample names from the VCF header's #CHROM line (bgzipped or plain).
def vcfSampleNames(vcf) {
    def input = vcf.toString().endsWith('.gz') ? new java.util.zip.GZIPInputStream(vcf.newInputStream()) : vcf.newInputStream()
    return input.withReader { reader ->
        def line
        while ((line = reader.readLine()) != null) {
            if (line.startsWith('#CHROM')) {
                return line.split('\t').drop(9) as List
            }
            if (!line.startsWith('#')) {
                break
            }
        }
        return []
    }
}

// The probands of a family for JOIN_PROBAND_CNV: the PED's affected (phenotype 2) samples with a
// father or mother named. A founder isn't one, even if affected. Several affected sibs are all
// returned.
def pedAffectedChildren(ped) {
    return ped.readLines()
        .findAll { line -> !line.startsWith('#') }
        .collect { line -> line.trim().split(/\s+/) }
        .findAll { cols -> cols.size() >= 6 && cols[5] == '2' && (cols[2] != '0' || cols[3] != '0') }
        .collect { cols -> cols[1] }
}

// The proband a family phenopacket (YAML or JSON) is about: its proband.subject.id, or null when
// it has none. Family-mode Exomiser runs on this child's JOIN_PROBAND_CNV VCF only: another
// affected child's VCF would be scored against this child's phenotypes.
def phenopacketProbandId(pheno) {
    def doc = new org.yaml.snakeyaml.Yaml().load(pheno.text)
    return doc instanceof Map ? doc.proband?.subject?.id : null
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
