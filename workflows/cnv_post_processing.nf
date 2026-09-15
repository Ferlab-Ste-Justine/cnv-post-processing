/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { EXOMISER } from '../modules/local/exomiser/main.nf'
include { EXOMISER as EXOMISER_SINGLE } from '../modules/local/exomiser/main.nf'
include { NORMALIZE_CNV } from '../modules/local/normalize_cnv/main.nf'
include { NORMALIZE_CNV as NORMALIZE_CNV_JOINT } from '../modules/local/normalize_cnv/main.nf'
include { TRUVARI_COLLAPSE as TRUVARI_COLLAPSE_PERSAMPLE } from '../modules/local/truvari_collapse/main.nf'
include { TRUVARI_COLLAPSE as TRUVARI_COLLAPSE_COHORT    } from '../modules/local/truvari_collapse/main.nf'
include { BCFTOOLS_SORT_INDEX as BCFTOOLS_SORT_INDEX_PERSAMPLE } from '../modules/local/bcftools_sort_index/main.nf'
include { BCFTOOLS_SORT_INDEX as BCFTOOLS_SORT_INDEX_COHORT    } from '../modules/local/bcftools_sort_index/main.nf'
include { BCFTOOLS_SORT_INDEX as BCFTOOLS_SORT_INDEX_JOINT     } from '../modules/local/bcftools_sort_index/main.nf'
include { BCFTOOLS_MERGE } from '../modules/local/bcftools_merge/main.nf'
include { BAM_VCF_DEPTH_GENOTYPE_REFINEMENT } from '../subworkflows/local/bam_vcf_depth_genotype_refinement/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { SLIVAR_EXPR } from '../modules/local/slivar_expr/main.nf'
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
    // Family mode: ch_input is family-level (multi-sample VCF + family-level phenopacket),
    // matching the sibling SNV pipeline's EXOMISER wiring (workflows/postprocessing.nf:308-357).
    def ch_input_for_exomiser = ch_input
        .filter { meta, _vcf, _tbi ->
            if (!meta.familyPheno) {
                log.warn("Skipping exomiser for family '${meta.id}': no familyPheno (phenopacket) provided in the samplesheet.")
            }
            meta.familyPheno
        }
        .map{
            meta, vcf, tbi ->
                def pheno = file(meta.familyPheno)
                def analysis_file = meta.sequencingType == "WES"? analysis_wes_path : analysis_wgs_path
                [meta, vcf, tbi, pheno, analysis_file]
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

    def reference_fasta = file(params.reference_fasta)
    def reference_fasta_fai = params.reference_fasta_fai? file(params.reference_fasta_fai) : file("${params.reference_fasta}.fai")

    def vep_cache = file(params.vep_cache)
    def slivar_js = file("${projectDir}/assets/cnv-slivar-functions.js")

    def ch_versions = Channel.empty()

    //
    // caller='DRAGEN_JOINT' rows are already a family-level, jointly-called, multi-sample VCF --
    // there's no per-sample VCF to normalize/collapse/merge, and joint calling doesn't have the
    // post-hoc-merge missing-genotype problem depth-based refinement exists to work around. Those
    // families skip stages 1-3 entirely (validated to be exactly one row per family in
    // PIPELINE_INITIALISATION's validateJointCallerExclusivity) and only get the same lightweight
    // prep (drop non-variant records, split multiallelic sites -- DRAGEN's own joint VCFs do carry
    // multiallelic <DEL>,<DUP> sites, unlike its per-sample output) before rejoining the regular
    // families at the "family-level VCF ready for VEP" point below.
    //
    def ch_samplesheet_branched = ch_samplesheet
        .branch { meta, _vcf, _cram ->
            joint: meta.caller == 'DRAGEN_JOINT'
            persample: true
        }

    //
    // STAGE 1 -- Per-sample normalize + truvari collapse (meta.id = "familyId.sample")
    //
    NORMALIZE_CNV(ch_samplesheet_branched.persample.map { meta, vcf, _cram -> [meta, vcf] })
    ch_versions = ch_versions.mix(NORMALIZE_CNV.out.versions)

    TRUVARI_COLLAPSE_PERSAMPLE(NORMALIZE_CNV.out.vcf.join(NORMALIZE_CNV.out.tbi))
    ch_versions = ch_versions.mix(TRUVARI_COLLAPSE_PERSAMPLE.out.versions)

    // truvari collapse's output isn't guaranteed sorted/indexed -- see BCFTOOLS_SORT_INDEX
    BCFTOOLS_SORT_INDEX_PERSAMPLE(TRUVARI_COLLAPSE_PERSAMPLE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_INDEX_PERSAMPLE.out.versions)

    //
    // STAGE 1b -- Exomiser (single-sample mode), on each sample's own cleaned per-sample VCF,
    // pre-VEP -- matches this pipeline's original (pre-family-redesign) behavior. Independent of
    // family mode below: every sample gets its own individual-level prioritization in addition to
    // its family's joint one, using its own per-sample `pheno` (distinct from the family-level
    // `familyPheno` exomiser (family mode) uses).
    //
    def ch_exomiser_single_input = BCFTOOLS_SORT_INDEX_PERSAMPLE.out.vcf
        .join(BCFTOOLS_SORT_INDEX_PERSAMPLE.out.tbi)
        .filter { meta, _vcf, _tbi ->
            if (!meta.pheno) {
                log.warn("Skipping exomiser (single-sample mode) for sample '${meta.id}': no pheno (phenopacket) provided in the samplesheet.")
                return false
            }
            true
        }
        .map { meta, vcf, tbi ->
            def pheno = file(meta.pheno)
            def analysis_file = meta.sequencingType == "WES" ? exomiser_analysis_wes_path : exomiser_analysis_wgs_path
            [meta, vcf, tbi, pheno, analysis_file]
        }

    EXOMISER_SINGLE(
        ch_exomiser_single_input,
        exomiser_data_dir,
        params.exomiser_genome,
        params.exomiser_data_version,
        [exomiser_local_frequency_file, exomiser_local_frequency_index_file],
        params.exomiser_remm_version ? [params.exomiser_remm_version, params.exomiser_remm_filename] : ["", ""],
        params.exomiser_cadd_version ? [params.exomiser_cadd_version, params.exomiser_cadd_snv_filename, params.exomiser_cadd_indel_filename] : ["", "", ""]
    )
    ch_versions = ch_versions.mix(EXOMISER_SINGLE.out.versions)

    //
    // STAGE 2 -- Family merge + cohort truvari collapse (fan-in: meta.id "familyId.sample" -> familyId)
    //
    def ch_grouped_by_family = BCFTOOLS_SORT_INDEX_PERSAMPLE.out.vcf
        .join(BCFTOOLS_SORT_INDEX_PERSAMPLE.out.tbi)
        .map { meta, vcf, tbi -> tuple(groupKey(meta.familyId, meta.sampleSize), meta, vcf, tbi) }
        .groupTuple()
        .map { _familyId, metas, vcfs, tbis ->
            // now that samples are grouped together, drop the per-sample-only meta fields and
            // reset id to the family level. Keep the sample-id list (needed to fan back out to
            // per-sample steps later, e.g. depth-based genotype refinement).
            def family_meta = metas[0].findAll { entry -> !["sample", "id", "caller", "cram"].contains(entry.key) }
            family_meta = family_meta + [id: family_meta.familyId, samples: metas.collect { it.sample }]
            [family_meta, vcfs.flatten(), tbis.flatten()]
        }
        .branch { meta, vcfs, tbis ->
            solo: meta.sampleSize == 1
                return [meta, vcfs[0], tbis[0]]
            family: meta.sampleSize > 1
                return [meta, vcfs, tbis]
        }

    BCFTOOLS_MERGE(ch_grouped_by_family.family)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    def ch_family_vcf = BCFTOOLS_MERGE.out.vcf
        .join(BCFTOOLS_MERGE.out.tbi)
        .mix(ch_grouped_by_family.solo)

    TRUVARI_COLLAPSE_COHORT(ch_family_vcf)
    ch_versions = ch_versions.mix(TRUVARI_COLLAPSE_COHORT.out.versions)

    // truvari collapse's output isn't guaranteed sorted/indexed -- see BCFTOOLS_SORT_INDEX
    BCFTOOLS_SORT_INDEX_COHORT(TRUVARI_COLLAPSE_COHORT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_INDEX_COHORT.out.versions)

    //
    // STAGE 3 -- Depth-based genotype refinement (mosdepth), family-level. Joint families never
    // appear in BCFTOOLS_SORT_INDEX_COHORT.out (they skipped stages 1-2 above), so this subworkflow
    // never attempts refinement for them regardless of what's in the cram-channel -- but the
    // cram-channel itself is still built from the persample branch only, for clarity.
    //
    BAM_VCF_DEPTH_GENOTYPE_REFINEMENT(
        BCFTOOLS_SORT_INDEX_COHORT.out.vcf.join(BCFTOOLS_SORT_INDEX_COHORT.out.tbi),
        ch_samplesheet_branched.persample
            .filter { meta, _vcf, cram ->
                if (!cram) {
                    log.warn("Skipping depth-based genotype refinement for sample '${meta.id}': no cram (alignment file) provided in the samplesheet.")
                }
                cram
            }
            .map { meta, _vcf, cram ->
                def crai = file("${cram}.crai")
                def bai = file("${cram}.bai")
                def index = crai.exists() ? crai : (bai.exists() ? bai : [])
                [[id: meta.sample, family: meta.familyId], file(cram), index]
            },
        reference_fasta,
        reference_fasta_fai
    )
    ch_versions = ch_versions.mix(BAM_VCF_DEPTH_GENOTYPE_REFINEMENT.out.versions)

    //
    // STAGE 3b -- Joint-called families: lightweight prep only (drop non-variant records, split
    // multiallelic sites), no collapse/merge/refinement needed. Rejoins the regular families right
    // below at the "family-level VCF ready for VEP" point.
    //
    def ch_joint_family_meta = ch_samplesheet_branched.joint
        .map { meta, vcf, _cram ->
            // Unlike the regular family_meta construction below (which drops 'caller' since it
            // becomes ambiguous once multiple samples with potentially different callers merge),
            // NORMALIZE_CNV_JOINT reuses NORMALIZE_CNV's own caller guard, so 'caller' must stay.
            def family_meta = meta.findAll { entry -> !["sample", "id", "cram"].contains(entry.key) }
            family_meta = family_meta + [id: family_meta.familyId, samples: [meta.sample]]
            [family_meta, vcf]
        }

    NORMALIZE_CNV_JOINT(ch_joint_family_meta)
    ch_versions = ch_versions.mix(NORMALIZE_CNV_JOINT.out.versions)

    BCFTOOLS_SORT_INDEX_JOINT(NORMALIZE_CNV_JOINT.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT_INDEX_JOINT.out.versions)

    def ch_family_vcf_for_vep = BAM_VCF_DEPTH_GENOTYPE_REFINEMENT.out.vcf
        .mix(BCFTOOLS_SORT_INDEX_JOINT.out.vcf.join(BCFTOOLS_SORT_INDEX_JOINT.out.tbi))

    //
    // STAGE 4 -- VEP annotation, family-level (AnnotSV/ClassifyCNV/vcfanno dropped per team
    // update 2026-09-04 -- VEP's own CSQ output is the final annotated VCF, nothing to merge)
    //
    // Pinned to the sibling SNV pipeline's exact vendored commits (vcf_annotate_ensemblvep
    // cfd937a6, ensemblvep/vep 6e3585d9, tabix/tabix 66665215) rather than nf-core/modules'
    // current HEAD: the freshly-vendored ensemblvep/vep doesn't compile under Nextflow 24.10.5
    // ("Variable prefix already defined in the process scope") -- confirmed directly, and
    // confirmed the sibling's older pinned commit (predating that upstream refactor) doesn't
    // have the issue. `cache` is a plain path here, not a [meta,cache] tuple (verified directly
    // against this exact module revision -- the subworkflow's own inline comment for ch_cache is
    // stale/inaccurate).
    VCF_ANNOTATE_ENSEMBLVEP(
        ch_family_vcf_for_vep.map { meta, vcf, _tbi -> [meta, vcf, []] },
        [[id: 'reference'], reference_fasta],
        params.vep_genome,
        params.vep_species,
        params.vep_cache_version,
        vep_cache,
        []
    )
    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)

    // Exomiser can run off either the VEP-annotated VCF or the pre-VEP family VCF (both are
    // family-level, indexed), matching the sibling's exomiser_start_from_vep toggle.
    def ch_exomiser_input = params.exomiser_start_from_vep
        ? VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
        : ch_family_vcf_for_vep

    EXOMISER_WORKFLOW(
        ch_exomiser_input,
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
    // STAGE 6 -- Slivar mode-of-inheritance classification, on the VEP-annotated family VCF.
    // The moi_cnv_*(fam) functions in assets/cnv-slivar-functions.js generalize to any family
    // shape the PED describes (trios, duos, singletons, sibships) via slivar's --family-expr --
    // each MoI's own structural gate (has_aff_parent, no_parents_in_fam, etc.) decides whether a
    // given family provides enough evidence for it, same as the sibling SNV pipeline's own
    // moi_*(fam) design. No family-shape pre-check needed here, matching that pipeline's
    // convention -- only familyPed presence is required.
    //
    def ch_slivar_input = VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi
        .filter { meta, _vcf, _tbi ->
            if (!meta.familyPed) {
                log.warn("Skipping slivar inheritance classification for family '${meta.id}': no familyPed (PED file) provided in the samplesheet.")
                return false
            }
            true
        }
        .map { meta, vcf, tbi -> [meta, vcf, tbi, file(meta.familyPed)] }

    SLIVAR_EXPR(ch_slivar_input, slivar_js)
    ch_versions = ch_versions.mix(SLIVAR_EXPR.out.versions)

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
