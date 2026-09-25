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
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_SLIVAR } from '../modules/local/bcftools_index/main.nf'
include { BAM_VCF_DEPTH_GENOTYPE_REFINEMENT } from '../subworkflows/local/bam_vcf_depth_genotype_refinement/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { VCF_ANNOTATE_ENSEMBLVEP as VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE } from '../subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf'
include { ENSEMBLVEP_DOWNLOAD } from '../modules/nf-core/ensemblvep/download/main.nf'
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
}

workflow CNV_POST_PROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    def exomiser_local_frequency_file = params.exomiser_local_frequency_path? file(params.exomiser_local_frequency_path) : []
    def exomiser_local_frequency_index_file = params.exomiser_local_frequency_index_path
        ? file(params.exomiser_local_frequency_index_path)
        : (params.exomiser_local_frequency_path ? file("${params.exomiser_local_frequency_path}.tbi") : [])
    def exomiser_data_dir = params.exomiser_data_dir? file(params.exomiser_data_dir) : []
    def exomiser_analysis_wes_path = params.exomiser_analysis_wes? file(params.exomiser_analysis_wes) : []
    def exomiser_analysis_wgs_path = params.exomiser_analysis_wgs? file(params.exomiser_analysis_wgs) : []

    def reference_fasta = file(params.reference_fasta)
    def reference_fasta_fai = params.reference_fasta_fai? file(params.reference_fasta_fai) : file("${params.reference_fasta}.fai")

    def slivar_js = file("${projectDir}/assets/cnv-slivar-functions.js")

    // Only carries the pinned nf-core components (which predate topic channels); every local
    // module reports its versions.yml via the `versions` topic instead (collated at the end).
    def ch_versions = Channel.empty()

    //
    // VEP cache: download it fresh via ENSEMBLVEP_DOWNLOAD when download_cache is set, otherwise
    // use the pre-installed directory at vep_cache. vep_install's --SPECIES needs the cache
    // flavor baked into the species name (e.g. "homo_sapiens_merged"), unlike VEP's own --species
    // at annotation time, which always stays the plain species name (--merged/--refseq select the
    // flavor there instead, see conf/modules.config's ENSEMBLVEP_VEP block) -- matches the sibling
    // SNV pipeline's own vep_species_download (workflows/postprocessing.nf:91-100).
    //
    def vep_species_download = params.vep_annotation ? "${params.vep_species}_${params.vep_annotation}" : params.vep_species

    def vep_cache
    if (params.download_cache) {
        def ch_ensemblvep_info = Channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, vep_species_download, params.vep_cache_version])
        ENSEMBLVEP_DOWNLOAD(ch_ensemblvep_info)
        vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect().map { _meta, cache -> [cache] }.first()
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions.first())
    } else {
        vep_cache = file(params.vep_cache)
    }

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

    TRUVARI_COLLAPSE_PERSAMPLE(NORMALIZE_CNV.out.vcf.join(NORMALIZE_CNV.out.tbi))

    // truvari collapse's output isn't guaranteed sorted/indexed -- see BCFTOOLS_SORT_INDEX
    BCFTOOLS_SORT_INDEX_PERSAMPLE(TRUVARI_COLLAPSE_PERSAMPLE.out.vcf)

    //
    // STAGE 1b -- VEP annotation, single-sample route (meta.id = "familyId.sample"). Feeds
    // EXOMISER_SINGLE below with a VEP-annotated VCF, mirroring the family-level VEP+Exomiser
    // wiring (STAGE 4) but scoped to each sample's own cleaned per-sample VCF rather than the
    // family-merged one. Unlike EXOMISER_WORKFLOW's exomiser_start_from_vep toggle, there's no
    // pre-VEP option here -- EXOMISER_SINGLE always starts from VEP.
    //
    // Restricted to true solo samples (no familyPheno and no familyPed) -- a sample with either
    // already gets its individual-level annotation/prioritization needs covered by the
    // family-level route (STAGE 4+) below, so this route would otherwise be redundant for it.
    //
    def ch_persample_vep_input = BCFTOOLS_SORT_INDEX_PERSAMPLE.out.vcf
        .filter { meta, _vcf -> !meta.familyPheno && !meta.familyPed }
        .map { meta, vcf -> [meta, vcf, []] }

    VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE(
        ch_persample_vep_input,
        [[id: 'reference'], reference_fasta],
        params.vep_genome,
        params.vep_species,
        params.vep_cache_version,
        vep_cache,
        []
    )
    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE.out.versions)

    //
    // STAGE 1c -- Exomiser (single-sample mode), on each solo sample's own VEP-annotated
    // per-sample VCF (same solo-only restriction as STAGE 1b above, since ch_exomiser_single_input
    // derives from its output). Independent of family mode below, for the samples it does cover:
    // uses its own per-sample `pheno` (distinct from the family-level `familyPheno` exomiser
    // (family mode) uses).
    //
    def ch_exomiser_single_input = VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE.out.vcf_tbi
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
        .filter { meta, _vcfs, _tbis ->
            // A true family of one (sampleSize==1, no familyPed/familyPheno) has nothing for the
            // family-level route (merge/cohort-collapse/refine/VEP/exomiser/slivar) to add over
            // the persample route above (STAGE 1b/1c) -- same gate that route already uses.
            // Dropped here rather than left to fall out downstream so it doesn't pay for a second
            // VEP annotation of the same sample's CNVs for no benefit. A documented singleton
            // (sampleSize==1 but with a real familyPed/familyPheno) still goes through the "solo"
            // branch below unaffected.
            !(meta.sampleSize == 1 && !meta.familyPed && !meta.familyPheno)
        }
        .branch { meta, vcfs, tbis ->
            solo: meta.sampleSize == 1
                return [meta, vcfs[0], tbis[0]]
            family: meta.sampleSize > 1
                return [meta, vcfs, tbis]
        }

    BCFTOOLS_MERGE(ch_grouped_by_family.family)

    def ch_family_vcf = BCFTOOLS_MERGE.out.vcf
        .join(BCFTOOLS_MERGE.out.tbi)
        .mix(ch_grouped_by_family.solo)

    TRUVARI_COLLAPSE_COHORT(ch_family_vcf)

    // truvari collapse's output isn't guaranteed sorted/indexed -- see BCFTOOLS_SORT_INDEX
    BCFTOOLS_SORT_INDEX_COHORT(TRUVARI_COLLAPSE_COHORT.out.vcf)

    //
    // STAGE 3 -- Depth-based genotype refinement (mosdepth), family-level. Joint families never
    // appear in BCFTOOLS_SORT_INDEX_COHORT.out (they skipped stages 1-2 above), so this subworkflow
    // never attempts refinement for them regardless of what's in the cram-channel -- but the
    // cram-channel itself is still built from the persample branch only, for clarity.
    //
    // True solos (sampleSize==1, no familyPed/familyPheno) are excluded here too, same predicate
    // as the ch_grouped_by_family filter above: they never reach BCFTOOLS_SORT_INDEX_COHORT.out
    // (their family route was dropped entirely), so leaving their cram in this channel would hand
    // BAM_VCF_DEPTH_GENOTYPE_REFINEMENT a family id with a cram-count but no matching family VCF --
    // an orphan its internal `join(..., remainder: true)` doesn't shape-pad correctly, which
    // crashes downstream with a MissingMethodException instead of a clean skip.
    //
    BAM_VCF_DEPTH_GENOTYPE_REFINEMENT(
        BCFTOOLS_SORT_INDEX_COHORT.out.vcf.join(BCFTOOLS_SORT_INDEX_COHORT.out.tbi),
        ch_samplesheet_branched.persample
            .filter { meta, _vcf, _cram -> !(meta.sampleSize == 1 && !meta.familyPed && !meta.familyPheno) }
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

    BCFTOOLS_SORT_INDEX_JOINT(NORMALIZE_CNV_JOINT.out.vcf)

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

    // Slivar's own output isn't indexed -- this is the pipeline's last file, so index it here
    // rather than leave that as a manual step for whoever consumes it next. Rewriting INFO tags
    // doesn't reorder records, so this is a plain re-index, not a re-sort (same reasoning as
    // BCFTOOLS_INDEX's other use after DEPTH_GENOTYPE_REFINE).
    BCFTOOLS_INDEX_SLIVAR(SLIVAR_EXPR.out.vcf)

    //
    // Collate and save software versions
    //
    topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'cnv-post-processing_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    versions = ch_versions   // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
