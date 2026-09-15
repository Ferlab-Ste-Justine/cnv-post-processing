//
// Depth-based genotype refinement: resolves missing (./.) genotypes at cohort-collapsed CNV
// sites into either an explicit hom-ref call (when depth confirms it) or leaves them missing
// (schema.txt's "missing -> hom-ref vs missing -> no-data" goal).
//
// Depth is inherently per-sample, so this fans a family-level cohort VCF back out to one
// mosdepth run per sample (each against the *same* shared site set, using that sample's own
// cram), computes each sample's per-site depth ratio (region depth / that sample's genome-wide
// mean depth), then feeds all of a family's ratio files together into one rewrite step.
//
// Originally designed around duphold, which does this same fold-change computation internally
// and annotates the VCF directly -- switched to mosdepth (2026-09-04) after duphold's ~6-year
// unmaintained biocontainer failed to decode a real CRAM 3.1 file, which is now the default
// output format for current samtools/DRAGEN. mosdepth is actively maintained and doesn't hit
// that problem, at the cost of implementing the ratio/rewrite logic here instead of getting it
// for free (see bin/refine_genotypes.py).
//
// A family is only refined if EVERY one of its samples has a cram in the samplesheet -- partial
// coverage skips refinement for the whole family (with a warning) rather than silently dropping a
// sample's genotype column.
//

include { VCF_TO_BED            } from '../../../modules/local/vcf_to_bed/main'
include { MOSDEPTH_RATIO        } from '../../../modules/local/mosdepth_ratio/main'
include { DEPTH_GENOTYPE_REFINE } from '../../../modules/local/depth_genotype_refine/main'
include { BCFTOOLS_INDEX        } from '../../../modules/local/bcftools_index/main'

workflow BAM_VCF_DEPTH_GENOTYPE_REFINEMENT {

    take:
    ch_family_vcf  // [meta(id:familyId, samples:[s1,s2,...]), vcf, tbi]
    ch_sample_cram // [meta(id:sample, family:familyId), cram, crai] -- pre-filtered to samples that have a cram
    ch_fasta       // path: reference genome fasta
    ch_fasta_fai   // path: reference genome fasta index

    main:
    def ch_versions = Channel.empty()

    //
    // Only refine families where every sample has a cram
    //
    def ch_cram_counts = ch_sample_cram
        .map { meta, _cram, _crai -> [meta.family, 1] }
        .groupTuple()
        .map { family, ones -> [family, ones.size()] }

    def ch_family_branched = ch_family_vcf
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .join(ch_cram_counts, remainder: true)
        .map { _family, meta, vcf, tbi, cram_count ->
            def has_all_crams = (cram_count ?: 0) == meta.samples.size()
            if (!has_all_crams) {
                log.warn("Skipping depth-based genotype refinement for family '${meta.id}': not every sample has a cram (alignment file) provided in the samplesheet.")
            }
            [meta, vcf, tbi, has_all_crams]
        }
        .branch { meta, vcf, tbi, has_all_crams ->
            refine: has_all_crams
                return [meta, vcf, tbi]
            skip: true
                return [meta, vcf, tbi]
        }

    //
    // Candidate site regions (shared across a family's samples), once per family
    //
    VCF_TO_BED(ch_family_branched.refine)
    ch_versions = ch_versions.mix(VCF_TO_BED.out.versions)

    //
    // Fan-out: family sites BED x per-sample cram -> N per-sample mosdepth runs. Each family has
    // its own sites BED, so it travels as part of the main per-invocation tuple (not a broadcast
    // value) -- otherwise every family's samples would incorrectly race for whichever family's
    // BED happened to be first in the channel.
    //
    def ch_sites_bed_by_family = VCF_TO_BED.out.bed.map { meta, bed -> [meta.id, bed] }

    def ch_mosdepth_input = ch_family_branched.refine
        .flatMap { meta, _vcf, _tbi -> meta.samples.collect { s -> [meta.id, [id: s, family: meta.id]] } }
        .combine(ch_sites_bed_by_family, by: 0)
        .map { _family, sample_meta, bed -> [sample_meta, bed] }
        .join(ch_sample_cram.map { meta, cram, crai -> [[id: meta.id, family: meta.family], cram, crai] })
        .map { sample_meta, bed, cram, crai -> [sample_meta, cram, crai, bed] }

    MOSDEPTH_RATIO(ch_mosdepth_input, ch_fasta, ch_fasta_fai)
    ch_versions = ch_versions.mix(MOSDEPTH_RATIO.out.versions)

    //
    // Fan-in: gather every sample's ratio BED for a family alongside that family's VCF
    //
    def ch_refine_input = ch_family_branched.refine
        .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
        .combine(
            MOSDEPTH_RATIO.out.ratio_bed
                .map { meta, bed -> [meta.family, bed] }
                .groupTuple(),
            by: 0
        )
        .map { _family, meta, vcf, tbi, beds -> [meta, vcf, tbi, beds] }

    DEPTH_GENOTYPE_REFINE(ch_refine_input)
    ch_versions = ch_versions.mix(DEPTH_GENOTYPE_REFINE.out.versions)

    // Rewriting genotypes doesn't reorder records, so this is a plain re-index, not a re-sort.
    BCFTOOLS_INDEX(DEPTH_GENOTYPE_REFINE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX.out.versions)

    def ch_output = BCFTOOLS_INDEX.out.vcf_tbi
        .mix(ch_family_branched.skip)

    emit:
    vcf      = ch_output // [meta, vcf, tbi] -- refined where possible, original cohort VCF otherwise
    versions = ch_versions
}
