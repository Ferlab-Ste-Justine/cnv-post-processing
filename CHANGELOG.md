# Ferlab-Ste-Justine/cnv-post-processing: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Changed`

- Pipeline version bumped to 2.0.0 (`manifest.version`, `.nf-core.yml`)
- Local-module software version collection moved to Nextflow topic channels (`topic: versions`, collated via `channel.topic("versions")`); only the pinned nf-core VEP components are still mixed into `ch_versions` manually
- CI Nextflow matrix is now 24.10.5 and 25.10.4, plus a non-blocking `latest-everything` canary in the nf-test workflow
- CI: every GitHub Action is pinned to a full commit SHA (release tag in a trailing comment), reusing Post-processing-Pipeline's SHAs where it pins the same action; `commit_lint.yml`'s `actions/checkout` moved from v2 to v4.2.2 in the process. The lint job is also pinned to Nextflow 24.10.5 instead of whatever release is latest
- [#XXX](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/XXX) BIOINFO-232: `EXOMISER` (and its `EXOMISER_SINGLE` alias) now uses the official `docker.io/exomiser/exomiser-cli:14.0.0-bash` image instead of the custom `ferlabcrsj/exomiser:2.8.1` one. The exomiser version reported in `versions.yml` is now read from the image's `exomiser-core` jar filename instead of a baked-in `/EXOMISER_VERSION.txt`, and the image's `/bin/bash` entrypoint is cleared for both processes in the `docker` and `podman` profiles.
- CI workflows harmonized with Post-processing-Pipeline: `ci-nf-test.yml` is replaced by nf-core's sharded `nf-test.yml` (affected tests split into up to 2 shards, a single `confirm-pass` result, full suite on releases and manual runs), `ci.yml` is renamed `ci-full-run.yml`, `linting.yml` now matches the reference (Python 3.13, PR and release triggers, `--release` lint only on PRs to `main`), and the S3 test-data download moved into the `.github/actions/copy-test-data` composite action. `.nf-core.yml` pins nf-core tools 4.1.0
- CI: new `ci-pr-title-lint.yml` checks that PR titles follow `<type>: <TICKET-123> <description>` (from Post-processing-Pipeline). `commit_lint.yml` stays for now; its removal is pending repository-setting changes, see `docs/journal.md`
- slivar mode-of-inheritance rules rewritten as a table (`assets/cnv-slivar-functions.js`): one call per affected child and record (`cnv_moi`), and, as in the SNV pipeline, one `--family-expr` per tag (`has_moi_call`) that tags the family with each of its children's calls, valued with its affected members; the calls are `denovo`, `denovo_suspect`, `inherited_maternal`, `inherited_paternal`, `recessive`, `recessive_lowconf`, `possible_denovo`, `possible_inherited_maternal`, `possible_inherited_paternal`, `inheritance_unresolved`, `mendelian_inconsistent` and `unknown_cn`, with the SNV pipeline's X-linked names on non-PAR chrX (`x_linked_recessive_denovo`, `x_linked_dominant_denovo`, `x_linked_recessive`). The call looks at the child and its parents only: there is no family-level tag, and who is affected elsewhere in the family is left to the reader (see `docs/journal.md`). These replace `de_novo_candidate`, `dominant_inherited`, `recessive_candidate`, `candidate` and `ambiguous`. Expected copy number is chromosome- and sex-aware (GRCh38 PARs; records crossing a PAR boundary, and chrM, aren't modelled), and parent of origin (`po_mother`, `po_father`, `po_both`, `po_ambiguous`, plus `po_denovo` and `po_possible_denovo`, SNV-style) is computed per sample (`--sample-expr`) by the same transmission rules the MoI tags use, for every child. Tested by `assets/cnv-slivar-moi-harness.js` against a golden file of named scenarios (`assets/moi_scenarios.expected.tsv`), which replaces `assets/tests/cnv_slivar_functions.test.js`
- CI: new `ci-slivar-harness.yml` runs the slivar rules harness on every PR and push to `main`, and fails if the committed scenario files are out of date
- The pipeline now stops at startup when a family has a `familyPed` and `vep_genome` isn't `GRCh38`, since slivar's chrX/chrY model uses GRCh38 PAR coordinates
- slivar parent-of-origin tags renamed: INFO `parent_of_origin_maternal`/`parent_of_origin_paternal`/`parent_of_origin_ambiguous` are now `po_mother`/`po_father`/`po_ambiguous`, and their functions in `assets/cnv-slivar-functions.js` (`moi_cnv_po_maternal`/`moi_cnv_po_paternal`/`moi_cnv_po_ambiguous`) carry the same names. Update any downstream filter on the old INFO keys
- slivar: a parent who is in the VCF but has no usable CN now counts as absent, as the SNV pipeline and Radiant treat a parent's no-call, so the child is called like a duo (e.g. kid CN 1, mother CN 2, father CN missing: `possible_denovo` instead of `unknown_cn`). `unknown_cn` is kept for a child without CN, or a carrier none of whose needed parents has a CN
- slivar: two parent-of-origin tags from the SNV pipeline, `po_possible_mother` and `po_possible_father`, for a child with one parent unknown (absent, or no CN) whose other parent carries a CNV it may have passed (duo D1: kid CN 1, mother CN 1). When that parent surely passed it, the child gets `po_mother` / `po_father` (duo D2: kid CN 0, mother CN 1), as in the SNV pipeline and Radiant. Before, these children got no parent-of-origin tag
- slivar: the de novo mode-of-inheritance tags are spelled `denovo`, like the parent-of-origin tags (`po_denovo`): `denovo`, `denovo_suspect`, `possible_denovo`, `x_linked_recessive_denovo`, `x_linked_dominant_denovo` (were `de_novo`, `de_novo_suspect`, `possible_de_novo`, `x_linked_recessive_de_novo`, `x_linked_dominant_de_novo` on this branch). The two X-linked names are no longer Radiant's `transmission_mode` values verbatim (see `docs/journal.md`)
- **Breaking:** The pipeline now checks each family's PED against the sample names in the family VCF before the slivar step (ported from snv-post-processing BIOINFO-240), and stops with a message naming the problem and its fix if (1) a VCF sample has no PED row, or (2) a sequenced sample's father or mother is named but isn't in the VCF. slivar exits normally in both cases. In the first it left the sample out of the family, and a PED matching none of the VCF's names gave no tags at all. In the second, when the absent parent was the father, slivar also dropped the child's link to the mother, so the child silently got no mode-of-inheritance or parent-of-origin tag. PEDs that named unsequenced parents and ran before now fail; set those parents to `0`. Rows for relatives who weren't sequenced are still allowed. The new `--skip_ped_validation` parameter (default `false`) turns these errors into warnings and lets the run continue, e.g. for tests. See `validatePedAgainstVcf()` and `vcfSampleNames()` in `subworkflows/local/utils_nfcore_cnv_post_processing_pipeline/main.nf`, both with nf-tests. `docs/usage.md` now also says the PED uses the VCFs' sample names, not the samplesheet's `sample` column

### `Fixed`

- A family's merged VCF now always lists its samples in name order. Before, `BCFTOOLS_MERGE` got the per-sample VCFs in whatever order their tasks finished, so the sample columns of the family VCF, and of every file built from it, could change from one run to the next, and `-resume` could re-run the merge and every family-level step after it
- `conf/modules.config`'s VEP `ext.args` no longer fails to parse under Nextflow's strict config syntax (the default from 26.04, "Variable declarations cannot be mixed with config statements"); its local variables now live inside an `ext.args` closure. Resolved args are unchanged for all three cache flavors

### `Dependencies`

- Supported Nextflow range is now `>=24.10.5, <26.0.0` (23.10.1 dropped)

## v1.0.0 - 2025-05-01

Initial release of Ferlab-Ste-Justine/cnv-post-processing, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#6](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/6) Add exomiser params
- [#7](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/7) Add exomiser into the workflow
- [#8](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/8) Add exomiser_outdir option
- [#9](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/9) Add extra files in pipeline_info folder to improve reproducibility and traceability

### `Fixed`

- [#10](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/10) Support sample id given as integer

### `Changed`

- [#5](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/pull/5) Replace autogenerated samplesheet with pipeline samplesheet

### `Fixed`

### `Dependencies`

### `Deprecated`
