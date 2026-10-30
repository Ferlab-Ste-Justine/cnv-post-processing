# CLAUDE.md

This file gives Claude Code the context it needs to work effectively in this repository.

## Project overview

`Ferlab-Ste-Justine/cnv-post-processing` is a Nextflow DSL2 pipeline for family-based
post-processing of Copy Number Variant (CNV) calls: it normalizes and collapses per-sample CNV
VCFs, merges samples into a family-level VCF (or reuses an already jointly-called one), resolves
missing genotypes from alignment depth, annotates with VEP, prioritizes with Exomiser, and tags
mode of inheritance with slivar. Standalone samples with no family/pedigree information get their
own lightweight per-sample VEP + Exomiser route instead of the family-level one.

The repo follows nf-core conventions (scaffolded with nf-core tools v3.2.0, see `.nf-core.yml`).
It is _not_ a published nf-core pipeline (`is_nfcore: false`).

Sibling of `Ferlab-Ste-Justine/Post-processing-Pipeline` (the SNV equivalent) — this pipeline
reuses that one's pinned `vcf_annotate_ensemblvep` subworkflow commits and mirrors its slivar
mode-of-inheritance architecture (per-sample segregation predicates + family-shape structural
gates), adapted for copy-number genotypes (`FORMAT/CN` instead of `GT`/`AD`/`GQ`).

Nextflow version range: `!>=24.10.5, <26.0.0` (`manifest.nextflowVersion` in `nextflow.config`;
fleet-wide decision, 23.10.1 dropped because topic-channel version collection needs >=24.04.0).
CI (`.github/workflows/nf-test.yml`, `.github/workflows/ci-full-run.yml`) matrixes `24.10.5` and
`25.10.4`; `nf-test.yml` also runs `latest-everything` as a non-blocking
(`continue-on-error`) forward-compat canary that may legitimately fail against the `<26.0.0` cap.
The 24.10.5 floor is a real constraint in practice, see "Pinned dependencies" below.

## High-level pipeline flow

The entry workflow is `main.nf` → `CNV_POST_PROCESSING` (`workflows/cnv_post_processing.nf`).
Unlike the SNV pipeline there's no `--step`/`--tools` gating — every samplesheet row runs the
route applicable to it, decided per-row/per-family purely from which optional columns
(`cram`, `pheno`, `familyPheno`, `familyPed`) are populated.

**Family route** (a family with `sampleSize > 1`, or a solo row that documents `familyPed`/`familyPheno`):

1. `NORMALIZE_CNV` (bcftools norm/view; also copies each call's ID into `FORMAT/OID` via
   `bin/add_original_id.awk`) → `TRUVARI_COLLAPSE_PERSAMPLE`, per sample
2. `BCFTOOLS_MERGE` across a family's per-sample VCFs → `TRUVARI_COLLAPSE_COHORT` →
   `RESTORE_COLLAPSED_FORMAT` (`bin/restore_collapsed_format.py`: copies back the FORMAT of
   absorbed calls truvari reads as unset, DRAGEN's `./1` and haploid `1`, from its removed VCF). Families
   supplied as an already jointly-called `caller=DRAGEN_JOINT` multi-sample VCF skip merge+collapse
   entirely (lightweight prep only: drop non-variant records, split multiallelics) and rejoin at
   step 4
3. `BAM_VCF_DEPTH_GENOTYPE_REFINEMENT` (mosdepth-based) — only if _every_ sample in the family has
   a `cram`; partial coverage skips refinement for the whole family with a warning. A resolved
   `./.` gets GT hom-ref _and_ `FORMAT/CN` = the expected CN (slivar reads CN, not GT), sex-aware
   on chrX/chrY from `familyPed`, mirroring `get_expected_cn()` in the slivar rules (keep the PAR
   coordinates in `bin/refine_genotypes.py` in sync with them). Records under 5 kb
   (`ext.min_size`) are never refined: depth there can't tell normal from a missed loss, which
   made false `denovo` calls
4. `VCF_ANNOTATE_ENSEMBLVEP`, family-level
5. `SLIVAR_EXPR` mode-of-inheritance classification — only for families with `familyPed`. Its
   output is the family's final VCF, so it gets an explicit indexing step of its own
   (`BCFTOOLS_INDEX_SLIVAR`) rather than relying on a downstream stage to need/produce one.
   Before it, in the workflow's slivar-input `.map` (not in a subworkflow, so tests can still feed
   slivar inconsistent PEDs), the run stops if the family VCF declares no `FORMAT/CN`
   (`vcfFormatIds()`), and `validatePedAgainstVcf` checks the PED against the family VCF's
   `#CHROM` names and stops the run on a mismatch, or only warns under `--skip_ped_validation`.
6. `JOIN_PROBAND_CNV` (`bin/join_proband_cnv.py`) — one VCF per proband (`pedAffectedChildren`:
   the PED's affected children), joining the proband's own input-VCF calls to the slivar records
   they're in through the decoded `FORMAT/OID`. Per-sample families only: `DRAGEN_JOINT` ones have
   no per-sample VCF or `OID`. The PED names VCF samples, so the proband's input VCF is found by
   `vcfSampleNames()`, not the samplesheet's `sample` column.
7. `EXOMISER_WORKFLOW` — only for families with `familyPheno`. Runs on the `JOIN_PROBAND_CNV` VCF
   of the phenopacket's proband (`phenopacketProbandId()`: `proband.subject.id`), output
   `<familyId>.<proband>.exomiser.*`; a family with no such VCF (`DRAGEN_JOINT`, no `familyPed`,
   proband not a sequenced affected PED child) falls back to the family VCF (pre-VEP, or VEP's
   under `--exomiser_start_from_vep`) through a `join(..., remainder: true)`. Unlike the SNV
   pipeline, which always runs Exomiser on the family VCF, in parallel with slivar.

**Solo route** (true singletons: `sampleSize==1`, no `familyPed`, no `familyPheno` — the family
route would add nothing over this for them, so they're dropped from `ch_grouped_by_family`
entirely rather than paying for a second VEP run):

1. `NORMALIZE_CNV` + `TRUVARI_COLLAPSE_PERSAMPLE` (same processes as the family route, shared)
2. `VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE`
3. `EXOMISER_SINGLE` — only if the sample row has a `pheno`

See `docs/schema.txt` for a visual ASCII diagram of this flow and `docs/images/cnv_post_processing_workflow.png`
for the nf-core-style graphic (the PNG predates the slivar stage).

## Repository layout

```
main.nf                                    # Entry point — PIPELINE_INITIALISATION, CNV_POST_PROCESSING, PIPELINE_COMPLETION
nextflow.config                            # Params, profiles, per-process resources, manifest
nextflow_schema.json                       # Authoritative parameter schema (use this, not the README)
nf-test.config                             # nf-test runner config (default profile "test"; add docker via --profile)
workflows/cnv_post_processing.nf           # Main CNV_POST_PROCESSING workflow (both routes above) + EXOMISER_WORKFLOW
subworkflows/local/                        # bam_vcf_depth_genotype_refinement, utils_nfcore_cnv_post_processing_pipeline
subworkflows/nf-core/                      # utils_nextflow_pipeline, utils_nfcore_pipeline, utils_nfschema_plugin, vcf_annotate_ensemblvep
modules/local/                             # normalize_cnv, truvari_collapse, bcftools_{merge,sort_index,index}, vcf_to_bed, mosdepth_ratio, depth_genotype_refine, slivar_expr, join_proband_cnv, restore_collapsed_format, exomiser
modules/nf-core/                           # ensemblvep/{vep,download}, tabix/tabix
conf/                                      # base.config, modules.config, test.config, test_full.config
assets/                                    # samplesheet.csv, schema_input.json, exomiser/ (default analysis YAMLs), cnv-slivar-functions.js (slivar rules) + cnv-slivar-moi-harness.js (its test bench; golden file moi_scenarios.expected.tsv)
bin/add_original_id.awk                    # Copies each call's ID into FORMAT/OID in NORMALIZE_CNV
bin/join_proband_cnv.py                    # Joins a proband's own calls to the slivar VCF through FORMAT/OID (JOIN_PROBAND_CNV)
bin/restore_collapsed_format.py            # Restores the FORMAT truvari's cohort collapse drops for ./1 and haploid calls (RESTORE_COLLAPSED_FORMAT)
bin/refine_genotypes.py                    # Depth-ratio-based genotype rewrite logic used by DEPTH_GENOTYPE_REFINE
docs/                                      # usage.md, output.md, reference_data.md, schema.txt (ASCII flow diagram)
scripts/                                   # run-test-suite.sh, run-smoke-tests.sh (see below)
```

## How to run

Typical invocation (from the README):

```bash
nextflow -c application.config run Ferlab-Ste-Justine/cnv-post-processing \
    -r v2.0.0 \
    --input samplesheet.csv \
    --outdir results \
    -params-file params.json
```

### Test dataset

Test data is expected under `data-test/` at the repo root; it is not checked in (`.gitignore`).
Sync it first:

```bash
aws s3 cp s3://ferlab-public-dataset/nextflow/cnv-post-processing/V2/data-test data-test --recursive
```

`data-test/samplesheet.csv` covers, in one file, a `DRAGEN` trio (`fam1`), the same trio as a
`DRAGEN_JOINT` multi-sample VCF (`fam1_joint`), and a true standalone sample (`ind1`) — exercising
the family route, the joint-caller shortcut, and the solo route in a single test run.

In CI, `nf-test.yml` and `ci-full-run.yml` both download it through the
`.github/actions/copy-test-data` composite action, the only place CI defines the S3 path. When the
dataset version changes, update it there, in `scripts/run-smoke-tests.sh`, and here.

### Test profile

```bash
nextflow run . -profile test,docker --outdir results
```

Make sure Docker is installed and running. To clean up outputs, run `nextflow clean -f`.

For a bundle of manual, eyeballed smoke runs (debug+test profile, plain test profile, and a VEP
refseq-cache-flavor check, each to its own persistent `--outdir`), run `scripts/run-smoke-tests.sh`
— it verifies `data-test/` is synced and Docker is running before starting, and leaves existing
output directories in place rather than wiping them (diff against a previous run, or
`nextflow clean -f`/remove them yourself for a clean slate).

### Tests (nf-test)

```bash
export NXF_FILE_ROOT=$PWD          # needed so nf-test finds test files
nf-test test --dryRun              # syntax check only, no execution
nf-test test --profile test,docker # full pipeline test
nf-test test modules/local/exomiser
```

nf-core upstream module/subworkflow tests are excluded via the `ignore` glob in `nf-test.config`.
`scripts/run-test-suite.sh` bundles the pre-push gate, harmonized with Post-processing-Pipeline's:
nf-core CLI version pin check, pre-commit, the slivar rules harness, a `-preview` launch under the
declared Nextflow floor, nf-test dry-run, full nf-test suite, and `nf-core pipelines lint --release`. All 12 local modules (`modules/local/*`) have their own module
test under `modules/local/<name>/tests/`.

**nf-test gotchas** (both hit, and fixed, while writing `modules/local/truvari_collapse/tests/`;
every other module test in this repo already avoids both, whether by luck or a previous author
already having learned this):

- A param set via a config file (e.g. `params.modules_testdata_base_path`) resolves fine in a
  test's `when` block but can come back `null` in its `then` block, throwing a
  `NullPointerException` -- a known nf-test bug (askimed/nf-test#288), not something to work around
  with `params { ... }` reassignment (doesn't help, same issue). Never reference `params.*` inside
  a `then` block -- use `${projectDir}`-relative paths instead (works reliably in both `when` and
  `then`).
- Calling `.linesGzip` from _inside_ a `with(process.out) { ... }` closure throws
  `NullPointerException: Cannot invoke "GZIPInputStream.close()" because "gzip" is null` --
  confirmed via a real `--debug` stack trace straight into
  `com.askimed.nf.test.lang.extensions.GlobalMethods.with(...)` wrapping the failing
  `PathExtension.getLinesGzip(...)` call. This isn't about which path is being read or whether the
  channel reference is bare vs fully-qualified (both fail identically) -- it's `with()` itself:
  being inside that closure at all breaks `.linesGzip`. Close `with(process.out) { }` right after
  the ordinary `.size()`/meta/filename assertions, then do every `.linesGzip` read afterward in
  plain `then` scope using fully-qualified `process.out.vcf...` (exactly the structure
  `modules/local/bcftools_sort_index/tests/main.nf.test` and `bcftools_index`'s already use).

### Linting

CI workflows live in `.github/workflows/`, harmonized with Post-processing-Pipeline (the fleet
reference): `linting.yml` (pre-commit + nf-core lint; plain lint on PRs to other branches,
`--release` lint on PRs to `main`), `nf-test.yml` (nf-core's sharded design: `nf-test-changes`
counts the tests affected since `HEAD^` and splits them into at most 2 shards, each shard runs on
every matrixed Nextflow version through the composite actions in `.github/actions/`, and
`confirm-pass` gives one result excluding `latest-everything`; releases and manual runs test
everything), `ci-full-run.yml` (a plain `nextflow run . -profile test,docker`, which doesn't depend
on nf-test's change detection), `ci-pr-title-lint.yml` (the PR title must look like
`<type>: <TICKET-123> <description>`, as in Post-processing-Pipeline), `ci-slivar-harness.yml` (runs
`assets/cnv-slivar-moi-harness.js` against its golden file and checks the generated `moi_scenarios.*` /
`trio_lookup_table.tsv` are committed up to date) and `commit_lint.yml` (Ferlab's
`action-commit-lint`, still the required check on `main`; its removal is pending, see
`docs/journal.md`). Files the nf-core template defines keep
their nf-core names (lint only recognizes those); this repo's own workflows get a `ci-` prefix.
Run lint locally with:

```bash
pre-commit run --all-files
nf-core pipelines lint --release
```

`.nf-core.yml` carries deliberate lint overrides:

- Several nf-core-template files are intentionally absent (e.g. `CODE_OF_CONDUCT.md`, nf-core
  logos, AWS CI workflows) because this is a Ferlab pipeline, not a published nf-core one. Don't
  reintroduce those files; update `.nf-core.yml` instead if lint behavior needs to change.
- `nextflow_config: [custom_config, ...]` silences the "outdated lines for loading custom
  profiles" failure — the `nf-core/configs` `includeConfig` lines are deliberately commented out
  in `nextflow.config` (see the comment there) for security reasons, so no external
  config/profile is ever auto-loaded from a remote source. `included_configs` (a related, separate
  lint test for the pipeline-specific institutional config line) is left as a visible warning
  rather than silenced, as a reminder of the same deviation.

## Samplesheet format

Authoritative schema: `assets/schema_input.json`. Required: `familyId`, `sample`,
`sequencingType` (`WES`/`WGS`), `caller` (`DRAGEN`/`DRAGEN_JOINT`/`GATK` — only `DRAGEN` and
`DRAGEN_JOINT` are currently supported by normalization), `vcf`. Optional: `cram` (enables
depth-based genotype refinement), `pheno` (per-sample phenopacket, enables single-sample
Exomiser), `familyPheno` (family-level phenopacket, identical across a family's rows, enables
family-mode Exomiser), `familyPed` (identical across a family's rows, enables slivar; its IDs
are the VCFs' `#CHROM` sample names, not the samplesheet's `sample` column, e.g. `fam1_joint`).
`caller=DRAGEN_JOINT` families must have exactly one row and no per-sample rows (enforced by
`validateJointCallerExclusivity` in `PIPELINE_INITIALISATION`).

## Working in this codebase

A few patterns and gotchas worth knowing before editing:

- **Channel shape convention.** Most VCF channels carry `[meta, vcf, tbi]`. After a `BCFTOOLS_*`
  call, the index is emitted separately and joined back: `out.vcf.join(out.tbi)`.
- **`meta` grows and shrinks across stages.** Per-sample meta (`id: "familyId.sample"`) carries
  `familyId`, `sampleSize`, `familyPed`, `familyPheno`, `pheno` from the samplesheet. At the
  family-merge fan-in (stage 2 in the flow above), sample-only fields (`sample`, `id`, `caller`,
  `cram`) are dropped and `id` resets to the family level, keeping a `samples: [...]` list for
  fanning back out later (e.g. depth-refinement's per-sample mosdepth runs). The fan-in sorts a
  family's samples by name first: `groupTuple` keeps arrival order, which varies between runs, and
  `bcftools merge` writes its sample columns in input order. Any new fan-in feeding an
  order-sensitive tool, or putting a per-sample list in `meta`, needs the same sort. (A `path`
  list alone doesn't: Nextflow hashes it regardless of order.)
- **The solo/family split is a filter, not just a branch.** True solos are dropped from
  `ch_grouped_by_family` entirely (a `.filter{...}` in `workflows/cnv_post_processing.nf`) _before_
  the `.branch{ solo / family }` that follows. Any other channel built from
  `ch_samplesheet_branched.persample` for a family-level step (e.g. the depth-refinement cram
  channel) needs that same exclusion applied independently — otherwise it hands that subworkflow a
  family id with no matching family VCF, an orphan that `BAM_VCF_DEPTH_GENOTYPE_REFINEMENT`'s
  internal `join(..., remainder: true)` doesn't null-pad correctly, crashing with a
  `MissingMethodException` instead of cleanly skipping. This exact class of bug was hit and fixed
  here (2026-09) by mirroring the solo-exclusion filter on that cram channel too.
- **`DRAGEN_JOINT` families take a shortcut**, not the normal per-sample path — see High-level
  flow. `NORMALIZE_CNV_JOINT` reuses `NORMALIZE_CNV`'s own caller guard, so unlike the regular
  family-merge meta construction, `caller` is deliberately kept in `family_meta` there.
- **Pinned dependencies — don't blindly `nf-core modules/subworkflows update`.** The
  `vcf_annotate_ensemblvep` subworkflow (and its `ensemblvep/vep`, `tabix/tabix` modules) are
  pinned to the SNV pipeline's exact vendored commits rather than nf-core/modules' current HEAD:
  the freshly-vendored `ensemblvep/vep` doesn't compile under Nextflow 24.10.5 ("Variable prefix
  already defined in the process scope"), confirmed directly. `utils_nextflow_pipeline` is
  similarly held back (its latest commit references `nextflow.script.types.VersionNumber`, which
  doesn't exist under Nextflow 24.10.5 either). The fleet's Nextflow-version decision has since
  been made (floor `24.10.5`, see top of this file), and since that floor still lacks both, these
  pins stay until the floor moves again — see `git log` on
  `subworkflows/nf-core/utils_nextflow_pipeline/` for context when that happens.
- **Software versions collection uses topic channels.** Every local module emits
  `path("versions.yml"), emit: versions, topic: versions`, and `CNV_POST_PROCESSING` collates
  `channel.topic("versions")` once at the end (harmonized with `quality-control-pipeline`).
  `main.nf` sets `nextflow.preview.topic = true` only when `nextflow.version < 25.04.0` — it's
  version-guarded because topics are GA from 25.04.0 and setting the flag unconditionally is a
  hard error on Nextflow >=25.10. The `emit: versions` is kept alongside the topic for the module
  tests: nf-test's module wrapper never runs `main.nf`, so the topic is inert there and the tests
  read `process.out.versions`. `ch_versions` now only carries the three pinned nf-core components
  that predate topics (`ENSEMBLVEP_DOWNLOAD`, `VCF_ANNOTATE_ENSEMBLVEP_PERSAMPLE`,
  `VCF_ANNOTATE_ENSEMBLVEP`), mixed in manually. New local processes must add `topic: versions`,
  not a `ch_versions.mix()` line.
- **Per-process resources** live in `nextflow.config` / `conf/base.config` under
  `process { withLabel/withName: ... }`; per-process publish paths and `ext.args` are in
  `conf/modules.config` (one exception: the Exomiser container is pinned directly in
  `modules/local/exomiser/main.nf`, not in `conf/modules.config`).
- **Adding an nf-core module:** use `nf-core modules install <tool>` so `modules.json` stays
  consistent. Local-only logic goes under `modules/local/`.
- **Schema and params stay in sync.** `nextflow_schema.json` drives the `nf-schema` plugin. When
  adding a param, update both `nextflow.config` defaults and the schema.
- **slivar rules: `assets/cnv-slivar-functions.js`.** One computation, two per-child outputs,
  with a map of each at the top of the file:
  - `get_transmission_outcome()` enumerates every way one copy from each parent adds up to the
    child's CN (a parent's copies are its haplotypes, CN 3 = `[2,1]`; a missing parent, i.e. not
    in the PED or no CN, could pass 0, 1 or 2) and names each passed copy normal, a hit (changed
    the child's way) or opposite. The outcome (`new`, `maternal`, `paternal`, `both`,
    `ambiguous`, the `duo_*` ones, or `undecided`) plus the child's deviation picks one
    `TRIO_TABLE` row, which gives both the MoI call and the parent of origin, so the two always
    agree. There are no parent labels (`same`/`half`) any more: a parent homozygous for a
    deletion (CN 0) is handled like any other (139 records of the real trio were
    `inheritance_unresolved` before). A parent with CN 4 or more can't be split, so its child is
    `undecided`: `is_mendelian_possible()` then gives `inheritance_unresolved` or
    `mendelian_inconsistent`.
  - `cnv_moi(sample, variant)` gives each affected carrier child its row's call, at most one per
    child. There is no family-level logic (no `dominant`, no segregation check across sibs): who is
    affected elsewhere in the family is left to whoever reads the tags. On non-PAR chrX, two calls
    take an X-linked name (`X_LINKED_CALLS`): a son's or a daughter's `denovo` → `x_linked_denovo`,
    a daughter's `recessive` → `x_linked_recessive`. A single inherited X copy keeps
    `inherited_*`, a son's like a daughter's: the SNV pipeline's sex-specific names
    (`x_linked_recessive_de_novo`, `x_linked_dominant_de_novo`, a son's `x_linked_recessive`) read
    dominant or recessive off the zygosity and who is affected, which doesn't fit a son's gain
    (CN 2) and which the calls don't look at. As in the SNV pipeline, the 14 mode-of-inheritance tags
    are `--family-expr` in `conf/modules.config`: `has_moi_call(fam, variant, tag)` is true when
    any child of the family gets that call, so a family with several affected children can carry
    several tags, and slivar values them with the family's affected members.
  - `cn_po_origin(sample, variant)` drives the 8 parent-of-origin tags, named as in the SNV
    pipeline (`po_mother`, `po_father`, `po_both`, `po_ambiguous`, `po_denovo`,
    `po_possible_denovo`, `po_possible_mother`, `po_possible_father`), as `--sample-expr` valued
    with exactly the children it resolves: its row's parent of origin, for every child, affected or
    not, duos included. `po_denovo` needs both parents known and normal (a parent carrying a CNV
    that can't explain the child's gets no tag: that's `mendelian_inconsistent`). In a duo, the
    known parent gets `po_mother`/`po_father` when it surely passed the CNV, `po_possible_*` when
    it may have. A founder never gets one (a male founder on chrX would otherwise see his
    "father", who passes nothing, as a known parent).
  - Expected CN is chromosome- and sex-aware: 2 on autosomes and PAR, on non-PAR chrX 1 for a male
    and 2 for a female, on non-PAR chrY 1 for a male. chrM, PAR-straddling records, unknown sex on
    non-PAR chrX/chrY and a female on non-PAR chrY aren't modelled (no tag). A parent who passes
    nothing (a father's X to a son, any mother's Y) is read as normal and their CN isn't needed
    (`does_parent_pass_nothing()`), even when absent from the PED or without CN: biology, not
    missing data (decided, harness X14-X16). The PARs are hard-coded GRCh38, hence the startup guard on
    `vep_genome` in `validateInputSamplesheet`; the PED sex column must be right.
  - The file is ES5, for slivar's JS engine (duktape). Reading an attribute slivar didn't set
    (`s.mom` on a founder, `INFO.END` on a record without one) throws there, and slivar then
    writes the record with no tag, silently: test parents with `has_parent()` (`in`) and read the
    CNV end from `variant.stop`. `has_cn()` tests `"CN" in s` the same way, and the workflow
    stops before `SLIVAR_EXPR` when the family VCF declares no `FORMAT/CN` (`vcfFormatIds()`).
  - `TRIO_TABLE` is hand-aligned under `// prettier-ignore`; keep it aligned.
    The worked-example notes inside `cnv_moi()`, `trio_call()` and `get_trio_row()` are kept on purpose.
  - Tests: `node assets/cnv-slivar-moi-harness.js` runs every named scenario against the golden file
    `assets/moi_scenarios.expected.tsv` (re-run with `--accept` once a change is understood) plus
    200,000 random families whose fake samples throw on unset attributes like slivar's, and checks that
    `conf/modules.config` has one slivar expression per call and parent of origin. It
    regenerates `assets/moi_scenarios.md` / `.tsv` and `trio_lookup_table.tsv`, which must be
    committed; `ci-slivar-harness.yml` checks both. The `SLIVAR_EXPR` module test runs real slivar:
    fixture row IDs name the harness scenario they reproduce, and `fam1_joint` covers the real
    trio's 55 chrX/chrY records. The pipeline test only reaches slivar on chr22 (VEP drops other
    chromosomes with the chr22-only test cache).
  - Open questions on the rules are in `docs/journal.md`.
- **`--no-version` on every bcftools view/norm/merge/annotate call.** bcftools appends a
  `##bcftools_<cmd>Command=...; Date=<now>` header by default, which breaks byte-for-byte
  reproducibility between otherwise-identical runs — every bcftools invocation in this pipeline's
  modules passes `--no-version` for that reason (`bcftools sort` has no such flag since it never
  adds this header in the first place).

## Outputs

`--outdir` is required. VEP and Exomiser outputs can be split to separate directories via
`--vep_outdir` and `--exomiser_outdir`. Only named/override publish paths are written by default;
intermediate per-sample normalize/collapse files are published only under `--publish_all`. Run
metadata (configs, timeline/report/trace/dag, collated software versions) is written to
`${outdir}/pipeline_info/`.

See `docs/output.md` for the full output layout (verified current as of this file's writing —
check it hasn't drifted from `conf/modules.config` before trusting it blindly on a later date).

## Pointers

- Parameter documentation: `nextflow_schema.json` (authoritative). README/usage docs intentionally
  avoid duplicating parameter details.
- Samplesheet schema: `assets/schema_input.json`.
- Reference data setup (Exomiser data dir/analysis files, VEP cache/reference fasta): `docs/reference_data.md`.
- Output layout: `docs/output.md`.
- Visual flow diagrams: `docs/schema.txt` (ASCII), `docs/images/cnv_post_processing_workflow.png`.
- slivar rules: the maps at the top of `assets/cnv-slivar-functions.js`; every named scenario's result in
  `assets/moi_scenarios.md`; an overview in `docs/cnv-slivar.txt`.
- Changelog and version history: `CHANGELOG.md`.
- PR template: `.github/PULL_REQUEST_TEMPLATE.md`.
