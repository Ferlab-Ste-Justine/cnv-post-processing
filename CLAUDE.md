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

Nextflow version range: `!>=23.10.1` (`manifest.nextflowVersion` in `nextflow.config`). CI
(`.github/workflows/ci-nf-test.yml`) matrixes `23.10.1`, `24.10.5`, and `latest-stable` — a real
constraint in practice, see "Pinned dependencies" below.

## High-level pipeline flow

The entry workflow is `main.nf` → `CNV_POST_PROCESSING` (`workflows/cnv_post_processing.nf`).
Unlike the SNV pipeline there's no `--step`/`--tools` gating — every samplesheet row runs the
route applicable to it, decided per-row/per-family purely from which optional columns
(`cram`, `pheno`, `familyPheno`, `familyPed`) are populated.

**Family route** (a family with `sampleSize > 1`, or a solo row that documents `familyPed`/`familyPheno`):

1. `NORMALIZE_CNV` (bcftools norm/view) → `TRUVARI_COLLAPSE_PERSAMPLE`, per sample
2. `BCFTOOLS_MERGE` across a family's per-sample VCFs → `TRUVARI_COLLAPSE_COHORT`. Families
   supplied as an already jointly-called `caller=DRAGEN_JOINT` multi-sample VCF skip merge+collapse
   entirely (lightweight prep only: drop non-variant records, split multiallelics) and rejoin at
   step 4
3. `BAM_VCF_DEPTH_GENOTYPE_REFINEMENT` (mosdepth-based) — only if _every_ sample in the family has
   a `cram`; partial coverage skips refinement for the whole family with a warning
4. `VCF_ANNOTATE_ENSEMBLVEP`, family-level
5. `EXOMISER_WORKFLOW` — only for families with `familyPheno`
6. `SLIVAR_EXPR` mode-of-inheritance classification — only for families with `familyPed`. Its
   output is the pipeline's last file, so it's the only stage that gets an explicit indexing step
   of its own (`BCFTOOLS_INDEX_SLIVAR`) rather than relying on a downstream stage to need/produce
   one.

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
nf-test.config                             # nf-test runner config (profile "test,docker")
workflows/cnv_post_processing.nf           # Main CNV_POST_PROCESSING workflow (both routes above) + EXOMISER_WORKFLOW
subworkflows/local/                        # bam_vcf_depth_genotype_refinement, utils_nfcore_cnv_post_processing_pipeline
subworkflows/nf-core/                      # utils_nextflow_pipeline, utils_nfcore_pipeline, utils_nfschema_plugin, vcf_annotate_ensemblvep
modules/local/                             # normalize_cnv, truvari_collapse, bcftools_{merge,sort_index,index}, vcf_to_bed, mosdepth_ratio, depth_genotype_refine, slivar_expr, exomiser
modules/nf-core/                           # ensemblvep/{vep,download}, tabix/tabix
conf/                                      # base.config, modules.config, test.config, test_full.config
assets/                                    # samplesheet.csv, schema_input.json, cnv-slivar-functions.js, exomiser/ (default analysis YAMLs)
bin/refine_genotypes.py                    # Depth-ratio-based genotype rewrite logic used by DEPTH_GENOTYPE_REFINE
docs/                                      # usage.md, output.md, reference_data.md, schema.txt (ASCII flow diagram)
scripts/                                   # run-test-suite.sh, run-smoke-tests.sh (see below)
```

## How to run

Typical invocation (from the README):

```bash
nextflow -c application.config run Ferlab-Ste-Justine/cnv-post-processing \
    -r v1.0.0 \
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
`scripts/run-test-suite.sh` bundles the pre-push gate: nf-test dry-run, full nf-test suite, and
`nf-core pipelines lint --release`. All 10 local modules (`modules/local/*`) have their own module
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

CI lint workflow: `.github/workflows/linting.yml` (nf-core lint + prettier/whitespace pre-commit
checks, config in `.pre-commit-config.yaml`). Run locally with:

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
family-mode Exomiser), `familyPed` (identical across a family's rows, enables slivar).
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
  fanning back out later (e.g. depth-refinement's per-sample mosdepth runs).
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
  doesn't exist under Nextflow 24.10.5 either) until a decision is made to bump the pipeline's
  supported Nextflow version — see `git log` on `subworkflows/nf-core/utils_nextflow_pipeline/` for
  context when that happens.
- **Per-process resources** live in `nextflow.config` / `conf/base.config` under
  `process { withLabel/withName: ... }`; per-process publish paths and `ext.args` are in
  `conf/modules.config` (one exception: the Exomiser container is pinned directly in
  `modules/local/exomiser/main.nf`, not in `conf/modules.config`).
- **Adding an nf-core module:** use `nf-core modules install <tool>` so `modules.json` stays
  consistent. Local-only logic goes under `modules/local/`.
- **Schema and params stay in sync.** `nextflow_schema.json` drives the `nf-schema` plugin. When
  adding a param, update both `nextflow.config` defaults and the schema.
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
- Changelog and version history: `CHANGELOG.md`.
- PR template: `.github/PULL_REQUEST_TEMPLATE.md`.
