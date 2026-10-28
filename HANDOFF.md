# Handoff: per-sample original call ID (FORMAT/OID) (BIOINFO-214)

Branch: `feat/BIOINFO-214-proband-representative-cnv`, branched from `feat/BIOINFO-214` at `cad99b4`.
Status: implemented, tested, committed and pushed. No PR opened yet.

This file is a temporary note for continuing the work on another machine (amd64). Delete it before the PR is merged.

## History of the branch

1. `d0c8afc`: the first approach, **proband priority**. It made the proband's call each cohort-collapse cluster's representative by rewriting QUAL and running truvari with `--keep maxqual`.
2. `bfb1e47`: reverted that approach. It changed which coordinates a record shows, but the other samples' original IDs were still lost.
3. The next commit: the current approach, **`FORMAT/OID`**, suggested by a colleague.

The branch name still says "proband-representative". Rename it, or pick a PR title that describes `OID`, e.g. `feat: BIOINFO-214 Keep each sample's original CNV call ID in FORMAT/OID`. PRs are squash-merged, and the title must match `<type>: <TICKET-123> <description>`.

## Problem

`truvari collapse` keeps one record per cluster, and that record keeps one `ID`: the leftmost call, often a parent's. After `TRUVARI_COLLAPSE_COHORT`, nothing in the record said which original call each other sample contributed.

That link existed only in truvari's removed VCF: `INFO/MatchId` on each absorbed call points to `INFO/CollapseId` on the kept record. Nothing downstream reads that file.

## Solution

`NORMALIZE_CNV` runs on one sample at a time, before `BCFTOOLS_MERGE`. It now copies each record's `ID` into that sample's `FORMAT/OID`.

When truvari folds a call into a kept record, it copies **all** of that call's FORMAT fields into samples that are still ref or missing in the kept record. So `OID` travels with the genotype. This was verified on a toy trio and on the real NA12878 trio.

`:` separates FORMAT values, and DRAGEN IDs contain it (`DRAGEN:LOSS:chr22:17287644-17294633`). So `OID` percent-encodes it as `%3A`, as VCF 4.3 requires, and `%` becomes `%25` first. To decode, reverse those replacements.

bcftools can't copy `ID` into FORMAT itself: `-c FMT/OID:=ID` fails with "Cannot rename". That's why the module uses a small awk step.

The full rationale is in `docs/journal.md`, under "Each sample keeps its own call ID in FORMAT/OID".

## Files changed (vs `cad99b4`)

| File                                                             | Change                                                                                                                                                                   |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `bin/add_original_id.awk`                                        | New. Works on a single-sample VCF from stdin: adds the `OID` header and sets `FORMAT/OID` = encoded `ID`. Exits 1 if the VCF doesn't have exactly one sample.            |
| `modules/local/normalize_cnv/main.nf`                            | When `task.ext.original_id` is set, pipes `\| bcftools view - \| add_original_id.awk` after norm and before the final `bcftools view -O z`.                              |
| `conf/modules.config`                                            | `withName: NORMALIZE_CNV { ext.original_id = true }`. Not set for `NORMALIZE_CNV_JOINT`: joint VCFs skip the collapses, and they're multi-sample, which the awk refuses. |
| `modules/local/normalize_cnv/tests/nextflow_original_id.config`  | New test config that sets `ext.original_id`.                                                                                                                             |
| `modules/local/normalize_cnv/tests/main.nf.test`                 | The per-sample DRAGEN test now uses that config, and asserts the `OID` header plus `OID == ID.replace(":", "%3A")` on every record.                                      |
| `CHANGELOG.md`, `CLAUDE.md`, `docs/output.md`, `docs/journal.md` | Documentation.                                                                                                                                                           |

## Verification done (on the Apple Silicon machine)

- **`normalize_cnv` module tests:** the per-sample and GATK-guard tests pass.
  - The `DRAGEN_JOINT` test fails only because `data-test/vcf/1463_NA12878_NA12891_NA12892.cnv.vcf.gz` is missing from S3 (see follow-ups). It was already failing before this change.
- **`truvari_collapse` module tests:** all pass. They were back to their original 3 tests after the revert.
- **Full pipeline on the NA12878 trio (`fam1`) plus `ind1`, with the samplesheet that has no joint row:** completed successfully.
  - `FORMAT/OID` is present on every record of `merged/`, `truvari/`, `depth_refined/`, `ensemblvep/` (family and single-sample) and `slivar/`.
  - Samples that depth refinement sets to `0/0` keep `OID=.`, as they should.
  - Example: at `chr22:17.29M`, the record's `ID` is NA12891's call `17287644-17294633`, and the `OID` of NA12878 and NA12892 is `17288730-17294633`.
  - A baseline run with `ext.original_id = false` gave a slivar VCF identical to this one once `FORMAT/OID` is removed, so nothing else changed.
  - Exomiser reports 0 variants in both runs. That's how the test data behaves, not an effect of `OID`.
- **Not verified on amd64 or under the real truvari biocontainer:** on arm64 the truvari biocontainer crashes, so the local runs used a custom arm64 image (`~/.nextflow/config` override). The awk step runs in the bcftools container, which is busybox awk.

## First steps on the amd64 machine

```bash
git fetch && git checkout feat/BIOINFO-214-proband-representative-cnv
export NXF_FILE_ROOT=$PWD
aws s3 cp s3://ferlab-public-dataset/nextflow/cnv-post-processing/V2/data-test data-test --recursive

nf-test test modules/local/normalize_cnv modules/local/truvari_collapse --profile +docker
```

You don't need any of the arm64 workarounds on amd64: no `~/.nextflow/config` override, no `DOCKER_DEFAULT_PLATFORM`, no custom truvari image.

To check `OID` end to end, run the pipeline test. It needs the full S3 dataset (see follow-ups). Without it, build a local samplesheet: copy `assets/samplesheet.csv` into `data-test/`, make its paths absolute, and drop the `fam1_joint` row. Then run:

```bash
nextflow run . -profile test,docker --input data-test/samplesheet_nojoint.csv --outdir results
```

## Known issues / follow-ups

1. **The S3 test dataset is incomplete.** `s3://.../cnv-post-processing/V2/data-test` has no `samplesheet.csv` and no `vcf/1463_NA12878_NA12891_NA12892.cnv.vcf.gz`. Because of that, `tests/default.nf.test` and the `DRAGEN_JOINT` `normalize_cnv` test fail. This predates this branch, and CI will hit it too.
   - `tests/default.nf.test.snap` has no VCF content md5s, so `OID` shouldn't require a snapshot update.
2. **One `OID` per sample per record.** If a sample has two calls in one cluster (e.g. a fragmented CNV), only the first call's `OID` is copied. The second survives only in the removed VCF.
3. **DEL/DUP collapse into each other (not fixed).** DRAGEN writes `SVTYPE=CNV` for both, and truvari treats both as the same unknown type, so a DUP can absorb a DEL at the same locus. A likely fix is to set `SVTYPE` from `ALT` before the collapse.
4. **Downstream consumers (Radiant, etc.)** have to decode `%3A` back to `:` to match `OID` against the original caller IDs.
5. **Delete `HANDOFF.md`** before merging.
