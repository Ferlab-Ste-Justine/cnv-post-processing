# Handoff: proband-priority cohort collapse (BIOINFO-214)

Branch: `feat/BIOINFO-214-proband-representative-cnv`, branched from `feat/BIOINFO-214` at `cad99b4`.
Status: implemented, tested, committed and pushed. No PR opened yet.

This file is a temporary note for continuing the work on another machine (amd64). Delete it before the PR is merged.

## Goal

In `TRUVARI_COLLAPSE_COHORT`, make the proband's CNV the representative of every collapse cluster it belongs to, instead of truvari's default (the leftmost call, often a parent's).

**Why:** truvari collapse compares each call only with the cluster's representative (we don't pass `--chain`). With the default `--keep first`, a parent's leftmost call can pull in both the proband's call and the other parent's call even when those two don't match each other. The proband's row then claims the other parent carries the CNV, which corrupts the inheritance tags downstream (slivar).

The full rationale is in `docs/journal.md`, under "The cohort collapse keeps the proband's call".

## How it works

truvari can't prefer a sample, so we use QUAL as the lever:

1. `bin/proband_priority.py mark` copies the real QUAL into `INFO/ORIG_QUAL`. It then sets QUAL to `1` on records where any proband has a non-ref GT, and to `0` on all others.
2. truvari collapse runs with `--keep maxqual`. Ties fall back to truvari's leftmost rule, so clusters without a proband call behave exactly as before.
3. `bin/proband_priority.py restore` puts the real QUAL back and removes `ORIG_QUAL` from both `collapsed` and `removed` VCFs.

Probands are the PED's **affected children**: phenotype `2` with a father or mother named. All affected sibs count. This is computed by `pedAffectedChildren()`.

A rejected alternative was swapping in the proband's call after a normal collapse. The cluster members would still have been matched against the parent's representative, so the bug would remain.

## Files changed

| File | Change |
|---|---|
| `bin/proband_priority.py` | New. Uses pysam. `mark` and `restore` subcommands. |
| `modules/local/truvari_collapse/main.nf` | When `task.ext.proband_priority && meta.probands`: runs `mark`, then collapse with `--keep maxqual`, then `restore`. Otherwise unchanged. |
| `conf/modules.config` | `TRUVARI_COLLAPSE_COHORT`: `ext.proband_priority = true`. Don't add `--keep` to its `ext.args`. |
| `subworkflows/local/utils_nfcore_cnv_post_processing_pipeline/main.nf` | New `pedAffectedChildren(ped)`. |
| `workflows/cnv_post_processing.nf` | Adds `meta.probands` on input to `TRUVARI_COLLAPSE_COHORT`. Strips it again before `BCFTOOLS_SORT_INDEX_COHORT`, so later tasks keep the same meta and `-resume` hashes. |
| `modules/local/truvari_collapse/tests/data/trio_proband.vcf.gz` (+ `.tbi`) | New fixture: DAD/KID/MOM on chr22, two clusters. |
| `modules/local/truvari_collapse/tests/nextflow_cohort_proband.config` | New test config: cohort args plus `ext.proband_priority = true`. |
| `modules/local/truvari_collapse/tests/main.nf.test` | 2 new tests on `trio_proband`. One uses the default collapse and shows the bug. The other uses priority and shows the fix. |
| `subworkflows/.../tests/main.function.pedAffectedChildren.nf.test` | New. Tests `valid.ped` and `famQ.ped`. |
| `CHANGELOG.md`, `CLAUDE.md`, `docs/output.md`, `docs/journal.md` | Documentation. |

## Verification done (on the Apple Silicon machine)

- **truvari_collapse module tests:** all 5 pass, including the 2 new ones and the NA12878 real-data test.
- **`pedAffectedChildren` function tests:** both pass.
- **Full pipeline on the real NA12878 trio (`fam1`):** completed, including VEP, Exomiser and slivar.
  - The `TRUVARI_COLLAPSE_COHORT (fam1)` log shows `Marked 22 proband record(s) for NA12878`, and the command line includes `--keep maxqual` and `restore`.
  - Compared with a default collapse of the same merged VCF, 2 of the 3 clusters now keep NA12878's coordinates instead of the father's: `LOSS chr22:17.29M` and `LOSS chr22:32.14M`. The third cluster (`GAIN chr22:23.71M`) was already NA12878's.
  - QUALs are restored (28, 36, 27), and there is no `ORIG_QUAL` in the published VCF.
- **Not verified:** the script under the real truvari biocontainer (`truvari:5.4.0`). On arm64 that image crashes (SIGILL) even with Rosetta, so the local run used a custom arm64 image built with `pip install truvari==5.4.0`. On amd64, the first thing to do is run the tests with the real container (see below).

## First steps on the amd64 machine

```bash
git fetch && git checkout feat/BIOINFO-214-proband-representative-cnv
export NXF_FILE_ROOT=$PWD
aws s3 cp s3://ferlab-public-dataset/nextflow/cnv-post-processing/V2/data-test data-test --recursive

# Uses the real truvari biocontainer; confirms bin/proband_priority.py works with its pysam
nf-test test modules/local/truvari_collapse --profile +docker
nf-test test subworkflows/local/utils_nfcore_cnv_post_processing_pipeline --profile +docker
```

You don't need any of the arm64 workarounds on amd64: no `~/.nextflow/config` override, no `DOCKER_DEFAULT_PLATFORM`, no custom truvari image.

## Known issues / follow-ups

1. **The S3 test dataset is incomplete.** `s3://.../cnv-post-processing/V2/data-test` has no `samplesheet.csv` and no `1463_NA12878_NA12891_NA12892.cnv.vcf.gz` (the `fam1_joint` VCF). Because of that, `tests/default.nf.test` (the pipeline test) fails at input validation. This predates this branch, and CI will hit it too.
   - Local workaround: copy `assets/samplesheet.csv` to `data-test/`, rewrite its paths to absolute, and drop the `fam1_joint` row. Then run `nextflow run . -profile test,docker --input data-test/samplesheet_nojoint.csv --outdir results`.
   - Real fix: re-upload the missing files to S3.
2. **DEL/DUP collapse into each other (not fixed here).** DRAGEN writes `SVTYPE=CNV` for both, and truvari reads both as the same unknown type, so a proband's DUP can absorb a parent's DEL at the same locus. A likely fix is to set `SVTYPE` from `ALT` before the collapse. This is documented under "Not covered" in the journal entry.
3. **Delete `HANDOFF.md`** before merging. The PR is squash-merged with a `<type>: <TICKET-123> <description>` title.
