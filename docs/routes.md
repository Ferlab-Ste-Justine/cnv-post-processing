# Ferlab-Ste-Justine/cnv-post-processing: Routes

This document describes, step by step, what the pipeline does to each kind of input: the per-sample route (each sample called on its own by DRAGEN, then merged and collapsed into a family), the joint route (a family called jointly by DRAGEN), and the standalone route. For each step, it says what changes in the VCF and why. [usage.md](usage.md) covers the samplesheet, [output.md](output.md) the published files, and [schema.txt](schema.txt) draws the flow. The reasons behind each choice, with the measurements that led to them, are in [journal.md](journal.md).

Terms used throughout:

- **Call**: one CNV the caller reported for one sample (one record of a per-sample VCF).
- **Record**: one line of a VCF. In a family VCF, one record holds one column per sample.
- **Carrier**: a sample whose genotype (`GT`) has an ALT allele at a record: `0/1`, `1/1`, `./1` (a DRAGEN gain, zygosity unknown), or a haploid `1` (e.g. a male's chrX).
- **Segment** (joint VCFs only): a stretch of one sample's genome with a constant copy number.

## Which route a row takes

The route is decided per row, from the samplesheet only:

| Samplesheet                                                               | Route                                                 |
| ------------------------------------------------------------------------- | ----------------------------------------------------- |
| `caller = DRAGEN_JOINT` (one row per family)                              | [Joint route](#joint-route-caller--dragen_joint)      |
| `caller = DRAGEN`, a family of two rows or more                           | [Per-sample route](#per-sample-route-caller--dragen)  |
| `caller = DRAGEN`, one row with a `familyPed` or a `familyPheno`          | Per-sample route, without the merge (a family of one) |
| `caller = DRAGEN`, one row with neither (no other row has its `familyId`) | [Standalone route](#standalone-route)                 |

On the per-sample route, a family with a **proband** gets its final VCF built around that proband's own calls. A proband is a row whose VCF sample is an affected child in the family's PED: phenotype `2`, with a father or mother named (`isPedProband()`). A family can have several. The joint route never has a proband: there is no per-sample VCF to build its calls from.

## Steps at a glance

| Step                                   | Per-sample route                                       | Joint route                            | Standalone route                   |
| -------------------------------------- | ------------------------------------------------------ | -------------------------------------- | ---------------------------------- |
| Normalize                              | `NORMALIZE_CNV`, each sample                           | `NORMALIZE_CNV_JOINT`, the family VCF  | `NORMALIZE_CNV`                    |
| `INFO/SVTYPE` from the ALT             | yes                                                    | yes                                    | yes                                |
| `FORMAT/OID` (each call's own ID)      | added                                                  | no (no collapse)                       | added                              |
| `FORMAT/FT` (each sample's filters)    | copied from each call's `FILTER`                       | DRAGEN's own                           | copied                             |
| Collapse within a sample               | `TRUVARI_COLLAPSE_PERSAMPLE`                           | no                                     | yes                                |
| Merge into a family VCF                | `BCFTOOLS_MERGE`                                       | no (already one VCF)                   | no                                 |
| Collapse across the family             | `TRUVARI_COLLAPSE_COHORT` + `RESTORE_COLLAPSED_FORMAT` | no                                     | no                                 |
| Fill missing genotypes                 | from each sample's own DRAGEN segments                 | no (every sample is already genotyped) | no                                 |
| Carrier filter (`FILTER_CNV_CARRIERS`) | on the family VCF, and on each proband VCF             | on the family VCF                      | on its VCF                         |
| VEP                                    | each proband's own calls, or the family VCF            | the family VCF                         | its VCF                            |
| slivar (needs `familyPed`)             | each proband VCF, or the family VCF                    | the family VCF                         | no                                 |
| Exomiser                               | family mode (needs `familyPheno`)                      | family mode (needs `familyPheno`)      | single-sample mode (needs `pheno`) |

## The input VCFs

The two kinds of DRAGEN CNV VCF are built differently, which is why the routes differ.

**A per-sample DRAGEN CNV VCF** has one sample. Each record is either a call (`ALT` `<DEL>` or `<DUP>`) or a reference segment (`ALT` `.`, normal copy number). `FILTER` holds the call's own filters (`PASS`, or e.g. `cnvLength;cnvQual`), `QUAL` its quality, and `INFO/SVTYPE` is `CNV` for losses and gains alike. `FORMAT` is `GT:SM:CN:BC:PE:LR`. A loss is `0/1` (CN 1) or `1/1` (CN 0); a gain is `./1`, since DRAGEN doesn't tell which haplotype carries the extra copies.

```
#CHROM POS      ID                                  REF ALT   QUAL FILTER INFO                                         FORMAT            NA12878
chr22  23996246 DRAGEN:GAIN:chr22:23996247-24007899 N   <DUP> 37   PASS   SVLEN=11653;SVTYPE=CNV;END=24007899;REFLEN=... GT:SM:CN:BC:PE:LR ./1:1.58407:3:7:2,3:11.3
```

**A joint DRAGEN CNV VCF** (`--pedigree-file`, one `--cnv-input` per sample) has every family member. Each sample is cut into segments of constant copy number, and a new record starts wherever any sample's segment changes. So:

- a record has one `POS`/`END` for all samples, and each sample's `GT` and `CN` are its state over those coordinates;
- each sample's other values (`SM`, `BC`, `QS`, `LR` and its filters, `FT`) are those of its whole segment, which can span several records. A sample without the CNV usually sits in a long normal segment, repeated on every record it spans;
- a sample without the CNV has a `CN` (usually 2) even though its `GT` is `./.`;
- `FORMAT` is `GT:SM:CN:BC:PE:QS:FT:DN:LR`: `QS` replaces `QUAL` (which is `.`), `FT` holds each sample's own filters, and `DN` is DRAGEN's de novo flag (not used by the pipeline);
- `FILTER` is only `PASS` (at least one sample's `FT` is `PASS`) or `SampleFT` (none is). The sample that passes may not carry the CNV;
- some sites are multiallelic (`<DEL>,<DUP>`: one sample lost a copy, another gained one), and `INFO/SVLEN` is declared `Number=.`.

Four consecutive records of the NA12878 trio's joint VCF, where the mother (NA12892) has no CNV:

| Record (chr1)         | Length   | NA12878 (child)                 | NA12891 (father)            | NA12892 (mother)          |
| --------------------- | -------- | ------------------------------- | --------------------------- | ------------------------- |
| 13,279,537-15,823,963 | 2.54 Mb  | CN 2, BC 2147, `PASS`           | CN 2, BC 2147, `PASS`       | CN 2, **BC 2711, `PASS`** |
| 15,823,963-15,825,729 | 1,766 bp | CN 2, BC 1, `cnvLength;cnvQual` | **CN 1**, BC 1, `cnvLength` | CN 2, **BC 2711, `PASS`** |
| 15,825,729-15,828,479 | 2,750 bp | **CN 1**, BC 2, `cnvLength`     | **CN 0**, BC 2, `cnvLength` | CN 2, **BC 2711, `PASS`** |
| 15,828,480-16,526,823 | 698 kb   | CN 2, BC 561, `PASS`            | CN 2, BC 561, `PASS`        | CN 2, **BC 2711, `PASS`** |

The mother's column is the same on all four records: her single 3.25 Mb normal segment (2,711 bins of about 1 kb), which passes. That `PASS` makes the third record's `FILTER` `PASS`, although both carriers' calls fail `cnvLength`. This is why the pipeline judges a record by its carriers' `FT`, not by `FILTER` (see [The carrier filter](#the-carrier-filter)).

## Per-sample route (caller = DRAGEN)

### 1. `NORMALIZE_CNV`, each sample

One pipe of bcftools and awk steps, in this order:

1. **Drop the reference segments** (`bcftools view -e 'ALT="."'`): they aren't CNVs.
2. **Sort, and split multiallelic sites** (`bcftools norm -m -any`). DRAGEN's per-sample VCFs have none; the step is shared with the joint route, whose VCFs do.
3. **Set `INFO/SVTYPE` from the ALT** (`bin/svtype_from_alt.awk`): `<DEL>` gives `DEL`, `<DUP>` gives `DUP`, in place of DRAGEN's `CNV`. truvari reads the variant type from `SVTYPE`; with `CNV` on every call, it could collapse a loss with a gain, leaving the absorbed sample's `GT` pointing at the wrong ALT. Other ALTs are left alone.
4. **Copy the call's `ID` into `FORMAT/OID`** (`bin/add_original_id.awk`, `ext.original_id`). The collapses later keep one `ID` per cluster of calls, often another sample's; `OID` keeps naming each sample's own call. `:` separates `FORMAT` values and DRAGEN IDs contain it, so it's percent-encoded as `%3A` (and `%` as `%25`).
5. **Copy the call's `FILTER` into `FORMAT/FT`** (`bin/add_sample_filter.awk`, `ext.sample_filter`): the field a joint VCF already gives each sample. Merged into a family VCF, each sample keeps its own call's filters, which `FILTER` can't hold.

The same call before and after:

```
before: ... <DUP> 37 PASS SVLEN=11653;SVTYPE=CNV;END=24007899;... GT:SM:CN:BC:PE:LR        ./1:1.58407:3:7:2,3:11.3
after:  ... <DUP> 37 PASS SVLEN=11653;SVTYPE=DUP;END=24007899;... GT:SM:CN:BC:PE:LR:OID:FT ./1:1.58407:3:7:2,3:11.3:DRAGEN%3AGAIN%3Achr22%3A23996247-24007899:PASS
```

### 2. `TRUVARI_COLLAPSE_PERSAMPLE`, each sample

Removes redundant calls within one sample: [truvari](https://github.com/ACEnglish/truvari) `collapse --pctseq 0.00 --pctsize 0.70 --pctovl 0.70 --sizemax -1`, i.e. calls of the same type that overlap each other by 70% and whose sizes are within 70% of each other. `--sizemax -1` lifts truvari's default 50 kb cap, which would leave the largest CNVs out. The kept call gets `INFO/NumCollapsed` / `NumConsolidated`.

### 3. Probands, and their own VEP run

A proband's collapsed VCF is annotated by VEP on its own (`VCF_ANNOTATE_ENSEMBLVEP_PROBAND`), all its calls, filtered or not. Each call keeps its own interval, `SVTYPE` (`DEL`/`DUP`) and `ID`, and gets a `CSQ` for its own interval. Step 8 joins the family's genotypes onto it.

### 4. `BCFTOOLS_MERGE`, the family

`bcftools merge` of the family's per-sample VCFs, in sample-name order, then sort and split multiallelic sites. Only identical calls (same `POS`, `END` and ALT) become one record; any other overlap stays as separate records for the next step. A sample without the call gets an empty column (`./.:.:.:...`: no `GT`, no `CN`, no `FT`). The merged record's `FILTER` combines the joined calls' filters.

A family of one (a single row with a `familyPed` or `familyPheno`) skips the merge.

### 5. `TRUVARI_COLLAPSE_COHORT`, then `RESTORE_COLLAPSED_FORMAT`

truvari `collapse --refdist 2000 --sizemin 1000 --pctovl 0.3 --pctsize 0.7 --pctseq 0 --sizemax -1` clusters the calls of different samples that describe the same CNV: same type, within 2 kb of each other, overlapping by 30% and with sizes within 70%. Calls under 1 kb aren't collapsed. truvari keeps the leftmost call of a cluster (the longest on a tie), with its `ID`, `POS`, `END`, `FILTER` and `INFO`, and copies each absorbed call's `FORMAT` into its sample's empty column ("consolidation"). The `OID` copied with it still names that sample's own call.

`RESTORE_COLLAPSED_FORMAT` (`bin/restore_collapsed_format.py`) fixes two ways the collapse loses a sample's call:

- truvari only copies an absorbed call's `FORMAT` when it reads its `GT` as set, and it reads DRAGEN's `./1` (a gain) and haploid `1` as unset. Those columns are copied back from truvari's removed-calls VCF (`-c`), and `NumConsolidated` counts them.
- When two calls of one sample fold into the same kept record, the record's single column for that sample holds one of them. The other gets a record of its own, as truvari writes a call it leaves alone, found by comparing `OID`s. With these thresholds that's rare (one call at exactly 70% of the kept call's length, the other starting where it ends).

### 6. `FILL_CN_FROM_SEGMENTS`, the missing genotypes

The merge leaves a sample without a call empty (`./.`, no `CN`), and slivar needs every parent's `CN`. DRAGEN's per-sample VCF has the answer for most of them: it covers the genome, with REF blocks (`ALT` `.`, `CN` 2, or 1 on a male's non-PAR chrX/chrY) between the sample's calls. Step 1 dropped them, so this step reads each sample's original per-sample VCF again, the samplesheet's `vcf`.

For each empty genotype on a record of 5 kb or more (`ext.min_size`), it counts the bases of the record each CN covers in the sample's REF blocks and own calls, whatever their `FILTER`. When one CN covers at least half of the record (`--cn_from_segment_overlap`, 0.5 by default; a tie stays empty), it writes that `CN`:

- the sample's expected CN, the one slivar's rules use (2 on autosomes and PAR; on non-PAR chrX 1 for a male and 2 for a female; on non-PAR chrY 1 for a male; the sex comes from the family's PED): `GT` becomes hom-ref, `0/0`, or a haploid `0` on a male's non-PAR chrX/chrY;
- any other CN: `CN` only, `GT` stays `./.`. It's usually the sample's own call, which the merge and collapse left in another record. slivar reads `CN` only; `GT` stays empty because no call of that sample has this record's ALT, so the carrier filter and Exomiser don't count it as a carrier.

With no expected CN (chrM, a record crossing a PAR boundary, non-PAR chrX/chrY without a known sex), the `CN` is written and `GT` stays `./.`. Only empty genotypes are touched: an existing call is never changed, and no CNV is called. The other `FORMAT` values of a filled sample stay missing, its `FT` included. Records under 5 kb stay empty: DRAGEN's per-sample caller misses many small CNVs, and a parent filled as normal there turns the child's inherited CNV into `denovo`.

Two family records after steps 4 to 6, from the test data (a collapsed `FORMAT`, `GT:CN:FT` only). The mother's 55 kb loss and the child's and father's 28 kb homozygous loss stayed apart in the collapse (they overlap by 51%, but the shorter is 51% of the longer's size, under `--pctsize 0.7`); each record's empty genotypes get the sample's own call's CN from the other:

```
#CHROM POS      ID                                  ALT   NA12878  NA12891  NA12892
chr22  23932072 DRAGEN:LOSS:chr22:23932073-23987577 <DEL> ./.:0:.  ./.:0:.  0/1:1:cnvLikelihoodRatio
chr22  23932072 DRAGEN:LOSS:chr22:23932073-23960275 <DEL> 1/1:0:PASS 1/1:0:PASS ./.:1:.
```

The depth-based fill had put the mother at `0/0` and `CN` 2 on the second record: her depth there reads normal for a loss she carries on a longer interval. slivar now sees one deleted copy passed by each parent.

**From depth instead (`--get_cn_by_computing_ratio`).** `BAM_VCF_DEPTH_GENOTYPE_REFINEMENT`, the fill before this one, kept as an option. It runs only when every sample of the family has a CRAM; otherwise the family VCF goes on unchanged, with a warning:

1. `VCF_TO_BED`: the family VCF's records as regions.
2. `MOSDEPTH_RATIO`, per sample: each region's mean depth divided by the sample's genome-wide mean depth.
3. `DEPTH_GENOTYPE_REFINE` (`bin/refine_genotypes.py`): for each empty genotype, the ratio is divided by (expected CN / 2). Between 0.8 and 1.2, the genotype becomes hom-ref and `CN` the expected copy number; otherwise it stays empty. Same 5 kb floor, same expected CN.

It only ever writes a normal CN, and can't tell a sample's own call in a neighbouring record from a missed CNV. On 7 full-genome trios, against the DRAGEN joint caller, it got 357 of 568 parents' missing genotypes on records of 10 kb or more right, against 533 for the segments, and slivar's mode of inheritance agreed with the joint VCF's on 86.4% of shared calls, against 97.7% (see [journal.md](journal.md)).

The result is published in `cn_filled/` (`<familyId>.cn_filled.vcf.gz`, or `<familyId>.refined.vcf.gz` from depth). It keeps every CNV, filtered or not.

### 7. `FILTER_CNV_CARRIERS`, the family VCF

Sets each record's `FILTER` from its carriers' `FT` and, under `--cnv_pass_only` (the default), drops the records no carrier passes; see [The carrier filter](#the-carrier-filter). Its output goes to family-level VEP (a family without a proband) and to family-mode Exomiser when it runs on the family VCF (step 10). The join (step 8) reads the unfiltered family VCF.

### 8. A family with a proband: `JOIN_PROBAND_CNV`, then the carrier filter

`JOIN_PROBAND_CNV` (`bin/join_proband_cnv.py`) builds one VCF per proband, from its own VEP-annotated calls (step 3). For each call, it finds the family record whose proband column's `OID` (decoded) and ALT are that call's, and writes:

- from the proband's call: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL` and `INFO` (with its `CSQ`);
- from the family record: `FORMAT` and every sample's column;
- `INFO/COHORT_ID` (percent-encoded like `OID`), `COHORT_POS` and `COHORT_END`, naming that family record, and `COHORT_NumCollapsed` / `COHORT_NumConsolidated` when it has them.

A proband call found in no family record, or in two, stops the run: every call should be in exactly one (step 5). The CNVs only the parents carry aren't in the proband VCF; they stay in `cn_filled/`.

`FILTER_CNV_CARRIERS` then runs on each proband VCF. The proband's own `FILTER` gives way to the family record's carriers: a proband call that fails is kept when a parent's call in the same family record passes.

A family without a proband instead gets VEP on its carrier-filtered family VCF (`VCF_ANNOTATE_ENSEMBLVEP`).

### 9. slivar

`SLIVAR_EXPR` runs on each proband VCF, or on the VEP-annotated family VCF of a family without a proband, when the family has a `familyPed`. Before it, the run stops if the VCF declares no `FORMAT/CN` (the rules read copy number only), or if the PED doesn't match the VCF's sample names (`validatePedAgainstVcf`; a warning only under `--skip_ped_validation`).

slivar reads each sample's `CN` and the PED, and adds 14 mode-of-inheritance tags as `INFO` keys valued with the family's affected members (`denovo`, `inherited_maternal`, `unknown_cn`...), and 8 parent-of-origin tags valued with the children they apply to (`po_mother`, `po_denovo`...). It also adds its own `impactful`, `genic` and `highest_impact_order` from VEP's `CSQ`. The rules are in `assets/cnv-slivar-functions.js`, every named scenario's result in `assets/moi_scenarios.md`. The output, indexed by `BCFTOOLS_INDEX_SLIVAR`, is the family's final VCF: `slivar/<familyId>.<proband>.cnv.slivar.vcf.gz` per proband, or `slivar/<familyId>.cnv.slivar.vcf.gz`.

A proband VCF record, abridged (`INFO` without `CSQ`):

```
chr22 10711935 DRAGEN:GAIN:chr22:10711936-10742589 N <DUP> . PASS
  INFO   SVTYPE=DUP;END=10742589;...;COHORT_ID=DRAGEN%3AGAIN%3Achr22%3A10711936-10742589;COHORT_POS=10711935;COHORT_END=10742589;inheritance_unresolved=NA12878;impactful;genic;highest_impact_order=8
  FORMAT GT:SM:CN:BC:PE:LR:OID:FT
  NA12878 ./1:1.80207:4:14:2,4:42.2:DRAGEN%3AGAIN%3Achr22%3A10711936-10742589:PASS
  NA12891 ./1:2.32743:5:14:0,3:75.5:DRAGEN%3AGAIN%3Achr22%3A10711936-10742589:PASS
  NA12892 ./.:.:.:.:.:.:.:.
```

### 10. Exomiser, family mode

Runs when the family has a `familyPheno`. When the phenopacket's proband (`proband.subject.id`) is one of the family's probands, on that proband's slivar VCF: `exomiser/<familyId>.<proband>.exomiser.*`. Otherwise (no proband, or a phenopacket about another sample), on the carrier-filtered family VCF, before VEP: `exomiser/<familyId>.exomiser.*`. Exomiser's analysis mode is `PASS_ONLY`: it reads `FILTER` only, which is why `FILTER` is set from the carriers before it.

## Joint route (caller = DRAGEN_JOINT)

### 1. `NORMALIZE_CNV_JOINT`

The same process as `NORMALIZE_CNV`, with different options (`conf/modules.config`):

1. **Drop the reference segments** (`ALT` `.`).
2. **Sort, and split multiallelic sites**: a `<DEL>,<DUP>` site becomes a `<DEL>` record and a `<DUP>` record, each sample's `GT` recoded to that record's allele (10 sites in the NA12878 trio).
3. **Drop `INFO/SVLEN`** (`ext.strip_svlen`): DRAGEN declares it `Number=.` in joint VCFs, so the split can leave several values on one record, and Exomiser crashes reading it. VEP and Exomiser take the size from `POS`, `END` and `SVTYPE`.
4. **Set `INFO/SVTYPE` from the ALT**, as on the per-sample route.

There is no `OID` (`ext.original_id = false`): the joint VCF isn't collapsed, so each record's `ID` is already DRAGEN's. There is no `FILTER` to `FT` copy (`ext.sample_filter = false`): each sample's `FT` is DRAGEN's own, and the record's `FILTER` (`PASS` or `SampleFT`) would overwrite it. The result is published in `joint/`.

### 2. What the joint route skips, and why

- **The collapses and the merge**: the family is already one VCF, one record per segment boundary, called by one caller.
- **Genotype fill**: joint calling gives every sample a `CN` on every record, carrier or not. The fill exists for the empty columns a merge leaves, which a joint VCF doesn't have.
- **The proband route**: there's no per-sample VCF to take a proband's own calls from, so a joint family's final VCF is family-level (`slivar/<familyId>.cnv.slivar.vcf.gz`) and its Exomiser runs on the family VCF.

### 3. `FILTER_CNV_CARRIERS`, then VEP, slivar and Exomiser

The carrier filter runs on the joint VCF, as on a per-sample route's family VCF, with DRAGEN's own `FT`. Its output goes to family-level VEP, then slivar (needs `familyPed`), and to family-mode Exomiser (needs `familyPheno`) before VEP, as `exomiser/<familyId>.exomiser.*`.

## Standalone route

A sample with no family (a single row, no `familyPed`, no `familyPheno`):

1. `NORMALIZE_CNV` and `TRUVARI_COLLAPSE_PERSAMPLE`, as on the per-sample route (with `OID` and `FT`).
2. `FILTER_CNV_CARRIERS_SOLO`: the sample is the only carrier of each call, so each call keeps its own `FILTER`, and under `--cnv_pass_only` only the `PASS` calls go on.
3. VEP: `ensemblvep/variants.<familyId>.<sample>.cnv.vep.vcf.gz`.
4. Exomiser in single-sample mode, when the row has a `pheno`: `exomiser/single_sample/<familyId>.<sample>.exomiser.*`.

No merge, genotype fill or slivar: there is no family to compare with.

## The carrier filter

`FILTER_CNV_CARRIERS` (`bin/filter_cnv_carriers.py`) runs on every VCF before VEP, slivar and Exomiser, on all three routes. It replaces each record's `FILTER`:

- `PASS` when at least one carrier's `FT` is `PASS`;
- otherwise the filters the carriers' `FT` name, e.g. `cnvLength;cnvQual;dinucQual`;
- `CarrierFT` when no carrier's `FT` names one (no carrier, or an unset `FT`).

Under `--cnv_pass_only` (the default), it then drops every record whose new `FILTER` isn't `PASS`. With `--cnv_pass_only false`, they stay, flagged; Exomiser still skips them.

What `FT` is depends on the route: on the per-sample route, each call's own `FILTER`, copied by `NORMALIZE_CNV`; in a joint VCF, DRAGEN's filters for each sample's segment. Either way, the rule only looks at carriers, so a sample without the CNV can't make a record pass.

On the five chr1 records of the NA12878 trio's joint VCF used by the module test:

| Record    | DRAGEN's `FILTER` | Carriers' `FT`                                  | New `FILTER` | Kept by default |
| --------- | ----------------- | ----------------------------------------------- | ------------ | --------------- |
| 4063702   | `SampleFT`        | all three `cnvLength`                           | `cnvLength`  | no              |
| 13225212  | `PASS`            | all three `PASS`                                | `PASS`       | yes             |
| 15823963  | `PASS`            | NA12891 `cnvLength` (the mother passes, no CNV) | `cnvLength`  | no              |
| 15825729  | `PASS`            | NA12878, NA12891 `cnvLength` (the same)         | `cnvLength`  | no              |
| 248647941 | `PASS`            | NA12878 `cnvLength`, NA12891 `PASS`             | `PASS`       | yes             |

DRAGEN sets `cnvLength` on every call under 10 kb by default (`--cnv-filter-length`), so the default drops nearly all of them: about 4 records in 5 on either route, in the full-genome trios. "Matching" calls depend on the route: on the joint route, carriers on the same record; on the per-sample route, calls truvari put in the same cluster.

## Fields the pipeline adds or changes

| Field                                                                                         | Set by                              | Routes                               | Meaning                                                                                                  |
| --------------------------------------------------------------------------------------------- | ----------------------------------- | ------------------------------------ | -------------------------------------------------------------------------------------------------------- |
| `INFO/SVTYPE`                                                                                 | `NORMALIZE_CNV`                     | all                                  | `DEL` / `DUP` from the ALT, in place of `CNV`                                                            |
| `INFO/SVLEN`                                                                                  | `NORMALIZE_CNV_JOINT`               | joint                                | removed                                                                                                  |
| `FORMAT/OID`                                                                                  | `NORMALIZE_CNV`                     | per-sample, standalone               | the sample's own call `ID`, `:` as `%3A`; `.` where the sample has no call                               |
| `FORMAT/FT`                                                                                   | `NORMALIZE_CNV` (DRAGEN's on joint) | all                                  | the sample's own filters: its call's `FILTER`, or its segment's on joint                                 |
| `INFO/NumCollapsed`, `NumConsolidated`, `CollapseId`                                          | truvari, `RESTORE_COLLAPSED_FORMAT` | per-sample, standalone               | calls collapsed into the record, sample columns copied into it, the cluster's id                         |
| `GT`, `FORMAT/CN` of an empty sample                                                          | `FILL_CN_FROM_SEGMENTS`             | per-sample                           | the CN of the sample's own segments; `GT` hom-ref when it's the expected CN (records of 5 kb or more)    |
| `FILTER`                                                                                      | `FILTER_CNV_CARRIERS`               | all                                  | `PASS` when a carrier's `FT` is `PASS`, the carriers' filters otherwise, `CarrierFT` when they name none |
| `INFO/CSQ`                                                                                    | VEP                                 | all                                  | consequences, for the record's own interval                                                              |
| `INFO/COHORT_ID`, `COHORT_POS`, `COHORT_END`, `COHORT_NumCollapsed`, `COHORT_NumConsolidated` | `JOIN_PROBAND_CNV`                  | per-sample, proband VCFs             | the family record the proband's call was joined to                                                       |
| mode-of-inheritance and `po_*` tags; `impactful`, `genic`, `highest_impact_order`             | slivar                              | per-sample, joint (with `familyPed`) | see [output.md](output.md#slivar) and `assets/moi_scenarios.md`                                          |

## How the two routes compare

On 7 full-genome CEPH trios run both ways (2026-11-03, before the carrier filter existed, so on every call):

- **CNVs of 5 kb or more**: the two routes are equivalent. A CNV only one parent carries is in the child's calls 47-48% of the time on both (about 50% expected), and matched calls get the same mode-of-inheritance tag 84% of the time.
- **CNVs under 5 kb**: DRAGEN's single-sample calling misses about half of the small inherited CNVs that joint calling finds (24.5% against 52.7% of the private parental CNVs found in the child). The gap is in the calls DRAGEN makes, before the pipeline runs; merging and collapsing can't recover them. The per-sample route also leaves half the parents' genotypes under 5 kb without a `CN`, since the genotype fill doesn't touch them.
- **The depth-based fill** (now `--get_cn_by_computing_ratio`) could read a parent as normal where joint calling finds a loss, e.g. near the chr21 centromere (segmental duplications), turning an inherited CNV into `denovo`. Filling from the samples' own segments, the default since, removes those.

When DRAGEN joint calling is available for a family, it's the better input, mainly for small CNVs. The carrier filter and the genotype fill, with its 5 kb floor, are discussed, with their evidence, in [journal.md](journal.md).
