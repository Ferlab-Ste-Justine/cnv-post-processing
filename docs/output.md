# Ferlab-Ste-Justine/cnv-post-processing: Output

## Introduction

This document describes the output produced by the pipeline. What each step does to the VCFs, route by route, is in [routes.md](routes.md).

By default, the directories listed below will be created in the output directory specified by the `outdir` parameter after the pipeline has finished. All paths are relative to the root of this directory.

For certain steps, it is possible to specify a completely separate output directory using dedicated parameters. Unless explicitly stated otherwise, this document assumes that the default results directory, as specified by the `outdir` parameter, is being used.

Only the files listed below are published by default. Intermediate, per-sample files from the normalization and per-sample collapsing steps (before family-level merging) are not published unless `--publish_all` is set, in which case they land under a subfolder named after the process's first word (e.g. `normalize/` for `NORMALIZE_CNV`, `join/` for `JOIN_PROBAND_CNV`).

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Family merge and collapse](#family-merge-and-collapse) - Merge a family's per-sample CNV calls and collapse redundant calls
- [Genotype fill](#genotype-fill) - Fill the missing genotypes the family merge leaves, from each sample's own DRAGEN segments (or alignment depth)
- [Carrier filter](#carrier-filter) - Set each CNV's `FILTER` from its carriers' filters, and keep the CNVs that pass
- [VEP annotation](#vep-annotation) - Annotate CNVs with Ensembl VEP
- [Exomiser](#exomiser) - Prioritize variants with Exomiser
- [Slivar](#slivar) - Classify mode of inheritance
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Family merge and collapse

<details markdown="1">
<summary>Output files</summary>

- `merged/`
  - `<familyId>.merged.vcf.gz`, `<familyId>.merged.vcf.gz.tbi`: a family's per-sample CNV VCFs merged into one multi-sample VCF (`bcftools merge`). Not produced for families supplied as an already jointly-called (`caller=DRAGEN_JOINT`) VCF.
- `truvari/`
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: the merged family VCF (or, for a true standalone sample, its own collapsed VCF) after cohort-level [truvari](https://github.com/ACEnglish/truvari) collapse of redundant/overlapping CNV calls. Each collapsed record keeps one original call's `ID`, but every sample's `FORMAT/OID` names that sample's own original call (`:` percent-encoded as `%3A`, e.g. `DRAGEN%3ALOSS%3Achr22%3A17288730-17294633`). To get the proband's own call ID for a record, read the proband sample's `OID` and decode it; it's `.` when the proband has no call there. The record's `ID`, `POS` and `END` stay those of the call truvari kept, which may be a parent's. Every sample with a call in the cluster keeps its own `FORMAT` (`GT`, `CN`, `OID`, and `FT`, its call's `FILTER`...), including calls truvari itself would have left out (DRAGEN's `./1` gains and haploid `1` calls, put back from truvari's removed VCF); `INFO/NumConsolidated` counts them. The later family-level VCFs (`cn_filled/`, `ensemblvep/`, `slivar/`) keep it.
- `joint/`
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: for `caller=DRAGEN_JOINT` families, the family's already jointly-called VCF after lightweight prep only (drop non-variant records, split multiallelic sites) -- no merge or collapse needed.

</details>

### Genotype fill

<details markdown="1">
<summary>Output files</summary>

- `cn_filled/`
  - `<familyId>.cn_filled.vcf.gz`, `<familyId>.cn_filled.vcf.gz.tbi`: the final family-level VCF. In a family with a proband, `JOIN_PROBAND_CNV` joins its `FORMAT` onto the proband's calls, and this is the only published file with the CNVs only the parents carry; in a family without one, VEP annotates it. A sample with no call in a record (`./.`, no `CN`, left by the merge) gets the copy number of its own DRAGEN segments, read from its original per-sample VCF (REF blocks and its own calls): the CN covering at least half of the record. When that CN is the sample's expected one, the one slivar's rules use (2 on autosomes and PAR; on non-PAR chrX 1 for a male and 2 for a female; on non-PAR chrY 1 for a male; the sex comes from the family's `familyPed`), `GT` becomes hom-ref (`0/0`, or `0` for a male's non-PAR chrX/chrY). Any other CN, usually the sample's own call that the collapse left in another record, is written as `FORMAT/CN` only, with `GT` left `./.`; so is a CN with no expected one (chrM, a record crossing a PAR boundary, non-PAR chrX/chrY without a known sex). Records shorter than 5 kb are never filled: DRAGEN's per-sample caller misses many small CNVs, and a parent wrongly filled as normal turns the child's inherited CNV into `denovo`. A sample with no call of its own on a record under 5 kb therefore stays `./.` with no CN, and slivar gives the child `possible_*` or `unknown_cn` there rather than `denovo`. Written before the [carrier filter](#carrier-filter): it keeps every CNV, with the caller's `FILTER`. Not produced for joint-called (`caller=DRAGEN_JOINT`) families, whose VCF gives every sample a CN.
  - `<familyId>.refined.vcf.gz`, `<familyId>.refined.vcf.gz.tbi`: the same, under `--get_cn_by_computing_ratio`, filled from read depth instead ([mosdepth](https://github.com/brentp/mosdepth)): a missing genotype becomes hom-ref with the expected CN where the depth ratio over the record reads normal, and stays missing otherwise. Only for families where every sample has a `cram` in the samplesheet; other families aren't filled, and have no file here.

</details>

### Carrier filter

Before VEP, slivar and Exomiser, each CNV's `FILTER` is set from its carriers' `FORMAT/FT`, each sample's own filters from the caller (`FILTER_CNV_CARRIERS`, `bin/filter_cnv_carriers.py`). A carrier is a sample whose genotype has the CNV (`0/1`, `1/1`, `./1`, haploid `1`). `FILTER` becomes:

- `PASS` when at least one carrier's `FT` is `PASS`;
- otherwise the filters its carriers' `FT` name, e.g. `cnvLength;cnvQual`;
- `CarrierFT` when they name none.

With `--cnv_pass_only` (the default), the CNVs that don't pass are then dropped. With `--cnv_pass_only false`, they stay, flagged in `FILTER`; Exomiser still skips them, since it only reads `PASS` variants. DRAGEN flags every call under 10 kb `cnvLength` by default, so the default drops nearly all of them.

In a jointly-called DRAGEN VCF, `FT` is DRAGEN's own, for each sample's segment, and the `FILTER` it replaces is `PASS` as soon as any sample passes, even one without the CNV (`SampleFT` otherwise). On the per-sample route, `NORMALIZE_CNV` copies each call's `FILTER` into `FORMAT/FT`, which the merge and the collapses carry into the family VCF. The step runs on every family VCF, on each proband VCF after `JOIN_PROBAND_CNV` (a proband's call then passes when a parent's call in the same family record does), and on each standalone sample's VCF. `merged/`, `truvari/`, `joint/` and `cn_filled/` are written before it. Its own output is only published under `--publish_all`, in `filter/`.

### VEP annotation

<details markdown="1">
<summary>Output files</summary>

- `ensemblvep/`
  - `variants.<familyId>.cnv.vep.vcf.gz`, `variants.<familyId>.cnv.vep.vcf.gz.tbi`: family-level [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)-annotated CNV VCF, for a family without a proband or a joint-called (`caller=DRAGEN_JOINT`) family, after the [carrier filter](#carrier-filter).
  - `variants.<familyId>.<sample>.cnv.vep.vcf.gz`, `variants.<familyId>.<sample>.cnv.vep.vcf.gz.tbi`: the same annotation run on one sample's own calls: a proband's, all of them (before `JOIN_PROBAND_CNV` adds the family's `FORMAT` and the [carrier filter](#carrier-filter), see [Slivar](#slivar)), or a standalone sample's, after the carrier filter.
- `cache/` (only when `--download_cache` is set)
  - The downloaded VEP cache, named `<vep_cache_version>_<vep_genome>/`.

</details>

By default, VEP output is saved in the `ensemblvep` subfolder within the main output directory. To save it to a different location, use the `vep_outdir` parameter; the downloaded cache's location can be overridden separately with `outdir_cache`.

### Exomiser

<details markdown="1">
<summary>Output files</summary>

- `exomiser/`
  - `<familyId>.<proband>.exomiser.genes.tsv`
  - `<familyId>.<proband>.exomiser.html`
  - `<familyId>.<proband>.exomiser.json`
  - `<familyId>.<proband>.exomiser.variants.tsv`
  - `<familyId>.<proband>.exomiser.vcf.gz`
  - `<familyId>.<proband>.exomiser.vcf.gz.tbi`
  - `<familyId>.exomiser.{genes.tsv,html,json,variants.tsv,vcf.gz,vcf.gz.tbi}`: same six reports, for a family run on its family VCF (see below).
  - `single_sample/<familyId>.<sample>.exomiser.{genes.tsv,html,json,variants.tsv,vcf.gz,vcf.gz.tbi}`: same six reports, for standalone samples run in single-sample mode instead of family mode.

</details>

The Exomiser analysis generates six reports, stored in the exomiser subfolder for each family (or, under `single_sample/`, each standalone sample). Each report uses the family or sample identifier as a filename prefix for easy identification.

In family mode, Exomiser runs on the slivar VCF of the proband the `familyPheno` phenopacket is about (its `proband.subject.id`), `slivar/<familyId>.<proband>.cnv.slivar.vcf.gz`: the proband's own calls and intervals, with every family member's genotypes for the pedigree analysis. Its reports are named `<familyId>.<proband>.exomiser.*`. Any other family (no proband, a `DRAGEN_JOINT` family, no `familyPed`, or a phenopacket about a sample who isn't one of the family's probands) runs on its family VCF instead, the genotype-filled or joint one after the [carrier filter](#carrier-filter), as `<familyId>.exomiser.*`. Exomiser's scores stay in its own reports; they aren't written into the slivar VCF.

A family gets one family-mode Exomiser run: a phenopacket names a single proband (`proband` isn't repeated in the Phenopacket schema). Other affected children (affected sibs, twins) are still used, through the phenopacket's pedigree, by Exomiser's segregation filters, and they still get their own slivar VCF, but no Exomiser report of their own. The `pheno` column is ignored on family rows.

For more details about these reports, refer to the [Exomiser documentation](https://exomiser.readthedocs.io/en/14.0.0/result_interpretation.html).

By default, exomiser output is saved in the `exomiser` subfolder within the main output directory. To save the output to a different location, use the `exomiser_outdir` parameter. In this case, the exomiser files will be written at the root of the specified location (with `single_sample/` still nested underneath it).

### Slivar

<details markdown="1">
<summary>Output files</summary>

- `slivar/`
  - `<familyId>.<proband>.cnv.slivar.vcf.gz`, `<familyId>.<proband>.cnv.slivar.vcf.gz.tbi`: one per proband, i.e. each sequenced affected child of the PED (phenotype `2`, with a father or mother named), `<proband>` being its VCF sample name: the family's final VCF, built around the proband. Each record is one of the proband's own calls, from its VEP-annotated per-sample VCF: `ID`, `POS`, `REF`, `ALT`, `QUAL` and `INFO` (`SVTYPE` as `DEL`/`DUP`, `END`, the caller's own fields, and VEP's `CSQ` for the call's own interval). `FORMAT`, with every family sample, comes from the family record whose `FORMAT/OID` names the call (`JOIN_PROBAND_CNV`), and `FILTER` from that `FORMAT`'s carriers (see [Carrier filter](#carrier-filter): with `--cnv_pass_only`, only the calls that pass are kept): `COHORT_ID` (percent-encoded like `OID`), `COHORT_POS` and `COHORT_END` name that record, and `COHORT_NumCollapsed` / `COHORT_NumConsolidated` give truvari's counts for it (`NumCollapsed` / `NumConsolidated` are the proband's own, from the per-sample collapse). slivar then adds the same tags as in the family VCF below. A family with a proband gets no `<familyId>.cnv.slivar.vcf.gz`: the CNVs only its parents carry are in `cn_filled/`. Not produced for joint-called (`caller=DRAGEN_JOINT`) families, which have no per-sample VCF. The join's output before slivar is only published under `--publish_all`, in `join/`.
  - `<familyId>.cnv.slivar.vcf.gz`, `<familyId>.cnv.slivar.vcf.gz.tbi`: for a family without a proband or a joint-called family, the family-level, VEP-annotated CNV VCF with mode-of-inheritance tags added to `INFO` by [slivar](https://github.com/brentp/slivar) (mode of inheritance, from each affected child's call, so a family with several affected children can carry several: `denovo`, `denovo_suspect`, `inherited_maternal`, `inherited_paternal`, `recessive`, `recessive_lowconf`, `possible_denovo`, `possible_inherited_maternal`, `possible_inherited_paternal`, `inheritance_unresolved`, `mendelian_inconsistent`, `unknown_cn`, and on non-PAR chrX `x_linked_denovo` (a son's or a daughter's de novo CNV) and `x_linked_recessive` (a daughter with both X copies hit), listing the family's affected members; and parent of origin, per child, affected or not: `po_mother`, `po_father`, `po_both`, `po_ambiguous`, `po_denovo` (both parents known and normal), `po_possible_denovo` (one parent normal, the other unknown), `po_possible_mother` / `po_possible_father` (one parent unknown; the other carries a CNV it may have passed), named as in the SNV pipeline, listing the children it applies to). The mode of inheritance looks at each child and its parents only: whether the transmitting parent is affected, or what sibs carry, is left to the reader. Records on chrM, or crossing a chrX/chrY pseudoautosomal-region boundary, get no tag. What each tag means, with worked examples, is in `assets/moi_scenarios.md` and at the top of `assets/cnv-slivar-functions.js`. Indexed here since nothing downstream does it otherwise. Only produced for families with a `familyPed` in the samplesheet.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt` and `pipeline_dag_*.html`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
  - Software versions collated across every tool the pipeline ran: `cnv-post-processing_software_versions.yml`.
  - Copy of the nextflow log file: `nextflow.log`
  - Configuration files passed at the command line: `configs/<CONFIG_FILE_NAME>`
  - Metadata useful for reproducibility and traceability: `metadata.txt`

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
