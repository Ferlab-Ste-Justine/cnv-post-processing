# Ferlab-Ste-Justine/cnv-post-processing: Output

## Introduction

This document describes the output produced by the pipeline.

By default, the directories listed below will be created in the output directory specified by the `outdir` parameter after the pipeline has finished. All paths are relative to the root of this directory.

For certain steps, it is possible to specify a completely separate output directory using dedicated parameters. Unless explicitly stated otherwise, this document assumes that the default results directory, as specified by the `outdir` parameter, is being used.

Only the files listed below are published by default. Intermediate, per-sample files from the normalization and per-sample collapsing steps (before family-level merging) are not published unless `--publish_all` is set, in which case they land under a subfolder named after the process (e.g. `normalize_cnv/`, `truvari_collapse_persample/`).

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Family merge and collapse](#family-merge-and-collapse) - Merge a family's per-sample CNV calls and collapse redundant calls
- [Depth-based genotype refinement](#depth-based-genotype-refinement) - Resolve missing genotypes at collapsed sites using alignment depth
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
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: the merged family VCF (or, for a true standalone sample, its own collapsed VCF) after cohort-level [truvari](https://github.com/ACEnglish/truvari) collapse of redundant/overlapping CNV calls. Each collapsed record keeps one original call's `ID`, but every sample's `FORMAT/OID` names that sample's own original call (`:` percent-encoded as `%3A`, e.g. `DRAGEN%3ALOSS%3Achr22%3A17288730-17294633`). To get the proband's own call ID for a record, read the proband sample's `OID` and decode it; it's `.` when the proband has no call there. The record's `ID`, `POS` and `END` stay those of the call truvari kept, which may be a parent's. Every sample with a call in the cluster keeps its own `FORMAT` (`GT`, `CN`, `OID`...), including calls truvari itself would have left out (DRAGEN's `./1` gains and haploid `1` calls, put back from truvari's removed VCF); `INFO/NumConsolidated` counts them. The later family-level VCFs (`depth_refined/`, `ensemblvep/`, `slivar/`) keep it.
- `joint/`
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: for `caller=DRAGEN_JOINT` families, the family's already jointly-called VCF after lightweight prep only (drop non-variant records, split multiallelic sites) -- no merge or collapse needed.

</details>

### Depth-based genotype refinement

<details markdown="1">
<summary>Output files</summary>

- `depth_refined/`
  - `<familyId>.refined.vcf.gz`, `<familyId>.refined.vcf.gz.tbi`: the final family-level VCF that annotation reads from. Missing (`./.`) genotypes at collapsed sites are resolved to an explicit hom-ref call where [mosdepth](https://github.com/brentp/mosdepth)-derived depth confirms it, and left missing otherwise. A resolved genotype also gets `FORMAT/CN` set to the sample's expected copy number, the one slivar's rules use: 2 on autosomes and PAR; on non-PAR chrX 1 for a male (GT `0`) and 2 for a female; on non-PAR chrY 1 for a male (GT `0`). The depth band is scaled by that expected CN, so a normal male's chrX counts as normal. The sex comes from the family's `familyPed`; without one, non-PAR chrX/chrY genotypes stay missing, as do chrM and records crossing a PAR boundary. Records shorter than 5 kb are never refined: at that size, depth can't reliably tell a normal sample from a one-copy loss the caller missed, and a parent wrongly read as normal turns the child's inherited CNV into `denovo`. A sample with no call of its own on a record under 5 kb therefore stays `./.` with no CN, and slivar gives the child `possible_*` or `unknown_cn` there rather than `denovo`. Only produced this way for families where every sample has a `cram` in the samplesheet; other families get an unrefined pass-through of the `truvari/`/`joint/` VCF under this same path instead.

</details>

### VEP annotation

<details markdown="1">
<summary>Output files</summary>

- `ensemblvep/`
  - `variants.<familyId>.cnv.vep.vcf.gz`, `variants.<familyId>.cnv.vep.vcf.gz.tbi`: family-level [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)-annotated CNV VCF.
  - `variants.<familyId>.<sample>.cnv.vep.vcf.gz`, `variants.<familyId>.<sample>.cnv.vep.vcf.gz.tbi`: for standalone samples, the same annotation run per sample instead of per family.
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
  - `single_sample/<sample>.exomiser.{genes.tsv,html,json,variants.tsv,vcf.gz,vcf.gz.tbi}`: same six reports, for standalone samples run in single-sample mode instead of family mode.

</details>

The Exomiser analysis generates six reports, stored in the exomiser subfolder for each family (or, under `single_sample/`, each standalone sample). Each report uses the family or sample identifier as a filename prefix for easy identification.

In family mode, Exomiser runs on the proband VCF of the child the `familyPheno` phenopacket is about (its `proband.subject.id`), `slivar/<familyId>.<proband>.cnv.proband.vcf.gz`: the proband's own calls and intervals, with every family member's genotypes for the pedigree analysis. Its reports are named `<familyId>.<proband>.exomiser.*`. A family with no such proband VCF (a `DRAGEN_JOINT` family, no `familyPed`, or a phenopacket proband who isn't one of the PED's sequenced affected children) runs on the family VCF instead, as `<familyId>.exomiser.*`: the depth-refined one, or the VEP-annotated one under `--exomiser_start_from_vep`. Exomiser's scores stay in its own reports; they aren't written into the proband VCF.

A family gets one family-mode Exomiser run: a phenopacket names a single proband (`proband` isn't repeated in the Phenopacket schema). Other affected children (affected sibs, twins) are still used, through the phenopacket's pedigree, by Exomiser's segregation filters, and they still get their own proband VCF from `JOIN_PROBAND_CNV`, but no Exomiser report of their own. The `pheno` column is ignored on family rows.

For more details about these reports, refer to the [Exomiser documentation](https://exomiser.readthedocs.io/en/14.0.0/result_interpretation.html).

By default, exomiser output is saved in the `exomiser` subfolder within the main output directory. To save the output to a different location, use the `exomiser_outdir` parameter. In this case, the exomiser files will be written at the root of the specified location (with `single_sample/` still nested underneath it).

### Slivar

<details markdown="1">
<summary>Output files</summary>

- `slivar/`
  - `<familyId>.cnv.slivar.vcf.gz`, `<familyId>.cnv.slivar.vcf.gz.tbi`: the family-level, VEP-annotated CNV VCF with mode-of-inheritance tags added to `INFO` by [slivar](https://github.com/brentp/slivar) (mode of inheritance, from each affected child's call, so a family with several affected children can carry several: `denovo`, `denovo_suspect`, `inherited_maternal`, `inherited_paternal`, `recessive`, `recessive_lowconf`, `possible_denovo`, `possible_inherited_maternal`, `possible_inherited_paternal`, `inheritance_unresolved`, `mendelian_inconsistent`, `unknown_cn`, and on non-PAR chrX `x_linked_denovo` (a son's or a daughter's de novo CNV) and `x_linked_recessive` (a daughter with both X copies hit), listing the family's affected members; and parent of origin, per child, affected or not: `po_mother`, `po_father`, `po_both`, `po_ambiguous`, `po_denovo` (both parents known and normal), `po_possible_denovo` (one parent normal, the other unknown), `po_possible_mother` / `po_possible_father` (one parent unknown; the other carries a CNV it may have passed), named as in the SNV pipeline, listing the children it applies to). The mode of inheritance looks at each child and its parents only: whether the transmitting parent is affected, or what sibs carry, is left to the reader. Records on chrM, or crossing a chrX/chrY pseudoautosomal-region boundary, get no tag. What each tag means, with worked examples, is in `assets/moi_scenarios.md` and at the top of `assets/cnv-slivar-functions.js`. Indexed here since nothing downstream does it otherwise. Only produced for families with a `familyPed` in the samplesheet.
  - `<familyId>.<proband>.cnv.proband.vcf.gz`, `<familyId>.<proband>.cnv.proband.vcf.gz.tbi`: one per proband, i.e. each affected child in the PED (phenotype `2`, with a father or mother named) who was sequenced, `<proband>` being its VCF sample name. For each record of the slivar VCF above where the proband has a call, the proband's own call: its `ID`, `POS`, `REF`, `ALT`, `QUAL`, `FILTER` and `INFO` (`END`, `SVLEN`, the caller's own fields) come from the proband's input VCF (the samplesheet's `vcf`), matched through its decoded `FORMAT/OID`. The slivar record adds the `INFO` keys the call doesn't have (the slivar tags, VEP's `CSQ`, truvari's `NumCollapsed`, ...) and `FORMAT` with every family sample; `COHORT_ID` (percent-encoded like `OID`), `COHORT_POS` and `COHORT_END` name that slivar record. Records where the proband has no call are left out. VEP's `CSQ` and the slivar tags describe the truvari cluster the proband's call belongs to, computed on the interval of the call truvari kept to represent it, which may be a parent's. The slivar tags come from each sample's copy number, so they apply to the proband's call as is. `CSQ`, however, lists the genes the kept call overlaps: a gene overlapped only by the part of the proband's own interval outside the kept call isn't in it. Not produced for joint-called (`caller=DRAGEN_JOINT`) families, which have no per-sample VCF.

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
