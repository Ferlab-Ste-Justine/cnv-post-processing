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
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: the merged family VCF (or, for a true standalone sample, its own collapsed VCF) after cohort-level [truvari](https://github.com/ACEnglish/truvari) collapse of redundant/overlapping CNV calls.
- `joint/`
  - `<familyId>.sorted.vcf.gz`, `<familyId>.sorted.vcf.gz.tbi`: for `caller=DRAGEN_JOINT` families, the family's already jointly-called VCF after lightweight prep only (drop non-variant records, split multiallelic sites) -- no merge or collapse needed.

</details>

### Depth-based genotype refinement

<details markdown="1">
<summary>Output files</summary>

- `depth_refined/`
  - `<familyId>.refined.vcf.gz`, `<familyId>.refined.vcf.gz.tbi`: the final family-level VCF that annotation reads from. Missing (`./.`) genotypes at collapsed sites are resolved to an explicit hom-ref call where [mosdepth](https://github.com/brentp/mosdepth)-derived depth confirms it, and left missing otherwise. Only produced this way for families where every sample has a `cram` in the samplesheet; other families get an unrefined pass-through of the `truvari/`/`joint/` VCF under this same path instead.

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
  - `<familyId>.exomiser.genes.tsv`
  - `<familyId>.exomiser.html`
  - `<familyId>.exomiser.json`
  - `<familyId>.exomiser.variants.tsv`
  - `<familyId>.exomiser.vcf.gz`
  - `<familyId>.exomiser.vcf.gz.tbi`
  - `single_sample/<sample>.exomiser.{genes.tsv,html,json,variants.tsv,vcf.gz,vcf.gz.tbi}`: same six reports, for standalone samples run in single-sample mode instead of family mode.

</details>

The Exomiser analysis generates six reports, stored in the exomiser subfolder for each family (or, under `single_sample/`, each standalone sample). Each report uses the family or sample identifier as a filename prefix for easy identification.

For more details about these reports, refer to the [Exomiser documentation](https://exomiser.readthedocs.io/en/14.0.0/result_interpretation.html).

By default, exomiser output is saved in the `exomiser` subfolder within the main output directory. To save the output to a different location, use the `exomiser_outdir` parameter. In this case, the exomiser files will be written at the root of the specified location (with `single_sample/` still nested underneath it).

### Slivar

<details markdown="1">
<summary>Output files</summary>

- `slivar/`
  - `<familyId>.cnv.slivar.vcf.gz`, `<familyId>.cnv.slivar.vcf.gz.tbi`: the family-level, VEP-annotated CNV VCF with mode-of-inheritance tags added to `INFO` by [slivar](https://github.com/brentp/slivar) (`de_novo_candidate`, `dominant_inherited`, `recessive_candidate`, `candidate`, `ambiguous`, `unknown_cn`, `parent_of_origin_maternal`, `parent_of_origin_paternal`, `parent_of_origin_ambiguous`). This is the pipeline's last file, indexed here since nothing downstream does it otherwise. Only produced for families with a `familyPed` in the samplesheet.

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
