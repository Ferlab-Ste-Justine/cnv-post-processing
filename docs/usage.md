# Ferlab-Ste-Justine/cnv-post-processing: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The Ferlab-Ste-Justine/cnv-post-processing pipeline post-processes Copy Number Variant (CNV) calls, from per-sample VCFs to family-level, annotated and prioritized results: normalization and redundant-call collapsing ([truvari](https://github.com/ACEnglish/truvari)), family-level merging (or a family VCF already jointly called by DRAGEN), a fill of the missing genotypes the merge leaves, from each sample's own DRAGEN segments (or alignment depth, with [mosdepth](https://github.com/brentp/mosdepth)), a filter on the caller's own per-sample filters, annotation with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html), prioritization with [Exomiser](https://exomiser.readthedocs.io/en/14.0.0/running.html), and mode-of-inheritance classification with [slivar](https://github.com/brentp/slivar). A family with an affected child gets its final VCF built around that child's own calls. Standalone samples with no family or pedigree information get their own lightweight VEP + Exomiser route. The flow is drawn in [schema.txt](schema.txt), each route is described step by step in [routes.md](routes.md), and the outputs in [output.md](output.md).

Exomiser requires clinical information for each sample or family, including details such as clinical signs, affected status, and gender. This information is provided in a separate file using the Phenopacket format. slivar's mode-of-inheritance classification requires a pedigree (PED) file. More details can be found in the samplesheet section.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet file format

Here is an example samplesheet file, covering a trio (`fam1`), the same trio supplied as an already jointly-called VCF (`fam1_joint`), and a standalone sample (`ind1`):

```csv title="samplesheet.csv"
familyId,sample,sequencingType,caller,vcf,cram,pheno,familyPheno,familyPed
fam1,NA12878,WGS,DRAGEN,data-test/vcf/NA12878.chr22.cnv.vcf.gz,data-test/cram/NA12878.chr22.cram,,data-test/pheno/trio_NA12878_91_92.yaml,data-test/ped/trio_NA12878_91_92.ped
fam1,NA12891,WGS,DRAGEN,data-test/vcf/NA12891.chr22.cnv.vcf.gz,data-test/cram/NA12891.chr22.cram,,data-test/pheno/trio_NA12878_91_92.yaml,data-test/ped/trio_NA12878_91_92.ped
fam1,NA12892,WGS,DRAGEN,data-test/vcf/NA12892.chr22.cnv.vcf.gz,data-test/cram/NA12892.chr22.cram,,data-test/pheno/trio_NA12878_91_92.yaml,data-test/ped/trio_NA12878_91_92.ped
fam1_joint,NA12878_91_92,WGS,DRAGEN_JOINT,data-test/vcf/1463_NA12878_NA12891_NA12892.cnv.vcf.gz,,,data-test/pheno/trio_NA12878_91_92.yaml,data-test/ped/trio_NA12878_91_92.ped
ind1,NA12878_ind,WGS,DRAGEN,data-test/vcf/NA12878.chr22.cnv.vcf.gz,data-test/cram/NA12878.chr22.cram,data-test/pheno/NA12878.yaml,,
```

| Column           | Required? | Description                                                                                                                                                                                                                                                                                                                                        |
| ---------------- | --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `familyId`       | Required  | Identifier grouping the samples that make up a family. A standalone sample with no relatives still needs a `familyId` (its own, used by no other row).                                                                                                                                                                                             |
| `sample`         | Required  | Identifier for the sample, unique within its family. It names the sample's own output files (`<familyId>.<sample>`), and can differ from the sample name in the VCF; the PED, the phenopackets and a proband's output files use the VCF's name.                                                                                                    |
| `sequencingType` | Required  | `WES` (whole exome) or `WGS` (whole genome). Selects Exomiser's analysis file (`--exomiser_analysis_wes` or `--exomiser_analysis_wgs`).                                                                                                                                                                                                            |
| `caller`         | Required  | The tool that produced the CNV VCF: `DRAGEN` for a sample's own DRAGEN CNV VCF, `DRAGEN_JOINT` for a DRAGEN jointly-called, multi-sample CNV VCF of the whole family. A `DRAGEN_JOINT` family has exactly one row, and no per-sample rows. `GATK` is accepted by the samplesheet but not supported yet: the run stops at normalization.            |
| `vcf`            | Required  | Path to the CNV VCF (`.vcf` or `.vcf.gz`). slivar reads copy number from its `FORMAT/CN`.                                                                                                                                                                                                                                                          |
| `cram`           | Optional  | Path to the sample's alignment file (`.cram` or `.bam`), with its index next to it (`<file>.crai` or `<file>.bai`). Only used with `--get_cn_by_computing_ratio`, to fill a family's missing genotypes from depth; a family is only filled that way when every one of its samples has one. Not used for `DRAGEN_JOINT` rows or standalone samples. |
| `pheno`          | Optional  | Path to a per-sample phenopacket (`.yaml`/`.yml`/`.json`). Enables Exomiser in single-sample mode, for a standalone sample only: it's ignored on the rows of a family, and not allowed on a `DRAGEN_JOINT` row.                                                                                                                                    |
| `familyPheno`    | Optional  | Path to a family phenopacket, identical on every row of the family (all rows or none). Enables Exomiser in family mode, see [Phenopacket files](#phenopacket-files-pheno-familypheno).                                                                                                                                                             |
| `familyPed`      | Optional  | Path to a PED file (`.ped`), identical on every row of the family (all rows or none). Enables slivar's mode-of-inheritance classification (GRCh38 only) and decides the family's probands, see [PED file](#ped-file-familyped).                                                                                                                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

A sample with neither `familyPed` nor `familyPheno` and no other family members (`familyId` used by no other row) is treated as a true standalone sample: it's routed through its own per-sample VEP + Exomiser (single-sample mode) steps instead of the family-level route, since the family-level route (merge, cohort collapse, genotype fill, family-level VEP/Exomiser, slivar) would have nothing to add for it. A single sample with a `familyPed` or a `familyPheno` takes the family route instead, without the merge.

The pipeline checks the samplesheet before any task runs, and stops if:

- the rows of a family don't all have the same `familyPed`, or the same `familyPheno` (an empty value counts as a value);
- a `DRAGEN_JOINT` family has more than one row, or a `DRAGEN_JOINT` row has a `pheno`;
- a family has a `familyPed` and `--vep_genome` isn't `GRCh38`: slivar's chrX/chrY model uses the GRCh38 pseudoautosomal regions.

The caller's filters decide which CNVs go on to VEP, slivar and Exomiser: by default (`--cnv_pass_only`), only the CNVs at least one carrier passes. DRAGEN flags every call under 10 kb `cnvLength`, so nearly all of them are left out. Set `--cnv_pass_only false` to keep them, flagged in `FILTER` (see [output.md](output.md#carrier-filter)).

### Phenopacket files (`pheno`, `familyPheno`)

Exomiser reads the patient's phenotypes (HPO terms) from a [Phenopacket](https://phenopacket-schema.readthedocs.io/) file, in YAML or JSON. Its sample identifiers must be the VCFs' sample names, not the samplesheet's `sample` column.

- `familyPheno`: a family phenopacket, with the proband's phenotypes (`proband`) and a `pedigree` listing the family members in the VCF. Its `proband.subject.id` decides which VCF Exomiser runs on: that proband's slivar VCF when it's one of the family's probands (see [PED file](#ped-file-familyped)), the family VCF otherwise. A family gets one Exomiser run, about that proband: other affected children only count through the pedigree. Exomiser reads this pedigree and slivar the PED, so keep their sexes and affected statuses consistent.
- `pheno`: a standalone sample's phenopacket, with no pedigree, or one listing only that sample.

Here is the test data's family phenopacket (`data-test/pheno/trio_NA12878_91_92.yaml`):

```yaml title="trio_NA12878_91_92.yaml"
---
id: family1
proband:
  subject:
    id: NA12878
    sex: FEMALE
  phenotypicFeatures:
    - type:
        id: HP:0001159
        label: Syndactyly

pedigree:
  persons:
    - individualId: NA12878
      paternalId: NA12891
      maternalId: NA12892
      sex: FEMALE
      affectedStatus: AFFECTED
    - individualId: NA12891
      sex: MALE
      affectedStatus: UNAFFECTED
    - individualId: NA12892
      sex: FEMALE
      affectedStatus: UNAFFECTED

metaData:
  resources:
    - id: hp
      name: human phenotype ontology
      url: http://purl.obolibrary.org/obo/hp.owl
      version: hp/releases/2019-11-08
      namespacePrefix: HP
      iriPrefix: "http://purl.obolibrary.org/obo/HP_"
  phenopacketSchemaVersion: 2.0
```

### PED file (`familyPed`)

A standard [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) file describing the family: whitespace-separated columns (family, sample, father, mother, sex `1`/`2`/`0`, phenotype `2` for affected), lines starting with `#` ignored. It's used by the slivar mode-of-inheritance classification step, and by the genotype fill for each sample's sex. It must be identical across every row of the same family, and its sample identifiers must be the VCFs' sample names (their `#CHROM` header line), not the samplesheet's `sample` column: for `DRAGEN` rows, the name in each per-sample VCF; for a `DRAGEN_JOINT` row, the names in the joint VCF. Its sex column must be correct: on non-PAR chrX and chrY, a normal male carries one copy, so slivar and the genotype fill use each sample's sex to decide what copy number is normal. A wrong sex misclassifies every call there, and a sex of `0` (unknown) makes slivar skip non-PAR chrX/chrY records for that family. Check sex against the data (e.g. chrX heterozygosity, or the CNV caller's own sex estimate) before running the pipeline.

The PED also decides the family's probands: its affected children (phenotype `2`, with a father or mother named) that have their own `DRAGEN` row. A family with a proband gets its final VCF built around that proband's own calls (`slivar/<familyId>.<proband>.cnv.slivar.vcf.gz`, see [output.md](output.md#slivar)) instead of a family-level one.

Before running slivar, the pipeline checks each family's PED against the sample names in the family VCF, and stops with an error naming the problem and its fix if:

- a sample in the VCF has no row in the PED. slivar would leave that sample out of the family, and a PED matching none of the VCF's names would give no tags at all;
- a sequenced sample's father or mother is named but isn't in the VCF. Write a parent who wasn't sequenced as `0`: when the named father is missing, slivar also drops the child's link to the mother, and a child with no parent in the family gets no mode-of-inheritance or parent-of-origin tag.

Rows for relatives who weren't sequenced can stay in the PED; they aren't checked. slivar exits normally in all these cases, so without the check the tags would be silently missing. To run anyway, for example for tests, set `--skip_ped_validation`: the same problems are then logged as warnings, and those families' slivar tags may be missing or wrong.

## Install reference data and prepare parameters file

To run the pipeline, you need to configure the required reference data and provide the necessary pipeline parameters

### Parameters File

Pipeline parameters can be passed using a JSON or YAML file, referred to as a parameters file. Using a parameters file is highly recommended for reproducibility and ease of configuration.

You can find more details about parameter files in the [nextflow documentation](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters). Additionally, the pipeline schema provides a complete list of available parameters and their descriptions.

### Reference data

The pipeline requires specific reference data to function correctly.
Detailed instructions on the reference data, its structure, and how to configure the associated parameters can be found in the [reference_data.md](reference_data.md) documentation.

## Running the pipeline

The pipeline needs Nextflow 24.10.5 or later, before 26.0.0 (`manifest.nextflowVersion`). Here is an example command to run the pipeline locally using the test profile with Docker, from a clone of the repository. The test data isn't in the repository: sync it to `data-test/` first.

```bash
aws s3 cp s3://ferlab-public-dataset/nextflow/cnv-post-processing/V2/data-test data-test --recursive
nextflow run . -profile test,docker --outdir results
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

Instead of using profiles, you can also pass pipeline settings via `.config` files and parameter files.
Here is an example command to run the pipeline in a production environment with a specific configuration and parameters:

```bash
nextflow -c application.config run Ferlab-Ste-Justine/cnv-post-processing \
    -r v2.0.0 \
    --input samplesheet.csv \
    --outdir results \
    -params-file params.json
```

Refer to the [nextflow documentation](https://www.nextflow.io/docs/latest/config.html#configuration-file) for details on the supported syntax for .config files.

For the [-params-file](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters) option, both `json` and `yaml` are supported. You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

By default, all pipeline outputs are saved in the directory specified by the `--outdir` parameter. However, you can customize the output locations for specific steps using other parameters (`--vep_outdir`, `--exomiser_outdir`), and `--publish_all` also publishes the intermediate files. For more details on that, see [output.md](output.md).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Ferlab-Ste-Justine/cnv-post-processing
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Ferlab-Ste-Justine/cnv-post-processing releases page](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/releases) and find the latest pipeline version (e.g. `v2.0.0`). Then specify this when running the pipeline with `-r` (one hyphen) - e.g. `-r v2.0.0`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

_NOTE: This section has been auto-generated by nf-core. The pipeline may not be compatible with all the profiles described in this section._

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> Use Docker or Singularity containers. Conda isn't supported: the pipeline's own modules have no Conda environment.

Unlike nf-core pipelines, this one doesn't load institutional profiles from [nf-core/configs](https://github.com/nf-core/configs) when it runs, so that no configuration is ever fetched from a remote source. Pass your institution's settings with `-c` instead.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Uses the test data in `data-test/` (see [Running the pipeline](#running-the-pipeline)), so needs no other parameters but `--outdir`
- `debug`
  - For troubleshooting: prints each task's host and dumps task hashes (to see why `-resume` re-runs a task)
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. If a task exits with a code that usually means it ran out of memory or time (130 to 145, or 104), it's resubmitted once with twice the resources; if it fails again, the pipeline stops (`conf/base.config`).

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website. The pipeline sets `docker.registry` (and Singularity's, Podman's...) to `quay.io`, so give a container's full address, e.g. `docker.io/ensemblorg/ensembl-vep:release_114.2`: a name without a registry is looked up on quay.io.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website. The defaults are in `conf/modules.config`. For example, the shortest record the genotype fill fills (5 kb), or, with `--get_cn_by_computing_ratio`, the depth fill's hom-ref depth-ratio band (0.8 to 1.2), can be changed in a config file passed with `-c` (the share of a record one copy number must cover in the sample's own segments is a parameter, `--cn_from_segment_overlap`, 0.5 by default):

```groovy
process {
    withName: FILL_CN_FROM_SEGMENTS {
        ext.min_size = 10000
    }
    withName: DEPTH_GENOTYPE_REFINE {
        ext.low      = 0.85
        ext.high     = 1.15
        ext.min_size = 10000
    }
}
```

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
