# Ferlab-Ste-Justine/cnv-post-processing: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction


The Ferlab-Ste-Justine/cnv-post-processing pipeline is designed for the post-processing of Copy Number Variants (CNVs). Currently, it integrates [Exomiser](https://exomiser.readthedocs.io/en/14.0.0/running.html) for variant prioritization. We might support additional tools and steps in the future.

Exomiser requires clinical information for each sample, including details such as clinical signs, affected status, and gender. This information is provided in a separate file using the Phenopacket format. More details can be found in the samplesheet section.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet file format

Here is an example samplesheet file:

```csv title="samplesheet.csv"
sample,sequencingType,vcf,pheno
NA07019,WES,data-test/vcf/NA07019.chr22.cnv.vcf.gz,data-test/pheno/NA07019.yml
NA07022,WES,data-test/vcf/NA07022.chr22.cnv.vcf.gz,data-test/pheno/NA07022.yml
NA07056,WGS,data-test/vcf/NA07056.chr22.cnv.vcf.gz,data-test/pheno/NA07056.yml
```

| Column           | Description                                                                                       |
| ---------------- | ------------------------------------------------------------------------------------------------- |
| `sample`         | Unique identifier for the sample being analyzed.                                                  |
| `sequencingType` | Use `WES` for whole exome sequencing and `WGS` for whole genome sequencing                        |
| `vcf`            | Path to the vcf file containing the CNV data. The file must have a `.vcf` or `.vcf.gz` extension. |
| `pheno`          | Path to the phenotype file in Phenopacket format. Both `.yaml` and `.json` formats are supported. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Phenopacket file example

Here is an example phenopacket file, in yaml format:

```yaml title="NA07019_pheno.yml"
---
id: NA07019_FAMILY
proband:
  subject:
    id: NA07019
    sex: FEMALE
  phenotypicFeatures:
    - type:
        id: HP:0001159
        label: Syndactyly
```

Notes about the phenopacket file (`pheno` column):
- It should include only sample-specific information. The pedigree section must either contain only the sample or be omitted entirely.
- The sample identifier in the phenopacket file must match the sample identifier in the corresponding VCF file.

## Running the pipeline

Here is an example command to run the pipeline locally using the test profile with Docker:
```bash
nextflow run . -profile test,docker --outdir results
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

Instead using profiles, you can also pass pipeline settings via `.config` files and parameter files.
Here is an example command to run the pipeline in a production environment with a specific configuration and parameters:
```bash
nextflow -c application.config run Ferlab-Ste-Justine/cnv-post-processing \
    -r v1.0.0 \
    --input samplesheet.csv \
    --outdir results \
    -params-file params.json
```

Refer to the [nextflow documentation](https://www.nextflow.io/docs/latest/config.html#configuration-file) for details on the supported syntax for .config files.

For the [-params-file](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters) option, both `json` and `yaml` are supported. You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Ferlab-Ste-Justine/cnv-post-processing
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Ferlab-Ste-Justine/cnv-post-processing releases page](https://github.com/Ferlab-Ste-Justine/cnv-post-processing/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

*NOTE: This section has been auto-generated by nf-core. The pipeline may not be compatible with all the profiles described in this section.*

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
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
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
