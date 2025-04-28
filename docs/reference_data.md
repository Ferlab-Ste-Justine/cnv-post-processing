# Ferlab-Ste-Justine/cnv-post-processing: Reference Data


Reference files must be correctly downloaded and specified through pipeline parameters. This document provides a comprehensive list of the required reference files and explains how to set the pipeline parameters appropriately.


## Exomiser reference data

The `exomiser_data_dir` parameter specifies the path to the directory containing the exomiser reference files. This directory will be passed to the exomiser tool via the exomiser option `--exomiser.data-directory`.

It's content should look like this:
```
2402_hg19/
2402_hg38/
2402_phenotype/
remm/
  ReMM.v0.3.1.post1.hg38.tsv.gz
  ReMM.v0.3.1.post1.hg38.tsv.gz.tbi
cadd/1.7/
  gnomad.genomes.r4.0.indel.tsv.gz
  gnomad.genomes.r4.0.indel.tsv.gz.tbi
  whole_genome_SNVs.tsv.gz
  whole_genome_SNVs.tsv.gz.tbi
```

- *2402_hg19/* and *2402_hg38/*: These folders contain data associated with the `hg19` and `hg38` genome assemblies, respectively. The number `2402` corresponds to the exomiser data version.
- *remm/*: This folder is required only if REMM is used as a pathogenicity source in the exomiser analysis. In this case, additional parameters must be provided to specify the REMM data version (here `0.3.1.post1`) and the name of the .tsv.gz file to be used within this folder. See below.
- *cadd/*: This folder is required only if CADD is used as a pathogenicity source in the exomiser analysis. Here `1.7` is the CADD data version. As for REMM, additionnal parameters must be provided. See below.

To prepare the exomiser data directory, follow the instructions in the [exomiser installation documentation](https://exomiser.readthedocs.io/en/14.0.0/installation.html#linux-install)

Exomiser allows the use of a custom file for frequency data sources, typically to reduce the priority of high-frequency variants caused by artifacts. To use this feature, specify the `LOCAL` frequency source in the exomiser analysis file. Then, provide the paths to your custom frequency file and its index using the parameters `exomiser_local_frequency_path` and `exomiser_local_frequency_index_path`. Note that the index file is required if using this feature.

Together with the `exomiser_data_dir` parameter, these parameters must be provided to exomiser and should match the reference data available
- `exomiser_genome`: The genome assembly version to be used by exomiser. Accepted values are `hg38` or `hg19`.
- `exomiser_data_version`: The exomiser data version. Example: `2402`.
- `exomiser_cadd_version`: The version of the CADD data to be used by exomiser (optional). Example: `1.7`.
- `exomiser_cadd_indel_filename`: The filename of the exomiser CADD indel data file (optional). Example: `gnomad.genomes.r4.0.indel.tsv.gz`
- `exomiser_cadd_snv_filename`: The filename of the exomiser CADD snv data file (optional). Example: `whole_genome_SNVs.tsv.gz`
- `exomiser_remm_version`: The version of the REMM data to be used by exomiser (optional). Example:`0.3.1.post1`
- `exomiser_remm_filename`: The filename of the exomiser REMM data file (optional). Example: `ReMM.v0.3.1.post1.hg38.tsv.gz`
- `exomiser_local_frequency_path`: Path to a custom frequency source file (optional).
- `exomiser_local_frequency_index_path`: Path to the index file (.tbi) of the custom frequency source file (optional).

## Exomiser analysis files
In addition to the reference data, exomiser requires an analysis file (.yml/.json) that contains, among others things, the variant frequency sources for prioritization of rare variants, variant pathogenicity sources to consider, the list of filters and prioritizers to apply, etc.

Typically, different analysis settings are used for whole exome sequencing (WES) and whole genome sequencing (WGS) data.
Defaults analysis files are provided for each sequencing type in the assets folder:
- assets/exomiser/default_exomiser_WES_analysis.yml
- assets/exomiser/default_exomiser_WGS_analysis.yml

You can override these defaults and provide your own analysis file(s) via parameters `exomiser_analyis_wes` and `exomiser_analysis_wgs`.
Note that the default analysis files do not include REMM or CADD pathogenicity sources.

The exomiser analysis file format follows  the `phenopacket` standard and is described in detail [here](https://exomiser.readthedocs.io/en/latest/advanced_analysis.html#analysis).
There are typically multiple sections in the analysis file. To be compatible with the way we run the exomiser command, your analysis file should contain only the `analysis` section.

## Reference data parameters summary

| Parameter name | Required? | Description |
| --- | --- | --- |
| `exomiser_data_dir` | _Required_ | Path to the exomiser reference data directory |
| `exomiser_genome` | _Required_ | Genome assembly version to be used by exomiser(`hg19` or `hg38`) |
| `exomiser_data_version` | _Required_ | Exomiser data version (e.g., `2402`) |
| `exomiser_cadd_version` | _Optional_ | Version of the CADD data to be used by exomiser (e.g., `1.7`) |
| `exomiser_cadd_indel_filename`|	_Optional_ | Filename of the exomiser CADD indel data file (e.g., `gnomad.genomes.r4.0.indel.tsv.gz`) |
| `exomiser_cadd_snv_filename`|	_Optional_ | Filename of the exomiser CADD snv data file (e.g., `whole_genome_SNVs.tsv.gz`) |
| `exomiser_remm_version` | _Optional_ | Version of the REMM data to be used by exomiser (e.g., `0.3.1.post1`)|
| `exomiser_remm_filename` | _Optional_	| Filename of the exomiser REMM data file (e.g., `ReMM.v0.3.1.post1.hg38.tsv.gz`) |
| `exomiser_local_frequency_path`| _Optional_ | Path to a custom frequency source file |
| `exomiser_local_frequency_index_path`| _Optional_ | Path to the index file (.tbi) of the custom frequency source file. Required if specifying `exomiser_local_frequency_path`. |
| `exomiser_analysis_wes` | _Optional_ | Path to the exomiser analysis file for WES data, if different from the default |
| `exomiser_analysis_wgs` | _Optional_ | Path to the exomiser analysis file for WGS data, if different from the default |

