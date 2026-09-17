# pQTL downstream analysis pipeline

Snakemake pipeline for downstream analysis of pQTL summary statistics.

This workflow is designed to process pQTL summary-statistics files, identify significant loci, select representative variants, annotate loci with gene/protein information, add array-version and literature annotations, collapse related signals, detect genomic hotspots, apply heterogeneity filtering, and append single-study information.

The pipeline can be run locally or on an HPC cluster using the included SLURM profile.


## Pipeline overview

The workflow is controlled by:

```text
workflow/Snakefile
```

The main configuration file is:

```text
config/config.<project>.yaml
```

The Snakefile includes the following rule files:

```text
workflow/rules/common.smk
workflow/rules/locus_breaker.smk
workflow/rules/mr_instruments.smk
workflow/rules/annotation.smk
workflow/rules/heterogeneity.smk
```

The main workflow steps are:

1. Read the list of pQTL summary-statistics files.
2. Run LocusBreaker to identify independent loci.
3. Collect loci across all proteins/analytes.
4. Select best SNPs for downstream MR-style analyses.
5. Map loci to proteins and genes.
6. Annotate loci with gene/protein information.
7. Add array-version annotation.
8. Run backward-literature annotation.
9. Collapse related or duplicated signals.
10. Detect genomic hotspots and lonespots.
11. Apply heterogeneity filtering.
12. Append single-study results.


## Repository structure

```text
pqtl_downstream/
├── config/
│   ├── config_backward_literature_believe.json
│   ├── config_backward_literature_metaanalysis.json
│   ├── config_believe.yaml
|   ├── config_metaanalysis.yaml
|   ├── config_example.yaml.yaml
│   ├── believe_filtered_harmonized_sumstats.txt
│   ├── metaanalysis_filtered_harmonized_sumstats.txt
│   ├── local.txt
│   └── example.txt
├── data/
├── slurm/
│   └── config.yaml
├── workflow/
│   ├── Snakefile
│   ├── envs/
│   │   ├── locus_breaker.yml
│   │   ├── r_environment.yml
│   │   └── single_studies.yml
│   ├── rules/
│   │   ├── common.smk
│   │   ├── locus_breaker.smk
│   │   ├── mr_instruments.smk
│   │   ├── annotation.smk
│   │   ├── heterogeneity.smk
│   │   └── collapsing.smk
│   └── scripts/
│       ├── LB/
│       ├── MR/
│       ├── backward_literature/
│       ├── collapsing/
│       ├── gp_annotation/
│       ├── heterogeneity/
│       ├── hotspot/
│       ├── single_studies/
│       └── version_mapping/
├── Makefile
├── environment.yml
├── submit.sbatch
├── dag.svg
├── dag.pdf
└── README.md
```


## Installation

Clone the repository:

```bash
git clone https://github.com/ht-diva/pqtl_downstream.git
cd pqtl_downstream
```

Create or update the Snakemake environment:

```
make dependencies
```

For development dependencies:

```
make dev-dependencies
```



## Configuration

The main configuration file defines the input files, output folder, genome build, column names, thresholds, annotation files, and downstream analysis parameters.

Project configurations are stored as:

```
config/config.<project>.yaml
```

Existing examples include:

- `config/config.example.yaml`
- `config/config.believe.yaml`
- `config/config.metaanalysis.yaml`

---

### Selecting a project

Use one of the Makefile helpers:

```
make project-believe
make project-metaanalysis
make project-example
```

The selected project name is stored in `.project`. If `.project` is absent, the Makefile uses the default project defined by `DEFAULT_PROJECT`.

To remove the current selection:

```
make clean-project
```

You can also bypass project selection and call Snakemake with an explicit configuration:

```
snakemake \
    --profile slurm \
    --snakefile workflow/Snakefile \
    --configfile config/config.believe.yaml
```

---

### Run options

Supported modes for `input:`:

| Mode | Description |
|---|---|
| `run_LB` | Run LocusBreaker directly from summary-statistics files. |
| `LB` | Start from an existing LocusBreaker output file. |

When using `LB` mode, provide:

```yaml
LB_file: "/path/to/collected_loci_excluding_mhc.csv"
```

---

### Summary-statistics list

```yaml
sumstats_list: "/config/<project>_filtered_harmonized_sumstats.txt"
```

This file should contain one summary-statistics file path per line.

Example:

```text
/path/to/protein_1.sumstats.tsv.gz
/path/to/protein_2.sumstats.tsv.gz
/path/to/protein_3.sumstats.tsv.gz
```

The pipeline derives the analyte/protein identifier from each file path.

---

### Output directories

```yaml
workspace_path: "results"
destination_path: "dest"
```

Pipeline outputs are written inside:

```text
workspace_path
```

Logs are written inside:

```text
workspace_path/logs/
```

---

### Summary-statistics column names

The default configuration expects the following columns:

| labels      | believe     | metaanalysis |
| ----------- | ----------- | ------------ |
| `chr_label` | `"CHR"`     | `"'##CHR'"`  |
| `pos_label` | `"POS"`     | `"POS"`      |
| `p_label`   | `"MLOG10P"` | `"MLOG10P"`  |

Update these values if your summary-statistics files use different column names.

---

### LocusBreaker thresholds

| Parameter | Description | Values (believe) | Values (metaanalysis) |
|---|---|---|---|
| `p1` | Primary genome-wide significance threshold. | `1.178411501296253e-11` | `1.256913021618904e-11` |
| `p2` | Secondary threshold used during locus construction. | `1e-06` | `1e-06` |
| `hole` | Genomic distance used to separate or merge loci. | `3000000` | `3000000` |

---

### Locus-selection settings

| Parameter | Description | Values (believe) | Values (metaanalysis) |
|---|---|---|---|
| `NLP12` | Toggle for NLP12-specific locus handling. | `0` | `0` |
| `MHC` | Toggle for MHC-region handling. | `0` | `0` |
| `build` | Genome build used for genomic coordinates. | `38` | `37` |

---

### Mapping and annotation files

| File | Description | Values (believe) | Values (metaanalysis) |
|---|---|---|---|
| `mapping_filepath` | Protein/analyte mapping file. | `"/exchange/healthds/pQTL/BELIEVE/Cis_trans_mapping/Build38_mapping_file/results/somascan_tss_ncbi_grch38_version_20251210.txt"` | `/exchange/healthds/pQTL/Reference_datasets_for_QC_proteomics/Cis_trans_mapping/somascan_tss_ncbi_grch37_ensembl_version_20241216.txt` |
| `gtf_file` | GTF annotation file used for gene/protein annotation. | `/exchange/healthds/pQTL/BELIEVE/Cis_trans_mapping/Build38_mapping_file/GCF_000001405.40_GRCh38.p14_genomic.gtf` | `/exchange/healthds/public_data/reference_genomes/GRCh37/GCF_000001405.25_GRCh37.p13_genomic.gtf` |
| `array_list_path` | Directory containing array-version annotation information. | `data` | `data` |

---

### Backward-literature configuration

```yaml
BL_config_file: "config/config_backward_literature_<project>.json"
```

This JSON file controls the backward-literature annotation step.

---

### Heterogeneity parameters

| Parameter | Description | Values (believe) | Values (metaanalysis) |
|---|---|---|---|
| `nef` | Effective sample-size or analysis-specific parameter used by downstream scripts. | `4243` | `3978` |
| `Isquare` | I-squared threshold used for heterogeneity filtering. | `90` | `90` |

---

### Hotspot finder settings

| Parameter | Description | Values (believe) | Values (metaanalysis) |
|---|---|---|---|
| `hotspot_window_size` | Window size used to detect hotspots. | `5000000` | `5000000` |
| `chr_col` | Chromosome column name. | `chr` | `chr` |
| `start_col` | Locus start-position column name. | `start` | `start` |
| `end_col` | Locus end-position column name. | `end` | `end` |
| `hotspot_threshold` | Minimum number of loci required to define a hotspot. | `50` | `50` |
| `lonespot_window_size` | Window size used to detect lonespots. | `10000000` | `10000000` |
| `lonespot_threshold` | Maximum number of loci allowed for lonespot definition. | `1` | `1` |

---

### Single studies

```yaml
single_studies:
  - CHRIS
  - INTERVAL
```

The pipeline creates study-specific appended outputs for each study listed here.



## Running the pipeline

### Recommended checks

Confirm the selected project:

```
make project-<project>
```

Inspect the planned jobs without executing them:

```
make dry-run
```

Generate a DAG:

```
make dag
```

This writes dag.svg.

### Submit the workflow

Submit the supplied SLURM wrapper:

```
sbatch submit.sbatch
```

Alternatively, start it through the Makefile from an appropriate execution environment:

```
make run
```

### Resume an interrupted workflow

Snakemake normally resumes from the existing outputs. To explicitly rerun incomplete jobs:

```
make rerun
```

### Unlock after an interrupted run

If the controlling Snakemake process was killed, the working directory may remain locked. First confirm that no other Snakemake process is using the same directory, then run:

```
make unlock
```

Never unlock a directory while another workflow is actively writing the same outputs.