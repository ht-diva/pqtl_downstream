# pQTL downstream analysis


# pQTL downstream analysis pipeline

Snakemake pipeline for downstream analysis of pQTL summary statistics.

This workflow is designed to process pQTL summary-statistics files, identify significant loci, select representative variants, annotate loci with gene/protein information, add array-version and literature annotations, collapse related signals, detect genomic hotspots, apply heterogeneity filtering, and append single-study information.

The pipeline can be run locally or on an HPC cluster using the included SLURM profile.

---

## Pipeline overview

The workflow is controlled by:

```text
workflow/Snakefile
```

The main configuration file is:

```text
config/config.yaml
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

---

## Repository structure

```text
pqtl_downstream/
├── config/
│   ├── config.yaml
│   ├── config_backward_literature.json
│   ├── believe_filtered_harmonized_sumstats.txt
│   ├── qced_filtered_meta_path.txt
│   ├── local.txt
│   ├── test.txt
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

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ht-diva/pqtl_downstream.git
cd pqtl_downstream
git checkout believe
```

## Configuration

Before running the pipeline, edit:

```text
config/config.yaml
```

This file defines the input files, output folder, genome build, column names, thresholds, annotation files, and downstream analysis parameters.

---

## Main configuration options

### Input mode

```yaml
input: "run_LB"
```

Supported modes:

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
sumstats_list: "/path/to/believe_filtered_harmonized_sumstats.txt"
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
workspace_path: "results_BELIEVE_LB"
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

```yaml
labels:
  chr_label: "CHR"
  pos_label: "POS"
  p_label: "MLOG10P"
```

Update these values if your summary-statistics files use different column names.

For example, if your file uses `chromosome`, `position`, and `pvalue`, change the config to:

```yaml
labels:
  chr_label: "chromosome"
  pos_label: "position"
  p_label: "pvalue"
```

---

### LocusBreaker thresholds

```yaml
thresholds:
  p1: 1.178411501296253e-11
  p2: 1e-06
  hole: 3000000
```

| Parameter | Description |
|---|---|
| `p1` | Primary genome-wide significance threshold. |
| `p2` | Secondary threshold used during locus construction. |
| `hole` | Genomic distance used to separate or merge loci. |

---

### Locus-selection settings

```yaml
loci_selection:
  NLP12: 0
  MHC: 0
  build: 38
```

| Parameter | Description |
|---|---|
| `NLP12` | Toggle for NLP12-specific locus handling. |
| `MHC` | Toggle for MHC-region handling. |
| `build` | Genome build used for genomic coordinates. |

---

### Mapping and annotation files

```yaml
mapping_filepath: "/path/to/somascan_tss_ncbi_grch38_version_20250104.txt"
gtf_file: "/path/to/GCF_000001405.40_GRCh38.p14_genomic.gtf"
array_list_path: "data"
```

| File | Description |
|---|---|
| `mapping_filepath` | Protein/analyte mapping file. |
| `gtf_file` | GTF annotation file used for gene/protein annotation. |
| `array_list_path` | Directory containing array-version annotation information. |

---

### Backward-literature configuration

```yaml
BL_config_file: "config/config_backward_literature.json"
```

This JSON file controls the backward-literature annotation step.

---

### Heterogeneity parameters

```yaml
params:
  nef: 4243
  Isquare: 90
```

| Parameter | Description |
|---|---|
| `nef` | Effective sample-size or analysis-specific parameter used by downstream scripts. |
| `Isquare` | I-squared threshold used for heterogeneity filtering. |

---

### Hotspot finder settings

```yaml
hotspot_finder:
  hotspot_window_size: 5000000
  chr_col: "chr"
  start_col: "start"
  end_col: "end"
  hotspot_threshold: 50
  lonespot_window_size: 10000000
  lonespot_threshold: 1
```

| Parameter | Description |
|---|---|
| `hotspot_window_size` | Window size used to detect hotspots. |
| `chr_col` | Chromosome column name. |
| `start_col` | Locus start-position column name. |
| `end_col` | Locus end-position column name. |
| `hotspot_threshold` | Minimum number of loci required to define a hotspot. |
| `lonespot_window_size` | Window size used to detect lonespots. |
| `lonespot_threshold` | Maximum number of loci allowed for lonespot definition. |

---

### Single studies

```yaml
single_studies:
  - CHRIS
  - INTERVAL
```

The pipeline creates study-specific appended outputs for each study listed here.

---

## Running the pipeline

The repository includes a `Makefile` with common commands.

To see available commands:

```bash
make
```

---

### Dry run

Before launching the full workflow, run:

```bash
make dry-run
```

This checks which jobs Snakemake would run without actually executing them.

---

### Run locally

To run the workflow locally:

```bash
make local-run
```

This runs Snakemake with local cores and software deployment through Conda and Apptainer/Singularity.

---

### Run on SLURM

To run the workflow using the SLURM profile:

```bash
make run
```

This uses:

```text
slurm/config.yaml
```

as the Snakemake profile.

---

### Submit with sbatch

A wrapper SLURM submission script is provided:

```bash
sbatch submit.sbatch
```

## Notes

This README describes the BELIEVE branch of the pQTL downstream pipeline.

Before running the workflow on a new system, update all absolute paths in `config/config.yaml` and run a dry run:

```bash
make dry-run
```