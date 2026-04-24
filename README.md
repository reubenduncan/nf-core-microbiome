# core-microbiome

A Nextflow pipeline for core microbiome detection from microbiome feature tables, producing prevalence/detection heatmap data across a range of thresholds.

## Introduction

This pipeline accepts a feature table (BIOM, TSV, or GTDB format) and a sample metadata CSV, then identifies core taxa using the `microbiome` R package. It sweeps a range of detection thresholds (log-spaced) at a fixed prevalence cutoff and records which taxa meet both thresholds at each step. The output supports downstream heatmap visualisation of core microbiome membership across threshold space.

An optional `--merge_parquet` flag consolidates all outputs into a single Parquet file.

## Quick start

```bash
nextflow run main.nf \
  --feature_table /path/to/feature_table.biom \
  --meta_table    /path/to/meta_table.csv \
  --groups_column Treatment \
  --which_level   Genus \
  --label         my_analysis
```

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--feature_table` | *(required)* | Path to feature table (BIOM, TSV, or GTDB format) |
| `--meta_table` | *(required)* | Path to sample metadata CSV (first column = sample IDs) |
| `--taxonomy_table` | `""` | Taxonomy TSV (required for `tsv`/`gtdb` input formats) |
| `--input_format` | `biom` | `biom` \| `tsv` \| `gtdb` |
| `--output_dir` | `results/` | Directory for output files |
| `--scripts_dir` | `/opt/ecology-scripts` | Path to R scripts (override for local runs) |

### Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min_library_size` | `5000` | Minimum per-sample read depth; samples below this are dropped |
| `--exclude_column` | `""` | Metadata column used to identify samples for exclusion |
| `--exclude_values` | `""` | Comma-separated values in `exclude_column` to remove |

### Grouping

| Parameter | Default | Description |
|---|---|---|
| `--groups_column` | `""` | Metadata column for the primary grouping variable |
| `--groups_paste_columns` | `""` | Comma-separated columns pasted together to form groups |
| `--type_column` | `""` | Optional additional metadata column for subsetting |
| `--type2_column` | `""` | Optional second additional metadata column for subsetting |
| `--type2_levels` | `""` | Comma-separated levels of `type2_column` to retain |

### Core microbiome

| Parameter | Default | Description |
|---|---|---|
| `--which_level` | `Phylum` | Taxonomic level for analysis (`Kingdom` \| `Phylum` \| `Class` \| `Order` \| `Family` \| `Genus` \| `Otus`) |
| `--label` | `analysis` | Label prepended to output file names |
| `--prevalence_min` | `0.85` | Minimum prevalence (fraction of samples) to define a taxon as core |
| `--what_detection` | `absolute` | Detection threshold type: `absolute` (read counts) \| `relative` (proportions) |
| `--detection_steps` | `20` | Number of log-spaced detection threshold steps to sweep |
| `--detection_range_min` | `-1` | Minimum detection threshold (`-1` = auto: 1 for absolute, 0.001 for relative) |
| `--detection_range_max` | `-1` | Maximum detection threshold (`-1` = auto: max(abundances)/10 for absolute, 0.2 for relative) |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `--merge_parquet` | `false` | Merge all output CSVs into a single Parquet file |

## Outputs

All files are written to `--output_dir`.

| File | Description |
|---|---|
| `Core_heatmap_{label}_{level}_{detection}.csv` | Heatmap matrix of core membership across detection thresholds (rows = taxa, columns = thresholds) |
| `Core_members_{label}_{level}_{detection}.csv` | List of taxa meeting the core definition at `--prevalence_min` |
| `core_microbiome_{label}.parquet` | Both CSVs merged with `analysis` and `table` metadata columns (`--merge_parquet` only) |

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- Docker (default) **or** a local R installation with: `optparse`, `vegan`, `dplyr`, `tidyr`, `stringr`, `phyloseq`, `microbiome`, `arrow`

## Running without Docker

```bash
nextflow run main.nf \
  -c nextflow.config \
  --scripts_dir   "$(pwd)" \
  --feature_table /path/to/table.biom \
  --meta_table    /path/to/meta.csv \
  --groups_column Treatment \
  --which_level   Genus \
  --label         my_analysis
```
