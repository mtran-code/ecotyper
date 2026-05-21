# EcoTyper Pipeline

This fork keeps the original EcoTyper R implementation but runs it from a reproducible Snakemake layout.

## Layout

- `config.yml`: run configuration.
- `Snakefile`: root workflow entrypoint.
- `workflow/Snakefile`: modular Snakemake rules.
- `workflow/scripts/`: R and Python scripts used by the workflow.
- `data/rawdata/`: local dataset path anchors, organized by dataset name. Do not commit datasets here.
- `data/procdata/`: generated staging files and rendered EcoTyper compatibility configs.
- `data/results/`: EcoTyper working files and final outputs.
- `docs/`: project documentation.

## Setup

Create the pixi environment and install the CRAN packages that are not available for `osx-arm64` from the configured conda channels:

```sh
pixi install
pixi run r-install-legacy
```

## Configure

Edit `config.yml`. For local datasets, point inputs at files under `data/rawdata/<dataset-name>/`.

For scRNA discovery, set `input.expression_format: h5ad` and `input.h5ad_file` to read sparse AnnData inputs directly. The workflow stages EcoTyper-compatible expression and annotation files under `data/procdata/h5ad/<run.name>/`, using `obs` columns configured by `input.h5ad_cell_type_column` and `input.h5ad_sample_column`.

Set `input.h5ad_max_cells_per_cell_type` to sample cells per cell type before matrix blocks are read from disk. The selected cells are recorded in `data/procdata/h5ad/<run.name>/selected_cells.tsv`.

For larger h5ad inputs, set `input.h5ad_direct_cell_state_inputs: true` to bypass the legacy global dense `data.txt` staging path. Direct mode reads the h5ad in backed mode, applies optional cell sampling and HVG selection first, then writes the per-cell-type EcoTyper intermediates consumed by state discovery. It supports external annotation joins through `input.h5ad_annotation_file`, `input.h5ad_annotation_id_column`, and `input.h5ad_obs_id_column`.

Direct h5ad mode does not support `pipeline.filter_genes: cell type specific`; use `pipeline.filter_genes: no filter` with `input.h5ad_feature_selection: hvg` and `input.h5ad_n_top_genes` instead.

The default run mode is `discovery_scRNA`. The base workflow also supports `discovery_bulk` and `discovery_presorted`.

## Run

Dry-run first:

```sh
pixi run dry-run
```

Run the workflow:

```sh
pixi run snakemake
```

Snakemake renders a legacy EcoTyper config under `data/procdata/config/` and writes final outputs under `data/results/<run.name>/`.

## NMF Backend

The default `pipeline.nmf_backend: r` uses the original R `NMF` implementation. Set `pipeline.nmf_backend: torch` or `auto` to run NMF updates through PyTorch; `auto` uses CUDA or Apple MPS when available and otherwise falls back to CPU. GPU-backed runs still write `estim.RData` restart files for downstream EcoTyper steps.
