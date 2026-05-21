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
