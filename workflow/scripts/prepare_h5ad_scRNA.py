import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def make_safe_ids(values):
    safe_values = []
    seen = {}
    for value in values:
        safe = re.sub(r"[^0-9A-Za-z_.]", ".", str(value))
        if re.match(r"^[0-9]", safe):
            safe = f"X{safe}"
        if not safe:
            safe = "cell"

        count = seen.get(safe, 0)
        seen[safe] = count + 1
        safe_values.append(safe if count == 0 else f"{safe}.{count + 1}")
    return safe_values


def make_safe_labels(values):
    return [re.sub(r"[^0-9A-Za-z_.]", "_", str(value)) or "label" for value in values]


def materialize_block(matrix):
    if sparse.issparse(matrix):
        return matrix.toarray()
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
        if sparse.issparse(matrix):
            return matrix.toarray()
    return np.asarray(matrix)


def optional_positive_int(value):
    if value in {None, "", "null", "NULL"}:
        return None
    value = int(value)
    return value if value > 0 else None


def is_missing(value):
    return value in {None, "", "null", "NULL"}


def stratified_sample_positions(obs, positions, strata_column, size, rng):
    strata = obs.iloc[positions][strata_column].astype(str).to_numpy()
    grouped = [
        group.to_numpy()
        for _stratum, group in pd.Series(positions).groupby(strata, sort=False)
    ]
    selected = []
    remaining = size
    for index, group in enumerate(sorted(grouped, key=len)):
        groups_left = len(grouped) - index
        target = int(np.ceil(remaining / groups_left))
        take = min(len(group), target)
        if take < len(group):
            group = rng.choice(group, size=take, replace=False)
        selected.extend(group.tolist())
        remaining -= take
    return np.asarray(selected, dtype=int)


def selected_cell_positions(
    obs,
    cell_type_column,
    max_cells_per_cell_type,
    sampling_strata_column,
    seed,
):
    if max_cells_per_cell_type is None:
        return np.arange(len(obs))

    if not is_missing(sampling_strata_column) and sampling_strata_column not in obs:
        raise ValueError(
            f"input.h5ad_sampling_strata_column '{sampling_strata_column}' not found"
        )

    rng = np.random.default_rng(seed)
    selected = []
    cell_types = obs[cell_type_column].astype(str)
    for _cell_type, positions in pd.Series(np.arange(len(obs))).groupby(
        cell_types.to_numpy(), sort=False
    ):
        positions = positions.to_numpy()
        if len(positions) > max_cells_per_cell_type:
            if is_missing(sampling_strata_column):
                positions = rng.choice(
                    positions, size=max_cells_per_cell_type, replace=False
                )
            else:
                positions = stratified_sample_positions(
                    obs,
                    positions,
                    sampling_strata_column,
                    max_cells_per_cell_type,
                    rng,
                )
        selected.extend(positions.tolist())

    return np.asarray(sorted(selected), dtype=int)


def matrix_layer(adata, layer):
    if layer in {None, "", "null", "NULL"}:
        return adata.X
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' was not found in the h5ad file")
    return adata.layers[layer]


def gene_names(adata, gene_symbol_column):
    if gene_symbol_column in {None, "", "null", "NULL"}:
        names = pd.Index(adata.var_names.astype(str))
    else:
        if gene_symbol_column not in adata.var:
            raise ValueError(f"var column '{gene_symbol_column}' was not found")
        names = pd.Index(adata.var[gene_symbol_column].astype(str))

    counts = {}
    deduplicated = []
    for name in names:
        count = counts.get(name, 0)
        counts[name] = count + 1
        deduplicated.append(name if count == 0 else f"{name}.{count + 1}")
    return pd.Index(deduplicated)


def write_annotation(obs, output_path, cell_ids, cell_type_column, sample_column):
    missing = [
        column for column in [cell_type_column, sample_column] if column not in obs
    ]
    if missing:
        raise ValueError("Missing h5ad obs column(s): " + ", ".join(missing))

    annotation = obs.copy()
    annotation.insert(0, "ID", cell_ids)
    annotation["CellType"] = make_safe_labels(annotation[cell_type_column])
    annotation["Sample"] = annotation[sample_column].astype(str).to_numpy()

    leading = ["ID", "CellType", "Sample"]
    extras = [column for column in annotation.columns if column not in leading]
    annotation = annotation[leading + extras]
    annotation.to_csv(output_path, sep="\t", index=False)
    return annotation


def write_expression(matrix, output_path, genes, cells, chunk_genes):
    n_obs, n_vars = matrix.shape
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write("Gene\t" + "\t".join(cells) + "\n")
        for start in range(0, n_vars, chunk_genes):
            stop = min(start + chunk_genes, n_vars)
            block = materialize_block(matrix[:, start:stop])
            if block.shape != (n_obs, stop - start):
                block = np.asarray(block).reshape(n_obs, stop - start)

            for offset, gene in enumerate(genes[start:stop]):
                values = block[:, offset]
                row = "\t".join(f"{value:.12g}" for value in values)
                handle.write(f"{gene}\t{row}\n")


h5ad_path = Path(snakemake.input["h5ad"])
matrix_output = Path(snakemake.output["matrix"])
annotation_output = Path(snakemake.output["annotation"])
selected_cells_output = Path(snakemake.output["selected_cells"])
layer = snakemake.params.get("layer")
cell_type_column = snakemake.params["cell_type_column"]
sample_column = snakemake.params["sample_column"]
gene_symbol_column = snakemake.params.get("gene_symbol_column")
chunk_genes = int(snakemake.params["chunk_genes"])
max_cells_per_cell_type = optional_positive_int(
    snakemake.params.get("max_cells_per_cell_type")
)
sampling_seed = int(snakemake.params["sampling_seed"])
sampling_strata_column = snakemake.params.get("sampling_strata_column")
sanitize_ids = as_bool(snakemake.params["sanitize_ids"])

if chunk_genes < 1:
    raise ValueError("h5ad_chunk_genes must be at least 1")

matrix_output.parent.mkdir(parents=True, exist_ok=True)
annotation_output.parent.mkdir(parents=True, exist_ok=True)
selected_cells_output.parent.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(h5ad_path, backed="r")
matrix = matrix_layer(adata, layer)
selected_positions = selected_cell_positions(
    adata.obs,
    cell_type_column,
    max_cells_per_cell_type,
    sampling_strata_column,
    sampling_seed,
)
selected_obs = adata.obs.iloc[selected_positions].copy()
source_cells = adata.obs_names[selected_positions].astype(str)
cells = make_safe_ids(source_cells) if sanitize_ids else list(source_cells)
genes = gene_names(adata, gene_symbol_column)
matrix = matrix[selected_positions, :]

annotation = write_annotation(
    selected_obs, annotation_output, cells, cell_type_column, sample_column
)
selected_manifest = annotation[["ID", "CellType", "Sample"]].copy()
selected_manifest.insert(0, "SourceID", source_cells)
selected_manifest.to_csv(selected_cells_output, sep="\t", index=False)
write_expression(matrix, matrix_output, genes, cells, chunk_genes)
