import json
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


def optional_positive_int(value):
    if value in {None, "", "null", "NULL"}:
        return None
    value = int(value)
    return value if value > 0 else None


def is_missing(value):
    return value in {None, "", "null", "NULL"}


def safe_id(value):
    safe = re.sub(r"[^0-9A-Za-z_.]", ".", str(value))
    if re.match(r"^[0-9]", safe):
        safe = f"X{safe}"
    return safe or "cell"


def safe_label(value):
    return re.sub(r"[^0-9A-Za-z_.]", "_", str(value)) or "label"


def make_unique(values):
    seen = {}
    out = []
    for value in values:
        count = seen.get(value, 0)
        seen[value] = count + 1
        out.append(value if count == 0 else f"{value}.{count + 1}")
    return out


def fractions_name(filter_genes):
    if filter_genes == "cell type specific":
        raise ValueError(
            "Direct h5ad cell-state input mode does not support "
            "pipeline.filter_genes: cell type specific. Use no filter and "
            "input.h5ad_feature_selection instead."
        )
    if filter_genes == "no filter":
        return "All_genes"
    n_genes = optional_positive_int(filter_genes)
    if n_genes is None:
        raise ValueError(f"Invalid pipeline.filter_genes value: {filter_genes}")
    return f"Top_{n_genes}"


def matrix_layer(adata, layer):
    if layer in {None, "", "null", "NULL"}:
        return adata.X
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' was not found in the h5ad file")
    return adata.layers[layer]


def materialize_block(matrix):
    if sparse.issparse(matrix):
        return matrix.toarray()
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
        if sparse.issparse(matrix):
            return matrix.toarray()
    return np.asarray(matrix)


def subset_block(matrix, rows, start, stop):
    block = matrix[rows, start:stop]
    block = materialize_block(block).astype(np.float32, copy=False)
    return np.nan_to_num(block, copy=False, nan=0.0, posinf=0.0, neginf=0.0)


def gene_names(adata, gene_symbol_column):
    if gene_symbol_column in {None, "", "null", "NULL"}:
        names = pd.Index(adata.var_names.astype(str))
    else:
        if gene_symbol_column not in adata.var:
            raise ValueError(f"var column '{gene_symbol_column}' was not found")
        names = pd.Index(adata.var[gene_symbol_column].astype(str))
    return pd.Index(make_unique(names.tolist()))


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
    rng = np.random.default_rng(seed)
    if (
        max_cells_per_cell_type is not None
        and not is_missing(sampling_strata_column)
        and sampling_strata_column not in obs
    ):
        raise ValueError(
            f"input.h5ad_sampling_strata_column '{sampling_strata_column}' not found"
        )

    selected = []
    cell_types = obs[cell_type_column].astype(str)
    for _cell_type, positions in pd.Series(np.arange(len(obs))).groupby(
        cell_types.to_numpy(), sort=False
    ):
        positions = positions.to_numpy()
        if (
            max_cells_per_cell_type is not None
            and len(positions) > max_cells_per_cell_type
        ):
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


def read_annotation(path):
    if path in {None, "", "null", "NULL"}:
        return None
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def obs_with_optional_annotation(
    adata,
    annotation_file,
    annotation_id_column,
    obs_id_column,
):
    obs = adata.obs.copy()
    obs["_h5ad_pos"] = np.arange(adata.n_obs)
    if annotation_file in {None, "", "null", "NULL"}:
        return obs.reset_index(drop=False).rename(columns={"index": "_h5ad_obs_name"})

    annotation = read_annotation(annotation_file)
    if annotation_id_column not in annotation.columns:
        raise ValueError(
            f"input.h5ad_annotation_id_column '{annotation_id_column}' "
            "was not found in input.h5ad_annotation_file"
        )

    if obs_id_column in {None, "", "null", "NULL"}:
        obs["_h5ad_join_id"] = adata.obs_names.astype(str)
    else:
        if obs_id_column not in obs.columns:
            raise ValueError(f"input.h5ad_obs_id_column '{obs_id_column}' not found")
        obs["_h5ad_join_id"] = obs[obs_id_column].astype(str)

    annotation = annotation.drop_duplicates(annotation_id_column).copy()
    annotation["_h5ad_join_id"] = annotation[annotation_id_column].astype(str)
    obs = obs.merge(
        annotation,
        on="_h5ad_join_id",
        how="inner",
        suffixes=("", "_annotation"),
    )
    if obs.empty:
        raise ValueError("External h5ad annotation join produced zero cells")
    return obs.reset_index(drop=True)


def gene_statistics(matrix, rows, chunk_genes):
    n_vars = matrix.shape[1]
    sums = np.zeros(n_vars, dtype=np.float64)
    sums_sq = np.zeros(n_vars, dtype=np.float64)
    nonzero = np.zeros(n_vars, dtype=np.int64)
    for start in range(0, n_vars, chunk_genes):
        stop = min(start + chunk_genes, n_vars)
        block = subset_block(matrix, rows, start, stop)
        sums[start:stop] = block.sum(axis=0)
        sums_sq[start:stop] = np.square(block, dtype=np.float64).sum(axis=0)
        nonzero[start:stop] = (block > 0).sum(axis=0)

    n_obs = max(len(rows), 1)
    means = sums / n_obs
    variances = np.maximum((sums_sq / n_obs) - np.square(means), 0)
    return nonzero, variances


def selected_gene_indices(
    matrix, rows, feature_selection, n_top_genes, min_cells, chunk
):
    nonzero, variances = gene_statistics(matrix, rows, chunk)
    valid = np.where((nonzero >= min_cells) & np.isfinite(variances))[0]
    valid = valid[variances[valid] > 0]
    if len(valid) == 0:
        raise ValueError("No genes passed direct h5ad feature filters")

    feature_selection = str(feature_selection).lower()
    if feature_selection in {"hvg", "variance", "top_variable"}:
        if n_top_genes is None:
            raise ValueError("input.h5ad_n_top_genes is required for hvg selection")
        ranked = valid[np.argsort(variances[valid])]
        return np.sort(ranked[-min(n_top_genes, len(ranked)) :])
    if feature_selection in {"all", "none"}:
        if n_top_genes is not None:
            ranked = valid[np.argsort(variances[valid])]
            return np.sort(ranked[-min(n_top_genes, len(ranked)) :])
        return valid

    raise ValueError(f"Unsupported input.h5ad_feature_selection: {feature_selection}")


def write_table(path, values, row_names, column_names):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(column_names) + "\n")
        for row_name, row_values in zip(row_names, values, strict=True):
            row = "\t".join(f"{value:.12g}" for value in row_values)
            handle.write(f"{row_name}\t{row}\n")


def scale_rows(values, groups=None):
    scaled = np.zeros(values.shape, dtype=np.float32)
    if groups is None:
        groups = np.repeat("all", values.shape[1])

    groups = pd.Series(groups.astype(str) if hasattr(groups, "astype") else groups)
    for group in groups.drop_duplicates():
        columns = np.where(groups.to_numpy() == group)[0]
        if len(columns) < 2:
            continue
        block = values[:, columns]
        means = block.mean(axis=1, keepdims=True)
        std = block.std(axis=1, ddof=1, keepdims=True)
        std[std == 0] = np.nan
        scaled[:, columns] = (block - means) / std
    return np.nan_to_num(scaled, copy=False, nan=0.0, posinf=0.0, neginf=0.0)


def cell_correlation(scaled):
    if scaled.shape[1] < 2:
        raise ValueError("At least two cells are required to compute correlations")
    corr = np.corrcoef(scaled.T)
    corr = np.nan_to_num(corr, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr


def prepare_cell_type(
    matrix,
    positions,
    selected_gene_positions,
    selected_gene_names,
    cell_ids,
    scale_groups,
    output_dir,
):
    blocks = []
    for start in range(0, len(selected_gene_positions), 512):
        stop = min(start + 512, len(selected_gene_positions))
        gene_idx = selected_gene_positions[start:stop]
        block = materialize_block(matrix[positions[:, None], gene_idx])
        blocks.append(block.astype(np.float32, copy=False))

    expression = np.concatenate(blocks, axis=1)
    expression = np.nan_to_num(expression, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    log_data = np.log2(expression + 1.0).T
    scaled = scale_rows(log_data, scale_groups)
    corr = cell_correlation(scaled)

    write_table(
        output_dir / "expression_full_matrix_log2.txt",
        log_data,
        selected_gene_names,
        cell_ids,
    )
    write_table(
        output_dir / "expression_full_matrix_scaled.txt",
        scaled,
        selected_gene_names,
        cell_ids,
    )
    write_table(
        output_dir / "expression_top_genes_scaled_filt.txt",
        scaled,
        selected_gene_names,
        cell_ids,
    )
    write_table(
        output_dir / "expression_top_genes_scaled.txt",
        corr,
        cell_ids,
        cell_ids,
    )


h5ad_path = Path(snakemake.input["h5ad"])
annotation_file = snakemake.input.get("annotation")
if isinstance(annotation_file, list):
    annotation_file = annotation_file[0] if annotation_file else None
complete_output = Path(snakemake.output["complete"])
annotation_output = Path(snakemake.output["annotation"])
selected_cells_output = Path(snakemake.output["selected_cells"])

discovery = snakemake.params["discovery_dataset_name"]
proc_dir = Path(snakemake.params["proc_dir"])
results_dir = Path(snakemake.params["results_dir"])
layer = snakemake.params.get("layer")
annotation_id_column = snakemake.params["annotation_id_column"]
obs_id_column = snakemake.params.get("obs_id_column")
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
feature_selection = snakemake.params["feature_selection"]
n_top_genes = optional_positive_int(snakemake.params.get("n_top_genes"))
min_cells_per_gene = int(snakemake.params["min_cells_per_gene"])
min_cells_per_cell_type = int(snakemake.params["min_cells_per_cell_type"])
filter_genes = snakemake.params["filter_genes"]
scale_by = snakemake.params.get("scale_by")

if chunk_genes < 1:
    raise ValueError("input.h5ad_chunk_genes must be at least 1")
if min_cells_per_gene < 0:
    raise ValueError("input.h5ad_min_cells_per_gene must be at least 0")
if min_cells_per_cell_type < 1:
    raise ValueError("input.h5ad_min_cells_per_cell_type must be at least 1")

fractions = fractions_name(filter_genes)
work_root = (
    results_dir
    / "ecotyper_work"
    / discovery
    / fractions
    / "Cell_States"
    / "discovery_cross_cor"
)
annotation_output.parent.mkdir(parents=True, exist_ok=True)
selected_cells_output.parent.mkdir(parents=True, exist_ok=True)
complete_output.parent.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(h5ad_path, backed="r")
obs_frame = obs_with_optional_annotation(
    adata,
    annotation_file,
    annotation_id_column,
    obs_id_column,
)
missing_columns = [
    column for column in [cell_type_column, sample_column] if column not in obs_frame
]
if missing_columns:
    raise ValueError("Missing h5ad obs column(s): " + ", ".join(missing_columns))
if scale_by not in {None, "", "null", "NULL"} and scale_by not in obs_frame.columns:
    raise ValueError(f"input.scale_by column '{scale_by}' was not found in h5ad obs")

matrix = matrix_layer(adata, layer)
obs_frame = obs_frame[
    obs_frame[cell_type_column].notna()
    & ~obs_frame[cell_type_column].astype(str).isin({"", "NA", "nan", "None"})
].reset_index(drop=True)
selected_frame_positions = selected_cell_positions(
    obs_frame,
    cell_type_column,
    max_cells_per_cell_type,
    sampling_strata_column,
    sampling_seed,
)
selected_obs = obs_frame.iloc[selected_frame_positions].copy()
selected_positions = selected_obs["_h5ad_pos"].to_numpy(dtype=int)
selected_obs["CellType"] = [
    safe_label(value) for value in selected_obs[cell_type_column]
]
selected_obs["Sample"] = selected_obs[sample_column].astype(str).to_numpy()

cell_type_counts = selected_obs["CellType"].value_counts()
kept_cell_types = cell_type_counts[cell_type_counts >= min_cells_per_cell_type].index
selected_mask = selected_obs["CellType"].isin(kept_cell_types).to_numpy()
selected_positions = selected_positions[selected_mask]
selected_obs = selected_obs.loc[selected_mask].reset_index(drop=True).copy()

source_cells = adata.obs_names[selected_positions].astype(str)
cell_ids = (
    make_unique([safe_id(value) for value in source_cells])
    if sanitize_ids
    else list(source_cells)
)
if "ID" in selected_obs.columns:
    selected_obs = selected_obs.drop(columns=["ID"])
selected_obs.insert(0, "ID", cell_ids)
annotation = selected_obs.copy()
leading = ["ID", "CellType", "Sample"]
internal_columns = {"_h5ad_pos", "_h5ad_join_id", "_h5ad_obs_name"}
extras = [
    column
    for column in annotation.columns
    if column not in leading and column not in internal_columns
]
annotation = annotation[leading + extras]
annotation.to_csv(annotation_output, sep="\t", index=False)

selected_manifest = annotation[["ID", "CellType", "Sample"]].copy()
selected_manifest.insert(0, "SourceID", source_cells)
selected_manifest.to_csv(selected_cells_output, sep="\t", index=False)

gene_index = selected_gene_indices(
    matrix,
    selected_positions,
    feature_selection,
    n_top_genes,
    min_cells_per_gene,
    chunk_genes,
)
all_gene_names = gene_names(adata, gene_symbol_column)
selected_gene_names = all_gene_names[gene_index].to_list()

summaries = []
for cell_type, group in annotation.groupby("CellType", sort=True):
    relative_positions = group.index.to_numpy()
    positions = selected_positions[relative_positions]
    ids = group["ID"].astype(str).to_list()
    scale_groups = None
    if scale_by not in {None, "", "null", "NULL"}:
        scale_groups = group[scale_by].astype(str).to_numpy()
    output_dir = work_root / cell_type
    prepare_cell_type(
        matrix,
        positions,
        gene_index,
        selected_gene_names,
        ids,
        scale_groups,
        output_dir,
    )
    summaries.append({"CellType": cell_type, "n_cells": len(group)})

summary = {
    "source_h5ad": str(h5ad_path),
    "discovery_dataset_name": discovery,
    "fractions": fractions,
    "n_cells": int(len(annotation)),
    "n_genes": int(len(gene_index)),
    "cell_types": summaries,
    "feature_selection": str(feature_selection),
    "n_top_genes": n_top_genes,
    "min_cells_per_gene": min_cells_per_gene,
    "min_cells_per_cell_type": min_cells_per_cell_type,
    "max_cells_per_cell_type": max_cells_per_cell_type,
    "sampling_strata_column": None
    if is_missing(sampling_strata_column)
    else str(sampling_strata_column),
}
complete_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
