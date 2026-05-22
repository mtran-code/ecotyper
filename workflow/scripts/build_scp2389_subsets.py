#!/usr/bin/env python3
"""Build donor-balanced SCP2389 h5ad subsets for EcoTyper discovery."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


RETAINED_LINEAGES = (
    "Malignant",
    "Myeloid",
    "Oligo",
    "Tcells",
    "Stromal",
    "Other_Immune",
)


@dataclass(frozen=True)
class CohortSpec:
    name: str
    label: str
    output_stem: str


COHORTS = (
    CohortSpec(
        name="scp2389_gbm_idhwt",
        label="GBM IDHwt",
        output_stem="scp2389_gbm_idhwt_ecotyper_discovery",
    ),
    CohortSpec(
        name="scp2389_idhmut",
        label="combined IDHmut glioma",
        output_stem="scp2389_idhmut_ecotyper_discovery",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create QCed, donor-balanced SCP2389 h5ad subsets from sparse "
            "Matrix Market counts without loading the full Seurat object."
        )
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        required=True,
        help="Path to the local SCP2389 bulk-download directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where subset h5ads and summaries will be written.",
    )
    parser.add_argument(
        "--max-cells-per-donor-cell-type",
        type=int,
        default=150,
        help="Maximum cells retained from each donor/cell-type stratum.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=500,
        help="Minimum detected genes per cell.",
    )
    parser.add_argument(
        "--min-umi",
        type=int,
        default=1000,
        help="Minimum UMI count per cell.",
    )
    parser.add_argument(
        "--max-mito",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage per cell.",
    )
    parser.add_argument(
        "--min-cells-per-gene",
        type=int,
        default=10,
        help="Drop genes detected in fewer selected cells than this threshold.",
    )
    parser.add_argument(
        "--normalization-scale",
        type=float,
        default=10000.0,
        help="Library-size normalization scale for h5ad X.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260522,
        help="Random seed for donor/cell-type downsampling.",
    )
    parser.add_argument(
        "--chunk-nnz",
        type=int,
        default=1_000_000,
        help="Number of selected nonzero entries buffered before flushing.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove existing output directory before rebuilding.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    elapsed = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{elapsed}] {message}", flush=True)


def read_metadata(path: Path) -> pd.DataFrame:
    numeric_columns = [
        "Genes_detected",
        "UMI",
        "Mitochondial_Percentage",
        "Malignant",
        "Myeloid",
        "Oligo",
        "Stromal",
        "Tcells",
        "Other_Immune",
        "Cycling",
        "Age",
        "Grade",
    ]
    metadata = pd.read_csv(path, sep="\t", skiprows=[1], low_memory=False)
    for column in numeric_columns:
        if column in metadata.columns:
            metadata[column] = pd.to_numeric(metadata[column], errors="coerce")
    return metadata


def cohort_mask(metadata: pd.DataFrame, spec: CohortSpec) -> pd.Series:
    if spec.name == "scp2389_gbm_idhwt":
        return metadata["Diagnosis"].eq("Glioblastoma") & metadata[
            "IDH_Mutation_Status"
        ].eq("WT")
    if spec.name == "scp2389_idhmut":
        return metadata["IDH_Mutation_Status"].eq("Mut") & metadata["Diagnosis"].isin(
            ["Astrocytoma", "Oligodendroglioma"]
        )
    raise ValueError(f"Unhandled cohort: {spec.name}")


def qc_mask(metadata: pd.DataFrame, args: argparse.Namespace) -> pd.Series:
    return (
        metadata["Genes_detected"].ge(args.min_genes)
        & metadata["UMI"].ge(args.min_umi)
        & metadata["Mitochondial_Percentage"].le(args.max_mito)
        & metadata["Total_Tumor_Annotation"].isin(RETAINED_LINEAGES)
    )


def sample_cohort(
    metadata: pd.DataFrame,
    mask: pd.Series,
    max_per_stratum: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    eligible = metadata.loc[mask].copy()
    selected_indices: list[int] = []
    strata = eligible.groupby(["donor_id", "Total_Tumor_Annotation"], sort=True)
    for _stratum, group in strata:
        indices = group.index.to_numpy()
        if len(indices) > max_per_stratum:
            indices = rng.choice(indices, size=max_per_stratum, replace=False)
        selected_indices.extend(indices.tolist())
    sampled = metadata.loc[selected_indices].copy()
    return sampled.sort_values("NAME").reset_index(drop=True)


def read_single_column(path: Path) -> list[str]:
    values: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            values.append(line.rstrip("\n\r").split("\t")[0])
    return values


def read_genes(path: Path) -> pd.DataFrame:
    genes = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "gene_symbol"])
    genes["gene_id"] = genes["gene_id"].astype(str)
    genes["gene_symbol"] = genes["gene_symbol"].astype(str)
    return genes


def make_unique(values: pd.Series) -> list[str]:
    seen: dict[str, int] = {}
    unique = []
    for value in values.astype(str):
        count = seen.get(value, 0)
        seen[value] = count + 1
        unique.append(value if count == 0 else f"{value}.{count + 1}")
    return unique


def flush_chunk(
    temp_dir: Path,
    cohort_name: str,
    chunk_index: int,
    rows: list[int],
    cols: list[int],
    counts: list[int],
) -> Path:
    path = temp_dir / cohort_name / f"chunk_{chunk_index:05d}.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        rows=np.asarray(rows, dtype=np.int32),
        cols=np.asarray(cols, dtype=np.int32),
        counts=np.asarray(counts, dtype=np.int32),
    )
    return path


def read_matrix_header(handle) -> tuple[int, int, int]:
    for raw in handle:
        line = raw.decode("utf-8")
        if line.startswith("%"):
            continue
        n_genes, n_cells, n_nonzero = [int(value) for value in line.split()]
        return n_genes, n_cells, n_nonzero
    raise ValueError("Matrix Market file ended before the size header")


def stream_selected_counts(
    matrix_path: Path,
    cohort_by_col: np.ndarray,
    new_col_by_col: np.ndarray,
    cohort_names: list[str],
    n_genes: int,
    temp_dir: Path,
    chunk_nnz: int,
) -> tuple[dict[str, list[Path]], dict[str, np.ndarray]]:
    buffers = {
        name: {"rows": [], "cols": [], "counts": [], "chunk_index": 0}
        for name in cohort_names
    }
    chunks = {name: [] for name in cohort_names}
    gene_counts = {name: np.zeros(n_genes, dtype=np.int32) for name in cohort_names}
    selected_entries = 0
    processed_entries = 0
    started = time.time()

    with gzip.open(matrix_path, "rb") as handle:
        header_genes, header_cells, n_nonzero = read_matrix_header(handle)
        if header_genes != n_genes:
            raise ValueError(
                f"Gene count mismatch: genes.tsv has {n_genes}, matrix has {header_genes}"
            )
        if header_cells + 1 != len(cohort_by_col):
            raise ValueError(
                f"Cell count mismatch: barcodes.tsv has {len(cohort_by_col) - 1}, "
                f"matrix has {header_cells}"
            )

        for raw in handle:
            processed_entries += 1
            gene_raw, cell_raw, count_raw = raw.split()
            old_col = int(cell_raw)
            cohort_code = cohort_by_col[old_col]
            if cohort_code >= 0:
                row = int(gene_raw) - 1
                name = cohort_names[int(cohort_code)]
                buffer = buffers[name]
                buffer["rows"].append(row)
                buffer["cols"].append(int(new_col_by_col[old_col]))
                buffer["counts"].append(int(count_raw))
                gene_counts[name][row] += 1
                selected_entries += 1
                if len(buffer["rows"]) >= chunk_nnz:
                    chunk_path = flush_chunk(
                        temp_dir,
                        name,
                        int(buffer["chunk_index"]),
                        buffer["rows"],
                        buffer["cols"],
                        buffer["counts"],
                    )
                    chunks[name].append(chunk_path)
                    buffer["rows"].clear()
                    buffer["cols"].clear()
                    buffer["counts"].clear()
                    buffer["chunk_index"] = int(buffer["chunk_index"]) + 1

            if processed_entries % 50_000_000 == 0:
                minutes = (time.time() - started) / 60
                pct = processed_entries / n_nonzero * 100
                log(
                    f"streamed {processed_entries:,}/{n_nonzero:,} matrix entries "
                    f"({pct:.1f}%); retained {selected_entries:,}; {minutes:.1f} min"
                )

    for name, buffer in buffers.items():
        if buffer["rows"]:
            chunk_path = flush_chunk(
                temp_dir,
                name,
                int(buffer["chunk_index"]),
                buffer["rows"],
                buffer["cols"],
                buffer["counts"],
            )
            chunks[name].append(chunk_path)

    log(
        f"completed matrix stream; retained {selected_entries:,} selected nonzero "
        f"entries from {processed_entries:,}"
    )
    return chunks, gene_counts


def build_sparse_matrices(
    chunk_paths: list[Path],
    gene_keep: np.ndarray,
    umi: np.ndarray,
    shape: tuple[int, int],
    normalization_scale: float,
) -> tuple[sparse.csr_matrix, sparse.csr_matrix]:
    gene_map = np.full(len(gene_keep), -1, dtype=np.int32)
    gene_map[np.where(gene_keep)[0]] = np.arange(int(gene_keep.sum()), dtype=np.int32)
    rows_blocks = []
    cols_blocks = []
    counts_blocks = []

    for path in chunk_paths:
        chunk = np.load(path)
        mapped_genes = gene_map[chunk["rows"]]
        keep = mapped_genes >= 0
        if not np.any(keep):
            continue
        rows_blocks.append(chunk["cols"][keep].astype(np.int32, copy=False))
        cols_blocks.append(mapped_genes[keep].astype(np.int32, copy=False))
        counts_blocks.append(chunk["counts"][keep].astype(np.int32, copy=False))

    if not rows_blocks:
        raise ValueError("No selected matrix entries passed gene filters")

    rows = np.concatenate(rows_blocks)
    cols = np.concatenate(cols_blocks)
    counts = np.concatenate(counts_blocks)
    safe_umi = np.maximum(umi[rows].astype(np.float32), 1.0)
    cp10k = counts.astype(np.float32) / safe_umi * np.float32(normalization_scale)
    matrix_shape = (shape[0], int(gene_keep.sum()))
    normalized = sparse.coo_matrix((cp10k, (rows, cols)), shape=matrix_shape).tocsr()
    counts_matrix = sparse.coo_matrix(
        (counts, (rows, cols)), shape=matrix_shape
    ).tocsr()
    return normalized, counts_matrix


def write_h5ad(
    output_path: Path,
    spec: CohortSpec,
    sampled: pd.DataFrame,
    genes: pd.DataFrame,
    gene_keep: np.ndarray,
    normalized: sparse.csr_matrix,
    counts: sparse.csr_matrix,
    args: argparse.Namespace,
    eligible_n: int,
) -> None:
    obs = sampled.copy()
    obs_names = obs["NAME"].astype(str).to_numpy()
    obs.index = pd.Index(obs_names, name=None)
    obs["ecotyper_cell_type"] = obs["Total_Tumor_Annotation"].astype(str)
    obs["ecotyper_sample"] = obs["donor_id"].astype(str)
    obs["scp2389_subset"] = spec.name
    obs["cohort_label"] = spec.label

    var = genes.loc[gene_keep].copy().reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"]), name=None)
    var["source_gene_index"] = np.where(gene_keep)[0].astype(np.int32)

    for column in obs.columns:
        if obs[column].dtype == object:
            obs[column] = obs[column].fillna("NA").astype(str)

    adata = ad.AnnData(X=normalized, obs=obs, var=var)
    adata.layers["counts"] = counts
    adata.uns["scp2389_subset"] = {
        "cohort": spec.name,
        "label": spec.label,
        "eligible_cells_after_qc": int(eligible_n),
        "selected_cells": int(adata.n_obs),
        "selected_genes": int(adata.n_vars),
        "normalization": f"raw counts per {args.normalization_scale:g} UMIs",
        "qc": {
            "min_genes": args.min_genes,
            "min_umi": args.min_umi,
            "max_mito": args.max_mito,
            "retained_lineages": list(RETAINED_LINEAGES),
            "excluded_lineages": ["Cycling"],
            "min_cells_per_gene": args.min_cells_per_gene,
            "max_cells_per_donor_cell_type": args.max_cells_per_donor_cell_type,
            "seed": args.seed,
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")


def write_summaries(
    output_dir: Path,
    summary_rows: list[dict[str, object]],
    sampled_frames: dict[str, pd.DataFrame],
) -> None:
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        output_dir / "scp2389_ecotyper_subset_summary.tsv", sep="\t", index=False
    )

    workbook = {}
    for name, sampled in sampled_frames.items():
        workbook[name] = {
            "cells_by_cell_type": sampled["Total_Tumor_Annotation"]
            .value_counts()
            .sort_index()
            .to_dict(),
            "donors_by_cell_type": sampled.groupby("Total_Tumor_Annotation")["donor_id"]
            .nunique()
            .sort_index()
            .to_dict(),
            "cells_by_diagnosis": sampled["Diagnosis"]
            .value_counts()
            .sort_index()
            .to_dict(),
            "cells_by_primary_recurrent": sampled["Primary_Recurrent_Status"]
            .value_counts()
            .sort_index()
            .to_dict(),
        }
    (output_dir / "scp2389_ecotyper_subset_summary.json").write_text(
        json.dumps(workbook, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    source_dir = args.source_dir
    output_dir = args.output_dir
    temp_dir = output_dir / "_tmp_triplets"
    raw_dir = source_dir / "expression" / "66f44f68d5b3cc83dd766a0f"
    metadata_path = source_dir / "metadata" / "GBM2_MetaData_Single_Cell_Portal.txt"
    barcodes_path = raw_dir / "barcodes.tsv"
    genes_path = raw_dir / "genes.tsv"
    matrix_path = raw_dir / "matrix.mtx.gz"

    if args.force and output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log("reading SCP2389 metadata and sparse matrix indexes")
    metadata = read_metadata(metadata_path)
    genes = read_genes(genes_path)
    barcodes = read_single_column(barcodes_path)
    barcode_to_col = {barcode: index + 1 for index, barcode in enumerate(barcodes)}
    if len(metadata) != len(barcodes):
        raise ValueError(
            f"Metadata/barcode row mismatch: metadata={len(metadata)}, "
            f"barcodes={len(barcodes)}"
        )

    rng = np.random.default_rng(args.seed)
    base_qc = qc_mask(metadata, args)
    sampled_frames: dict[str, pd.DataFrame] = {}
    eligible_counts: dict[str, int] = {}
    cohort_by_col = np.full(len(barcodes) + 1, -1, dtype=np.int8)
    new_col_by_col = np.full(len(barcodes) + 1, -1, dtype=np.int32)
    cohort_names = [spec.name for spec in COHORTS]
    umi_by_cohort: dict[str, np.ndarray] = {}

    for cohort_code, spec in enumerate(COHORTS):
        mask = base_qc & cohort_mask(metadata, spec)
        eligible_counts[spec.name] = int(mask.sum())
        sampled = sample_cohort(
            metadata,
            mask,
            args.max_cells_per_donor_cell_type,
            rng,
        )
        old_cols = sampled["NAME"].map(barcode_to_col)
        if old_cols.isna().any():
            missing = sampled.loc[old_cols.isna(), "NAME"].head().to_list()
            raise ValueError(
                f"{spec.name} selected barcodes missing from matrix: {missing}"
            )
        old_cols_array = old_cols.to_numpy(dtype=np.int32)
        sampled["_matrix_col"] = old_cols_array
        sampled = sampled.sort_values("_matrix_col").reset_index(drop=True)
        sampled["_subset_col"] = np.arange(len(sampled), dtype=np.int32)

        cohort_by_col[sampled["_matrix_col"].to_numpy(dtype=np.int32)] = cohort_code
        new_col_by_col[sampled["_matrix_col"].to_numpy(dtype=np.int32)] = sampled[
            "_subset_col"
        ].to_numpy(dtype=np.int32)
        sampled_frames[spec.name] = sampled
        umi_by_cohort[spec.name] = sampled["UMI"].to_numpy(dtype=np.float32)

        log(
            f"{spec.name}: {eligible_counts[spec.name]:,} eligible cells after QC; "
            f"{len(sampled):,} selected across {sampled['donor_id'].nunique()} donors"
        )

    chunks, gene_counts = stream_selected_counts(
        matrix_path,
        cohort_by_col,
        new_col_by_col,
        cohort_names,
        len(genes),
        temp_dir,
        args.chunk_nnz,
    )

    summary_rows: list[dict[str, object]] = []
    for spec in COHORTS:
        sampled = sampled_frames[spec.name]
        gene_keep = gene_counts[spec.name] >= args.min_cells_per_gene
        log(
            f"{spec.name}: building h5ad with {len(sampled):,} cells, "
            f"{int(gene_keep.sum()):,} genes, and {len(chunks[spec.name])} chunks"
        )
        normalized, counts = build_sparse_matrices(
            chunks[spec.name],
            gene_keep,
            umi_by_cohort[spec.name],
            (len(sampled), int(gene_keep.sum())),
            args.normalization_scale,
        )
        output_path = output_dir / f"{spec.output_stem}.h5ad"
        write_h5ad(
            output_path,
            spec,
            sampled.drop(columns=["_matrix_col", "_subset_col"]),
            genes,
            gene_keep,
            normalized,
            counts,
            args,
            eligible_counts[spec.name],
        )
        summary_rows.append(
            {
                "cohort": spec.name,
                "label": spec.label,
                "h5ad": str(output_path),
                "eligible_cells_after_qc": eligible_counts[spec.name],
                "selected_cells": len(sampled),
                "selected_genes": int(gene_keep.sum()),
                "donors": sampled["donor_id"].nunique(),
                "biosamples": sampled["biosample_id"].nunique(),
            }
        )
        log(f"{spec.name}: wrote {output_path}")

    write_summaries(output_dir, summary_rows, sampled_frames)
    shutil.rmtree(temp_dir, ignore_errors=True)
    log("SCP2389 EcoTyper subset build complete")


if __name__ == "__main__":
    main()
