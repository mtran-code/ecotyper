from pathlib import Path

import yaml


def optional(value):
    return None if value in {"", "null", "NULL"} else value


def legacy_common(config, output_dir):
    input_cfg = config["input"]
    pipeline_cfg = config["pipeline"]
    return {
        "Input": {
            "Discovery dataset name": input_cfg["discovery_dataset_name"],
            "Expression matrix": optional(input_cfg.get("expression_matrix")),
            "Annotation file": optional(input_cfg.get("annotation_file")),
            "Annotation file column to scale by": optional(input_cfg.get("scale_by")),
            "Annotation file column(s) to plot": input_cfg.get("plot_columns", []),
        },
        "Output": {
            "Output folder": output_dir,
        },
        "Pipeline settings": {
            "Pipeline steps to skip": pipeline_cfg.get("steps_to_skip", []),
            "Filter genes": pipeline_cfg.get("filter_genes", "cell type specific"),
            "Number of threads": pipeline_cfg.get("threads", 10),
            "Number of NMF restarts": pipeline_cfg.get("nmf_restarts", 5),
            "Maximum number of states per cell type": pipeline_cfg.get(
                "max_states_per_cell_type", 20
            ),
            "Cophenetic coefficient cutoff": pipeline_cfg.get(
                "cophenetic_cutoff", 0.975
            ),
            "Jaccard matrix p-value cutoff": pipeline_cfg.get(
                "jaccard_p_value_cutoff", 1
            ),
            "Minimum number of states in ecotypes": pipeline_cfg.get(
                "min_states_in_ecotypes", 3
            ),
        },
    }


def render_discovery_scrna(config, output_dir):
    legacy = legacy_common(config, output_dir)
    legacy["Input"]["Expression type"] = config["input"].get("expression_type", "scRNA")
    return legacy


def render_discovery_bulk(config, output_dir):
    legacy = legacy_common(config, output_dir)
    input_cfg = config["input"]
    pipeline_cfg = config["pipeline"]
    legacy["Input"].update(
        {
            "CIBERSORTx username": optional(input_cfg.get("cibersortx_username")),
            "CIBERSORTx token": optional(input_cfg.get("cibersortx_token")),
            "Cell type fractions": optional(input_cfg.get("cell_type_fractions")),
        }
    )
    legacy["Pipeline settings"].update(
        {
            "Filter non cell type specific genes": pipeline_cfg.get(
                "filter_non_cell_type_specific_genes", False
            ),
            "CIBERSORTx fractions Singularity path": optional(
                pipeline_cfg.get("cibersortx_fractions_singularity")
            ),
            "CIBERSORTx hires Singularity path": optional(
                pipeline_cfg.get("cibersortx_hires_singularity")
            ),
        }
    )
    return legacy


def render_discovery_presorted(config, output_dir):
    legacy = legacy_common(config, output_dir)
    input_cfg = config["input"]
    pipeline_cfg = config["pipeline"]
    legacy["Input"].pop("Expression matrix", None)
    legacy["Input"]["Expression matrices"] = optional(
        input_cfg.get("expression_matrices")
    )
    legacy["Pipeline settings"]["Filter non cell type specific genes"] = (
        pipeline_cfg.get("filter_non_cell_type_specific_genes", False)
    )
    return legacy


RENDERERS = {
    "discovery_scRNA": render_discovery_scrna,
    "discovery_bulk": render_discovery_bulk,
    "discovery_presorted": render_discovery_presorted,
}

mode = snakemake.params["mode"]
run_name = snakemake.params["run_name"]
results_dir = Path(snakemake.params["results_dir"])
output_dir = str(results_dir / run_name)

config = dict(snakemake.config)

legacy = {"default": RENDERERS[mode](config, output_dir)}

output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)
with open(output_path, "w", encoding="utf-8") as handle:
    yaml.safe_dump(legacy, handle, sort_keys=False)
