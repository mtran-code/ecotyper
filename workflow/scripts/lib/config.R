library(data.table)

configure_nmf_backend <- function(config) {
  settings = config$"Pipeline settings"

  backend = settings$"NMF backend"
  if (is.null(backend) || is.na(backend)) {
    backend = "r"
  }

  torch_device = settings$"NMF torch device"
  if (is.null(torch_device) || is.na(torch_device)) {
    torch_device = "auto"
  }

  max_iter = settings$"NMF max iterations"
  if (is.null(max_iter) || is.na(max_iter)) {
    max_iter = 500
  }

  tolerance = settings$"NMF tolerance"
  if (is.null(tolerance) || is.na(tolerance)) {
    tolerance = 0.0001
  }

  Sys.setenv(
    ECOTYPER_NMF_BACKEND = as.character(backend),
    ECOTYPER_NMF_TORCH_DEVICE = as.character(torch_device),
    ECOTYPER_NMF_MAX_ITER = as.character(max_iter),
    ECOTYPER_NMF_TOLERANCE = as.character(tolerance)
  )
}

check_discovery_configuration <- function(config) {
  input_mat = config$Input$"Expression matrix"
  discovery = config$Input$"Discovery dataset name"
  annotation = config$Input$"Annotation file"
  output_dir = file.path("data/procdata/datasets/discovery", discovery)
  CSx_singularity_path_fractions = config$"Pipeline settings"$"CIBERSORTx fractions Singularity path"
  CSx_singularity_path_hires = config$"Pipeline settings"$"CIBERSORTx hires Singularity path"

  if (!file.exists(input_mat)) {
    stop(paste0("Input format error: Input file '", input_mat, "' is missing!"))
  }

  mat = fread(input_mat, sep = "\t", nrows = 5)
  if (ncol(mat) < 2) {
    stop(paste0(
      ncol(mat),
      " columns detected in file '",
      input_mat,
      "'. Please make sure that the file is tab-delimited!"
    ))
  }

  dir.create(output_dir, recursive = T, showWarning = F)
  system(paste0(
    "ln -sf '",
    normalizePath(input_mat),
    "' '",
    file.path(output_dir, "data.txt"),
    "'"
  ))

  if (!is.null(annotation)) {
    if (!file.exists(annotation)) {
      stop(paste0(
        "Input format error: Annotation file '",
        annotation,
        "' is missing!"
      ))
    }
    system(paste0(
      "ln -sf '",
      normalizePath(annotation),
      "' '",
      file.path(output_dir, "annotation.txt"),
      "'"
    ))
  }

  if (
    !is.na(CSx_singularity_path_fractions) &&
      !is.null(CSx_singularity_path_fractions) &&
      !file.exists(CSx_singularity_path_fractions)
  ) {
    stop(paste0(
      "CIBERSORTx fractions Singularity path provided does not exist:",
      CSx_singularity_path_fractions
    ))
  }
  if (
    !is.na(CSx_singularity_path_hires) &&
      !is.null(CSx_singularity_path_hires) &&
      !file.exists(CSx_singularity_path_hires)
  ) {
    stop(paste0(
      "CIBERSORTx hires Singularity path provided does not exist:",
      CSx_singularity_path_hires
    ))
  }
}

check_discovery_configuration_scRNA <- function(config) {
  input_mat = config$Input$"Expression matrix"
  discovery = config$Input$"Discovery dataset name"
  annotation = config$Input$"Annotation file"
  use_prepared_cell_state_inputs = as.logical(
    config$"Pipeline settings"$"Use prepared cell state inputs"
  )
  if (is.na(use_prepared_cell_state_inputs)) {
    use_prepared_cell_state_inputs = F
  }
  p_value_cutoff = as.numeric(as.character(
    config$"Pipeline settings"$"Jaccard matrix p-value cutoff"
  ))
  output_dir = file.path("data/procdata/datasets/discovery", discovery)

  if (!use_prepared_cell_state_inputs && !file.exists(input_mat)) {
    stop(paste0("Input format error: Input file '", input_mat, "' is missing!"))
  }

  if (
    length(p_value_cutoff) == 0 ||
      is.na(p_value_cutoff) ||
      p_value_cutoff <= 0 ||
      p_value_cutoff > 1
  ) {
    stop(paste0(
      "The p-value cutoff in field 'Jaccard matrix p-value cutoff' needs to be a number in the iterval (0,1]. Value provided:",
      p_value_cutoff,
      "."
    ))
  }

  dir.create(output_dir, recursive = T, showWarning = F)

  if (!use_prepared_cell_state_inputs) {
    mat = fread(input_mat, sep = "\t", nrows = 5)
    if (ncol(mat) < 2) {
      stop(paste0(
        ncol(mat),
        " columns detected in file '",
        input_mat,
        "'. Please make sure that the file is tab-delimited!"
      ))
    }
    system(paste0(
      "ln -sf '",
      normalizePath(input_mat),
      "' '",
      file.path(output_dir, "data.txt"),
      "'"
    ))
  }

  if (!file.exists(annotation)) {
    stop(paste0(
      "Input format error: Annotation file '",
      annotation,
      "' is missing! This file needs to be provided for discovery in scRNA-seq data!"
    ))
  }

  ann = read.delim(annotation)
  if (!all(c("ID", "CellType", "Sample") %in% colnames(ann))) {
    stop(paste0(
      "Input format error: Annotation file '",
      annotation,
      "' does not contain columns: 'ID', 'CellType' or 'Sample'! All three columns are required for discovery in scRNA-seq data!"
    ))
  }

  if (!all(as.character(ann$ID) == make.names(as.character(ann$ID)))) {
    stop(paste0(
      "Input format error: The values in column 'ID' of the annotation file '",
      annotation,
      "' contain special characters (e.g. ' ', '-'), or start with digits. Please make sure that this column and the column names of the expression matrix do not contain values modified by the R function 'make.names'."
    ))
  }

  if (!use_prepared_cell_state_inputs && !all(colnames(mat)[-1] %in% ann$ID)) {
    stop(paste0(
      "Input format error: The following ids present in the column names of the expression matrix are missing from the annotation file (column 'ID'): '",
      paste(
        colnames(mat)[-1][!colnames(mat)[-1] %in% ann$ID],
        collapse = "', '"
      ),
      "'."
    ))
  }

  if (
    normalizePath(annotation) !=
      normalizePath(file.path(output_dir, "annotation.txt"))
  ) {
    system(paste0(
      "ln -sf '",
      normalizePath(annotation),
      "' '",
      file.path(output_dir, "annotation.txt"),
      "'"
    ))
  }
}

check_discovery_configuration_presorted <- function(config) {
  input_path = config$Input$"Expression matrices"
  discovery = config$Input$"Discovery dataset name"
  annotation = config$Input$"Annotation file"

  if (config$"Pipeline settings"$"Filter genes" == "cell type specific") {
    fractions = "Cell_type_specific_genes"
  } else {
    if (config$"Pipeline settings"$"Filter genes" == "no filter") {
      fractions = "All_genes"
    } else {
      n_genes = as.integer(as.numeric(
        config$"Pipeline settings"$"Filter genes"
      ))
      fractions = paste0("Top_", n_genes)
    }
  }

  output_dir = file.path("data/procdata/datasets/discovery", discovery)
  csx_dir = file.path("data/results/cibersortx/hires", discovery, fractions)
  dir.create(output_dir, recursive = T, showWarning = F)
  dir.create(csx_dir, recursive = T, showWarning = F)

  if (!is.null(annotation)) {
    if (!file.exists(annotation)) {
      stop(paste0(
        "Input format error: Annotation file '",
        annotation,
        "' is missing!"
      ))
    }
    system(paste0(
      "ln -sf '",
      normalizePath(annotation),
      "' '",
      file.path(output_dir, "annotation.txt"),
      "'"
    ))
  }

  classes = NULL
  for (file in list.files(input_path)) {
    cell_type = gsub(".txt$", "", file)
    input_mat = file.path(input_path, file)
    if (!file.exists(input_mat)) {
      stop(paste0(
        "Input format error: Input file '",
        input_mat,
        "' is missing!"
      ))
    }

    mat = fread(input_mat, sep = "\t")
    if (ncol(mat) < 2) {
      stop(paste0(
        ncol(mat),
        " columns detected in file '",
        input_mat,
        "'. Please make sure that the file is tab-delimited!"
      ))
    }
    classes = rbind(classes, data.frame(x = cell_type))
    system(paste0(
      "ln -sf '",
      normalizePath(input_mat),
      "' '",
      file.path(csx_dir, paste0(cell_type, ".txt")),
      "'"
    ))
  }
  write.table(
    t(classes),
    file.path(csx_dir, "classes.txt"),
    sep = "\t",
    row.names = F,
    col.names = F
  )
}
