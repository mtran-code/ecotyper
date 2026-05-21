suppressPackageStartupMessages({
  library(doParallel)
  library(NMF)
})

run_torch_nmf <- function(input_file, output_dir, n_clusters, seed) {
  python = Sys.which("python")
  if (python == "") {
    stop(
      "Python was not found on PATH. Run through pixi or set PATH before using torch NMF."
    )
  }

  script = file.path("run_nmf_torch.py")
  status = system2(
    python,
    args = c(
      script,
      "--input",
      input_file,
      "--output-dir",
      output_dir,
      "--rank",
      as.character(n_clusters),
      "--seed",
      as.character(seed),
      "--device",
      Sys.getenv("ECOTYPER_NMF_TORCH_DEVICE", "auto"),
      "--max-iter",
      Sys.getenv("ECOTYPER_NMF_MAX_ITER", "500"),
      "--tolerance",
      Sys.getenv("ECOTYPER_NMF_TOLERANCE", "0.0001")
    )
  )
  if (!identical(status, 0L)) {
    stop("Torch NMF backend failed. See the Python traceback above.")
  }

  W = as.matrix(read.delim(
    file.path(output_dir, "torch_W.txt"),
    row.names = 1,
    check.names = F
  ))
  H = as.matrix(read.delim(
    file.path(output_dir, "torch_H.txt"),
    row.names = 1,
    check.names = F
  ))
  storage.mode(W) = "double"
  storage.mode(H) = "double"

  model = nmfModel(W = W, H = H)
  methods::new(
    "NMFfit",
    fit = model,
    method = paste0("torch:", Sys.getenv("ECOTYPER_NMF_TORCH_DEVICE", "auto")),
    seed = as.character(seed),
    residuals = NA_real_
  )
}

args = c(
  "discovery",
  "scRNA_CRC_Park",
  "scRNA_specific_genes",
  "CD4.T.cells",
  "7",
  "2"
)
args = commandArgs(T)
dataset_type = args[1]
dataset = args[2]
fractions = args[3]
cell_type = args[4]
n_clusters = as.integer(as.character(args[5]))
restart = as.integer(as.character(args[6]))

input_dir = file.path(
  "../../data/results/ecotyper_work",
  dataset,
  fractions,
  "Cell_States",
  dataset_type,
  cell_type
)

if (!file.exists(file.path(input_dir, "expression_top_genes_scaled.txt"))) {
  stop(paste0("No input data for cell type: ", cell_type))
}

output_dir = file.path(
  "../../data/results/ecotyper_work",
  dataset,
  fractions,
  "Cell_States",
  dataset_type,
  cell_type,
  n_clusters,
  "restarts",
  restart
)
dir.create(output_dir, recursive = T, showWarning = F)

input_file = file.path(input_dir, "expression_top_genes_scaled.txt")

cat(paste0(
  "Running NMF on '",
  cell_type,
  "' (number of states = ",
  n_clusters,
  ", restart ",
  restart,
  ")...\n"
))
seed = 1234 + restart
nmf_backend = tolower(Sys.getenv("ECOTYPER_NMF_BACKEND", "r"))
if (nmf_backend %in% c("torch", "gpu", "auto")) {
  estim.r <- run_torch_nmf(input_file, output_dir, n_clusters, seed)
} else {
  raw_data = read.delim(input_file)
  data = posneg(as.matrix(raw_data))
  estim.r <- nmf(
    data,
    n_clusters,
    nrun = 1,
    method = "brunet",
    seed = seed,
    .opt = 'P1'
  )
}
save(estim.r, file = file.path(output_dir, "estim.RData"))
