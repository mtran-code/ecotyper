packages = c("NMF", "HiClimR")

missing_packages <- function() {
  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
}

install_with_pak <- function(packages) {
  if (!requireNamespace("pak", quietly = TRUE)) {
    return(FALSE)
  }

  tryCatch(
    {
      pak::pkg_install(packages)
      TRUE
    },
    error = function(error) {
      message("pak installation failed: ", conditionMessage(error))
      FALSE
    }
  )
}

missing = missing_packages()
if (length(missing) == 0) {
  quit(save = "no")
}

if (!install_with_pak(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

missing = missing_packages()
if (length(missing) > 0) {
  stop("Failed to install R package(s): ", paste(missing, collapse = ", "))
}
