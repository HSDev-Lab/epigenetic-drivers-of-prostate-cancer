#!/usr/bin/env Rscript

NEURONAL_CUTOFF <- 0.8
RHOGTPASE_CUTOFF <- 1.0
MEVOTYPE_LEVELS <- c("MEvotype2", "MEvotype1")
MEVOTYPE_COLOURS <- c("MEvotype2" = "#5CAE63", "MEvotype1" = "#986EAC")

package_path <- function(script_dir) {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

load_programme_data <- function(input_file) {
  if (!file.exists(input_file)) stop("Missing input file: ", input_file)
  dat <- data.table::fread(input_file)
  required <- c(
    "patient", "methylation_sample", "mevotype", "neuronal_score",
    "rho_score", "rfs_years", "rfs_event", "grade_group", "advanced_stage"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  stopifnot(
    nrow(dat) == 163L,
    !anyDuplicated(dat$patient),
    sum(dat$rfs_event) == 43L,
    all(dat$mevotype %in% MEVOTYPE_LEVELS),
    !anyNA(dat$neuronal_score),
    !anyNA(dat$rho_score)
  )
  dat[, mevotype := factor(mevotype, levels = MEVOTYPE_LEVELS)]
  dat[, neuronal_high := neuronal_score >= NEURONAL_CUTOFF]
  dat[, rho_high := rho_score >= RHOGTPASE_CUTOFF]
  dat[, high_programme_count := as.integer(neuronal_high) + as.integer(rho_high)]
  data.table::setorder(dat, patient)
  dat
}
