#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}

scripts <- c("01_generate_figure1.R", "02_generate_figure2.R")
for (script in scripts) {
  status <- system2("Rscript", file.path(script_dir, script))
  if (!identical(status, 0L)) stop("Failed: ", script)
}

package_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
writeLines(capture.output(sessionInfo()), file.path(package_dir, "sessionInfo.txt"))
message("Generated both final figures and sessionInfo.txt")
