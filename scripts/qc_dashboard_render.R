#!/usr/bin/env Rscript

# Render the interactive QC dashboard (Rmd -> HTML).
#
# CLI adaptation of upstream snpArcher modules/qc/scripts/qc_dashboard_render.R
# (which read snakemake params/output globals).  The rendering logic and the
# downstream qc_dashboard_interactive.Rmd are unchanged.
#
# Usage:
#   qc_dashboard_render.R <qc_dir> <n_clusters> <google_api_key> <has_qc_report> <rmd_path> <dashboard_output>

render_qcplots <- function(qc_dir, nClusters, GMKey, has_qc_report, rmd_path, output_path) {
    workd <- getwd()

    script.in <- rmd_path
    script.out <- gsub(".Rmd", ".html", script.in)

    rmarkdown::render(script.in,
                      params = list(qc_dir = qc_dir,
                                    nClusters = nClusters,
                                    GMKey = GMKey,
                                    has_qc_report = has_qc_report),
                      knit_root_dir = workd)

    copy_successful <- file.copy(script.out, output_path)

    if (copy_successful) {
        file.remove(script.out)
    } else {
        cat("snpArcher: Failed to move the qc dashboard html.\n")
    }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("Usage: qc_dashboard_render.R <qc_dir> <n_clusters> <google_api_key> <has_qc_report> <rmd_path> <dashboard_output>")
}

render_qcplots(args[1], as.integer(args[2]), args[3], as.logical(args[4]), args[5], args[6])
