#!/usr/bin/env Rscript

# Render the interactive QC dashboard (Rmd -> HTML).
#
# CLI adaptation of upstream snpArcher modules/qc/scripts/qc_dashboard_render.R
# (which read snakemake params/output globals).  v2.3.0 threads min_depth into
# the Rmd params (the About panel quotes it) and renders straight to the
# declared output with knit intermediates in a private directory: rendering
# beside the Rmd wrote into the workflow's own scripts directory, so projects
# run from one clone collided and a read-only install failed.
#
# Usage:
#   qc_dashboard_render.R <qc_dir> <n_clusters> <min_depth> <google_api_key> <has_qc_report> <rmd_path> <dashboard_output>

render_qcplots <- function(qc_dir, nClusters, minDepth, GMKey, has_qc_report, rmd_path, output_path) {
    workd <- getwd()

    script.in <- rmd_path

    intermediates <- tempfile("qc_dashboard_")
    dir.create(intermediates)
    on.exit(unlink(intermediates, recursive = TRUE), add = TRUE)

    rmarkdown::render(script.in,
                      params = list(qc_dir = qc_dir,
                                    nClusters = nClusters,
                                    minDepth = minDepth,
                                    GMKey = GMKey,
                                    has_qc_report = has_qc_report),
                      knit_root_dir = workd,
                      output_dir = dirname(output_path),
                      output_file = basename(output_path),
                      intermediates_dir = intermediates)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Usage: qc_dashboard_render.R <qc_dir> <n_clusters> <min_depth> <google_api_key> <has_qc_report> <rmd_path> <dashboard_output>")
}

render_qcplots(args[1], as.integer(args[2]), as.numeric(args[3]), args[4], as.logical(args[5]), args[6], args[7])
