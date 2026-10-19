#!/usr/bin/env Rscript
# Recreate candidate figures from an already verified R run WITHOUT refitting PCA.
if (!file.exists(file.path("R","03_figures.R"))) stop("Run from project root")
source(file.path("R","03_figures.R"))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0L) {
  run_dir <- normalizePath(args[1L],mustWork=TRUE)
} else {
  run_id_file <- file.path("runs","LATEST_VERIFIED_PUBLIC.txt")
  if (!file.exists(run_id_file)) stop("No verified run yet. Execute Rscript run_all.R")
  run_id <- trimws(readLines(run_id_file,warn=FALSE)[1L])
  if (!grepl("^PUBLIC_CANDIDATE_R_[0-9]{8}_[0-9]{6}_[0-9]+$",run_id)) stop("Invalid last-run identifier")
  run_dir <- normalizePath(file.path("runs",run_id),mustWork=TRUE)
}
if (!file.exists(file.path(run_dir,"PUBLIC_CANDIDATE_VERIFICATION_PASS.txt"))) {
  stop("Cannot render unverified R run")
}
revision_id <- paste0("figure_revision_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
revision_subdir <- file.path("figure_revisions",revision_id)
render_figures_from_outputs(run_dir, revision_subdir)
message("New figure revision created (without modifying validated plots or rerunning PCA): ",
        file.path(run_dir,revision_subdir))
