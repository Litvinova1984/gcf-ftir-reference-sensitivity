# Packaging-only I/O. No PCA, projections, statistics, or changes to input data.
# Uses the unchanged file_sha256(), write_csv(), and write_sha_manifest() helpers.

verify_run_manifest <- function(folder, required = c("run_log.txt", "PUBLIC_CANDIDATE_VERIFICATION_PASS.txt")) {
  folder <- normalizePath(folder, mustWork = TRUE)
  manifest <- file.path(folder, "RUN_SHA256.csv")
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(m), c("relative_path", "sha256")) || !nrow(m) ||
      anyNA(m) || anyDuplicated(m$relative_path) ||
      any(!grepl("^[0-9a-fA-F]{64}$", m$sha256))) stop("Invalid run manifest")
  if (any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", m$relative_path)))
    stop("Run manifest must use relative paths inside the run")
  files <- list.files(folder, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  rel <- substring(files, nchar(folder) + 2L)
  rel <- rel[rel != "RUN_SHA256.csv"]
  if (!setequal(rel, m$relative_path)) stop("Run manifest coverage differs from actual files")
  if (!all(required %in% m$relative_path)) stop("Required final files absent from run manifest")
  actual <- vapply(file.path(folder, m$relative_path), file_sha256, character(1))
  bad <- which(actual != tolower(m$sha256))
  if (length(bad)) stop("RUN SHA256 MISMATCH: ", paste(m$relative_path[bad], collapse = ", "))
  invisible(nrow(m))
}

finalize_verified_run <- function(root, out, run_id, logline) {
  # Call ONLY after all numerical verification steps have succeeded.
  # This routine seals files; it does not itself verify scientific calculations.
  root <- normalizePath(root, mustWork = TRUE)
  out <- normalizePath(out, mustWork = TRUE)
  if (basename(out) != run_id || dirname(out) != file.path(root, "runs"))
    stop("Run path and run identifier do not agree")
  if (file.exists(file.path(out, "RUN_FAILED.txt"))) stop("Cannot finalize a failed run")
  pass_file <- file.path(out, "PUBLIC_CANDIDATE_VERIFICATION_PASS.txt")
  manifest <- file.path(out, "RUN_SHA256.csv")
  completed <- FALSE
  on.exit({
    if (!completed) {
      unlink(pass_file)
      # A partial/stale manifest must not be presented as a sealed successful run.
      if (file.exists(manifest) && !dir.exists(manifest)) unlink(manifest)
    }
  }, add = TRUE)
  writeLines(c("R CANDIDATE FULL DOWNSTREAM REPRODUCTION PASS",
               paste("Run", run_id),
               "Start from normalized adult 18x793 and pediatric 25x793 CSVs",
               "PASS: 12 baseline numerical comparisons + full K=3 diagnostics + metadata",
               "NOT COVERED: upstream DPT/OriginPro or pediatric interpolation/SNV",
               "NOT A PUBLICATION APPROVAL: author declarations and release approval remain separate"),
             pass_file)
  logline("SUCCESS: verified candidate, NOT public release authorization")
  # The pass flag and final log line are now complete. Nothing below changes a
  # file INSIDE the run directory after its manifest is written.
  write_sha_manifest(out, manifest)
  verify_run_manifest(out)
  # This convenience pointer is OUTSIDE the sealed run directory.
  writeLines(run_id, file.path(root, "runs", "LATEST_VERIFIED_PUBLIC.txt"))
  completed <- TRUE
  invisible(out)
}
