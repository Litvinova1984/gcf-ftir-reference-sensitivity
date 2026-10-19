#!/usr/bin/env Rscript
# Packaging regression only: checks files, R syntax, and finalization on a temporary fixture.
# Does NOT call compute_results(), fit_adult(), project_fixed(), or any statistical function.

if (!file.exists("PACKAGE_SHA256.csv") || !file.exists("R/07_release_io.R"))
  stop("Run from the root of the unpacked RC2 repository")
root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "R", "00_core.R"))
source(file.path(root, "R", "02_verify.R"))
source(file.path(root, "R", "07_release_io.R"))
run_id <- paste0("PACKAGING_CHECK_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
out <- file.path(root, "runs", run_id)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lines <- character()
record <- function(s) { lines <<- c(lines, s); message(s) }
expect_error <- function(expr) {
  failed <- FALSE
  tryCatch(force(expr), error = function(e) { failed <<- TRUE })
  if (!failed) stop("Negative packaging test unexpectedly accepted invalid files")
  invisible(TRUE)
}

passed <- tryCatch({
  m <- utils::read.csv(file.path(root, "PACKAGE_SHA256.csv"), stringsAsFactors = FALSE,
                       check.names = FALSE)
  if (!identical(names(m), c("relative_path", "sha256", "size_bytes")) || !nrow(m) ||
      anyNA(m) || anyDuplicated(m$relative_path) ||
      any(!grepl("^[0-9a-fA-F]{64}$", m$sha256))) stop("Invalid package manifest")
  if (any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", m$relative_path))) stop("Invalid manifest path")
  paths <- file.path(root, m$relative_path)
  if (any(!file.exists(paths)) || any(file.info(paths)$isdir)) stop("Package has missing files")
  if (any(file.info(paths)$size != m$size_bytes)) stop("Package file size mismatch")
  actual <- vapply(paths, file_sha256, character(1))
  bad <- which(actual != tolower(m$sha256))
  if (length(bad)) stop("PACKAGE SHA256 MISMATCH: ", paste(m$relative_path[bad], collapse = ", "))
  record(paste("PASS package hashes:", nrow(m), "files"))
  verify_candidate_input_hashes(root)
  record("PASS all five original data/grid/demographic input hashes")
  rfiles <- m$relative_path[grepl("[.]R$", m$relative_path)]
  invisible(lapply(file.path(root, rfiles), function(f) parse(file = f)))
  record(paste("PASS actual R parser:", length(rfiles), "scripts"))
  dependency <- "verification/frozen_R_20261018/PRIVATE_pediatric_Q_decomposition_K3.csv"
  ignore <- readLines(file.path(root, ".gitignore"), warn = FALSE)
  if (!(dependency %in% m$relative_path) || !(paste0("!", dependency) %in% ignore))
    stop("Required verified checkpoint or its exact gitignore exception is missing")
  record("PASS required verification CSV and narrowly scoped gitignore exception")

  fixture <- tempfile("mdpi_packaging_fixture_")
  dir.create(file.path(fixture, "runs", "fixture"), recursive = TRUE)
  f <- file.path(fixture, "runs", "fixture")
  tryCatch({
    # A tiny NON-SCIENTIFIC fixture to test the exact I/O routine used by run_all.R.
    writeLines("packaging fixture; no participant or research results", file.path(f, "fixture.txt"))
    writeLines("START packaging-only fixture", file.path(f, "run_log.txt"))
    logger <- function(s) cat(s, "\n", file = file.path(f, "run_log.txt"), append = TRUE)
    finalize_verified_run(fixture, f, "fixture", logger)
    verify_run_manifest(f)
    if (!identical(readLines(file.path(fixture, "runs", "LATEST_VERIFIED_PUBLIC.txt")), "fixture"))
      stop("Latest-run pointer did not update after verification")
    record("PASS final log and pass flag both included in complete RUN_SHA256")
    # Appending to a sealed log must be detected, unlike the original bug.
    cat("unexpected late append\n", file = file.path(f, "run_log.txt"), append = TRUE)
    expect_error(verify_run_manifest(f))
    record("PASS negative test: late log changes rejected")
    write_sha_manifest(f, file.path(f, "RUN_SHA256.csv"))
    writeLines("unlisted file", file.path(f, "unexpected.txt"))
    expect_error(verify_run_manifest(f))
    unlink(file.path(f, "unexpected.txt"))
    record("PASS negative test: unlisted run file rejected")
    unlink(file.path(f, "PUBLIC_CANDIDATE_VERIFICATION_PASS.txt"))
    expect_error(verify_run_manifest(f))
    record("PASS negative test: missing pass flag rejected")
    writeLines("simulated failure", file.path(f, "RUN_FAILED.txt"))
    expect_error(finalize_verified_run(fixture, f, "fixture", logger))
    record("PASS negative test: failed run cannot be sealed as successful")
    # An I/O error during finalization must not leave a successful pass flag.
    g <- file.path(fixture, "runs", "io_error_fixture")
    dir.create(g)
    writeLines("START", file.path(g, "run_log.txt"))
    bad_logger <- function(s) stop("simulated log write error")
    expect_error(finalize_verified_run(fixture, g, "io_error_fixture", bad_logger))
    if (file.exists(file.path(g, "PUBLIC_CANDIDATE_VERIFICATION_PASS.txt")))
      stop("I/O failure left an invalid pass flag")
    record("PASS negative test: I/O failure removes incomplete pass flag")
  }, finally = unlink(fixture, recursive = TRUE))
  writeLines(lines, file.path(out, "CHECK_LOG.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(c("PACKAGING CHECK PASS", "No PCA or statistical calculations executed.",
               "This does not grant publication approval."), file.path(out, "PACKAGING_CHECK_PASS.txt"))
  write_sha_manifest(out, file.path(out, "RUN_SHA256.csv"))
  verify_run_manifest(out, required = c("CHECK_LOG.txt", "sessionInfo.txt", "PACKAGING_CHECK_PASS.txt"))
  writeLines(run_id, file.path(root, "runs", "LATEST_PACKAGING_CHECK.txt"))
  TRUE
}, error = function(e) {
  message("PACKAGING CHECK FAILED: ", conditionMessage(e))
  unlink(file.path(out, "PACKAGING_CHECK_PASS.txt"))
  writeLines(c(lines, conditionMessage(e)), file.path(out, "PACKAGING_CHECK_FAILED.txt"))
  FALSE
})
if (!passed) quit(status = 1L)
message("PACKAGING CHECK PASS - no PCA, no new analysis. Results: ", out)
