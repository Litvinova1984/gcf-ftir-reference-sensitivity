# Strict R-versus-frozen-Python, per-record/per-model numerical verifier.
# The frozen tables are checkpoints ONLY, never PCA training or projection input.

file_sha256 <- function(file) {
  if (!file.exists(file)) stop("Missing file for SHA-256: ", file)
  sha_bin <- Sys.which("shasum")
  if (nzchar(sha_bin)) {
    output <- suppressWarnings(system2(sha_bin, args = c("-a", "256", shQuote(file)),
                                      stdout = TRUE, stderr = TRUE))
  } else {
    sha_bin <- Sys.which("sha256sum")
    if (!nzchar(sha_bin)) stop("Neither shasum nor sha256sum found; SHA-256 gate cannot run")
    output <- suppressWarnings(system2(sha_bin, args = shQuote(file),
                                      stdout = TRUE, stderr = TRUE))
  }
  if (length(output) < 1L || !grepl("^[0-9a-fA-F]{64} ",output[1L])) {
    stop("Invalid SHA-256 tool response for: ",file, " ; ",paste(output,collapse="; "))
  }
  tolower(substr(output[1L],1L,64L))
}

verify_candidate_input_hashes <- function(root) {
  f <- file.path(root,"verification","PUBLIC_INPUT_SHA256.csv")
  m <- utils::read.csv(f,check.names=FALSE,stringsAsFactors=FALSE)
  if (!identical(names(m),c("relative_path","sha256")) || anyDuplicated(m$relative_path))
    stop("Invalid normalized input hash manifest")
  for (i in seq_len(nrow(m))) {
    actual <- file_sha256(file.path(root,m$relative_path[i]))
    if (!identical(actual,tolower(m$sha256[i])))
      stop("INPUT SHA256 MISMATCH: ",m$relative_path[i])
  }
  message("PASS: ",nrow(m)," candidate inputs SHA256 verified")
  invisible(nrow(m))
}

compare_table <- function(root, out_dir, frozen_name, generated_name, keys,
                          restrict = NULL, columns = NULL) {
  frozen_path <- file.path(root,"verification","frozen_R_20261017",frozen_name)
  generated_path <- file.path(out_dir,generated_name)
  if (!file.exists(frozen_path) || !file.exists(generated_path)) {
    stop("Missing cross-language comparison pair: ",frozen_name," / ",generated_name)
  }
  f <- utils::read.csv(frozen_path,check.names=FALSE,stringsAsFactors=FALSE)
  r <- utils::read.csv(generated_path,check.names=FALSE,stringsAsFactors=FALSE)
  if (!is.null(restrict)) {
    f <- restrict(f)
    r <- restrict(r)
  }
  if (!all(keys %in% names(f)) || !all(keys %in% names(r))) {
    stop("Missing key columns: ",frozen_name)
  }
  # Robust composite key with no locale/lexical coercion ambiguities for current IDs.
  key_for <- function(df) do.call(paste,c(lapply(df[keys],as.character),sep="|"))
  if (identical(keys, "wavenumber_cm_1")) {
    if (nrow(f)!=nrow(r) ||
        max(abs(f$wavenumber_cm_1-r$wavenumber_cm_1)) > 1e-8) {
      stop("Full-precision wave grid/order differs from frozen source: ",frozen_name)
    }
    fk <- seq_len(nrow(f))
  } else {
    fk <- key_for(f)
    rk <- key_for(r)
    if (anyDuplicated(fk) || anyDuplicated(rk) ||
        length(fk)!=length(rk) || !setequal(fk,rk)) stop("Row IDs/model keys do not match: ",frozen_name)
    r <- r[match(fk,rk),,drop=FALSE]
  }
  if (is.null(columns)) columns <- setdiff(names(f), keys)
  if (!all(columns %in% names(r))) stop("Missing numeric or text fields in ",generated_name)
  maxima <- numeric()
  for (col in columns) {
    a <- f[[col]]; b <- r[[col]]
    if (is.numeric(a) && is.numeric(b)) {
      if (anyNA(a) || anyNA(b) || any(!is.finite(a)) || any(!is.finite(b))) {
        stop("Nonfinite value in comparison: ",generated_name,":",col)
      }
      maxdiff <- max(abs(a-b))
      # Strict absolute gates; principal angles allow harmless microdegree SVD drift.
      tol <- if (grepl("angle",col,ignore.case=TRUE)) 1e-4 else 1e-6
      if (maxdiff > tol) {
        i <- which.max(abs(a-b))
        stop("NUMERICAL MISMATCH ",generated_name," field ",col,
             " key ",fk[i]," ; frozen=",format(a[i],digits=16),
             " ; R=",format(b[i],digits=16),
             " ; max abs difference=",format(maxdiff,digits=12)," ; tolerance=",tol)
      }
      maxima[col] <- maxdiff
    } else {
      if (!identical(as.character(a),as.character(b))) {
        stop("Non-numeric field mismatch: ",generated_name,":",col)
      }
    }
  }
  data.frame(comparison=frozen_name,rows=nrow(f),fields=length(columns),
             maximum_abs_discrepancy=if(length(maxima)) max(maxima) else 0)
}

verify_all <- function(root, out_dir) {
  # All 25 children, full 793-bin SNV, 17 adult LOO, 5 reference models * 25,
  # all K=2:5 * models, PCA mean, loading vectors, and baseline eigenvalues.
  jobs <- list(
   list("pediatric_on_adult_grid_SNV.csv","pediatric_on_adult_grid_SNV.csv",c("sample_id")),
   list("adult17_PCA_center.csv","adult17_PCA_center.csv",c("wavenumber_cm_1")),
   list("adult17_PCA_loadings_PC1_PC5.csv","adult17_PCA_loadings_PC1_PC5.csv",c("wavenumber_cm_1")),
   list("pediatric_external_scores_K3.csv","pediatric_external_scores_K3.csv",c("pediatric_sample_id")),
   list("pediatric_model_17_vs_18_K3.csv","pediatric_model_17_vs_18_K3.csv",c("pediatric_sample_id")),
   list("adult_loo_reference_K3.csv","adult_loo_reference_K3.csv",c("anonymous_adult_index")),
   list("adult_leave_one_out_influence_aggregated.csv","adult_leave_one_out_influence_aggregated.csv",c("reference_row_anonymous")),
   list("component_sensitivity_K2_K5.csv","component_sensitivity_K2_K5.csv",c("K")),
   list("model_comparison_K3.csv","model_comparison_K3.csv",c("model")),
   list("pediatric_per_model_K3.csv","pediatric_per_model_K3.csv",c("model","pediatric_id")),
   list("model_sensitivity_K2_K5.csv","model_sensitivity_K2_K5.csv",c("model","K")))
  results <- do.call(rbind,lapply(jobs,function(j) {
    compare_table(root,out_dir,j[[1]],j[[2]],j[[3]])
  }))
  # Historical eigenvalues include all components; this R pipeline exports first 5.
  ev <- compare_table(root,out_dir,"adult_reference_eigenvalues_17_18_PC1_PC5.csv",
                      "adult_reference_eigenvalues_17_18_PC1_PC5.csv",
                      c("model","component"))
  results <- rbind(results,ev)
  # Extra defensive checks: all models have every pediatric ID exactly once.
  x <- utils::read.csv(file.path(out_dir,"pediatric_per_model_K3.csv"),stringsAsFactors=FALSE)
  expected_models <- c("adult17_main","adult16_no_A04","adult16_no_A07",
                       "adult15_no_A04_A07","archival18_mixed_age")
  if (!setequal(unique(x$model), expected_models)) stop("Unexpected sensitivity model inventory")
  for (nm in expected_models) {
    ids <- as.character(x$pediatric_id[x$model==nm])
    if (length(ids)!=25L || anyDuplicated(ids) || !setequal(ids, as.character(1:25))) {
      stop("Missing pediatric projections for reference model ",nm)
    }
  }
  z <- utils::read.csv(file.path(out_dir,"model_sensitivity_K2_K5.csv"))
  for (nm in expected_models) {
    if (!identical(sort(z$K[z$model==nm]),2:5)) stop("Incomplete K=2:5 series: ",nm)
  }
  write_csv(results,file.path(out_dir,"cross_language_verification.csv"))
  message("PASS: ",nrow(results)," verified-checkpoint-vs-R tables verified by keys and numeric values")
  invisible(results)
}

write_sha_manifest <- function(folder, destination) {
  files <- list.files(folder,recursive=TRUE,full.names=TRUE,all.files=FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[normalizePath(files)!=normalizePath(destination,mustWork=FALSE)]
  files <- sort(files)
  rel <- substring(files,nchar(normalizePath(folder))+2L)
  data <- data.frame(relative_path=rel,sha256=vapply(files,file_sha256,character(1)))
  write_csv(data,destination)
  invisible(data)
}
