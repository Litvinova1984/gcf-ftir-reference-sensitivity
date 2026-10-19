# Creates ALL analytical stages in one R execution, without using frozen Python
# matrices as inputs. Frozen Python tables are only read by R/02_verify.R.

compute_results <- function(root, out_dir) {
  d <- read_inputs(root)
  if (!identical(dim(d$pediatric_snv), c(25L,793L))) stop("Pediatric matrix size drift")
  if (!identical(dim(d$adult17), c(17L,793L))) stop("Adult model size drift")
  k_main <- 3L # exploratory working specification; 90% explanation was retrospective
  main <- fit_adult(d$adult17)
  mixed18 <- fit_adult(d$adult18)
  p_main <- project_fixed(main, d$pediatric_snv, k_main)
  p_mixed <- project_fixed(mixed18, d$pediatric_snv, k_main)

  # Explicit link between independent stages: THIS R-derived intermediate is
  # written first and used by every reference refit below; no Python intermediate.
  write_csv(data.frame(sample_id = d$ids, d$pediatric_snv, check.names = FALSE),
            file.path(out_dir, "pediatric_on_adult_grid_SNV.csv"))
  # The public input is already aligned + SNV; pre-SNV absorbance is not distributed.
  write_csv(data.frame(wavenumber_cm_1 = d$wn), file.path(out_dir, "adult_wavenumbers_exact.csv"))
  saveRDS(list(wn = d$wn, mean_adult17 = main$center, V_adult17 = main$V,
               eigenvalues_adult17 = main$eigenvalues, K_main = k_main),
          file.path(out_dir, "adult17_fixed_model_parameters.rds"))

  # Main reference outputs and baseline cohort sensitivities.
  write_csv(data.frame(wavenumber_cm_1 = d$wn, adult17_SNV_mean = main$center),
            file.path(out_dir, "adult17_PCA_center.csv"))
  write_csv(data.frame(wavenumber_cm_1 = d$wn, main$V, check.names = FALSE),
            file.path(out_dir, "adult17_PCA_loadings_PC1_PC5.csv"))
  write_csv(data.frame(pediatric_sample_id = d$ids,
                       PC1 = p_main$scores[,1], PC2 = p_main$scores[,2],
                       PC3 = p_main$scores[,3], Q_residual = p_main$Q,
                       Q_fraction = p_main$Q_fraction, T2 = p_main$T2),
            file.path(out_dir, "pediatric_external_scores_K3.csv"))
  write_csv(data.frame(pediatric_sample_id = d$ids,
                       Q_fraction_17 = p_main$Q_fraction, Q_fraction_18 = p_mixed$Q_fraction,
                       T2_17 = p_main$T2, T2_18 = p_mixed$T2),
            file.path(out_dir, "pediatric_model_17_vs_18_K3.csv"))
  ev <- do.call(rbind, lapply(list(adult17 = main, mixed18 = mixed18), function(m) {
    data.frame(component = 1:5, eigenvalue = m$eigenvalues,
               variance_pct = m$variance_pct,
               cumulative_pct = cumsum(m$variance_pct))
  }))
  ev$model <- rep(c("adult17", "mixed18"), each = 5L)
  ev <- ev[,c("model","component","eigenvalue","variance_pct","cumulative_pct")]
  write_csv(ev, file.path(out_dir, "adult_reference_eigenvalues_17_18_PC1_PC5.csv"))

  loo_main <- adult_loo_stats(d$adult17, k_main)
  write_csv(data.frame(anonymous_adult_index = loo_main$table$reference_row_anonymous,
                       Q_residual_LOO = loo_main$table$Q_residual_LOO,
                       Q_fraction_LOO = loo_main$table$Q_fraction_LOO,
                       T2_LOO = loo_main$table$T2_LOO),
            file.path(out_dir, "adult_loo_reference_K3.csv"))
  influence <- do.call(rbind, lapply(seq_len(17L), function(i) {
    child_i <- project_fixed(loo_main$models[[i]], d$pediatric_snv, k_main)
    data.frame(reference_row_anonymous = sprintf("A%02d",i),
      subspace_max_angle_deg = max_principal_angle(main,loo_main$models[[i]],k_main),
      median_child_q_fraction = median(child_i$Q_fraction),
      max_change_child_q_fraction = max(abs(child_i$Q_fraction-p_main$Q_fraction)),
      rho_child_q_fraction_vs_primary = spearman(child_i$Q_fraction,p_main$Q_fraction),
      median_child_T2 = median(child_i$T2),
      max_change_child_T2 = max(abs(child_i$T2-p_main$T2)),
      rho_child_T2_vs_primary = spearman(child_i$T2,p_main$T2))
  }))
  write_csv(influence,file.path(out_dir,"adult_leave_one_out_influence_aggregated.csv"))
  write_csv(do.call(rbind,lapply(2L:5L,function(k) {
    child_k <- project_fixed(main,d$pediatric_snv,k)
    data.frame(K=k,adult_variance_pct=sum(main$variance_pct[seq_len(k)]),
      child_qfrac_median=median(child_k$Q_fraction),
      child_qfrac_min=min(child_k$Q_fraction),
      child_qfrac_max=max(child_k$Q_fraction),
      child_t2_median=median(child_k$T2))
  })),file.path(out_dir,"component_sensitivity_K2_K5.csv"))
  pediatric_loo <- do.call(rbind,lapply(seq_len(17L),function(i) {
    z <- project_fixed(loo_main$models[[i]], d$pediatric_snv, k_main)
    data.frame(reference_row_anonymous=sprintf("A%02d",i),pediatric_id=d$ids,
               Q=z$Q,Q_fraction=z$Q_fraction,T2=z$T2)
  }))
  write_csv(pediatric_loo,file.path(out_dir,"pediatric_LOO_per_adult_K3.csv"))

  # Dedicated comparative models, INCLUDING the post-hoc joint A04/A07 deletion.
  subsets <- list(adult17_main=d$adult17,
                  adult16_no_A04=d$adult17[-4L,,drop=FALSE],
                  adult16_no_A07=d$adult17[-7L,,drop=FALSE],
                  adult15_no_A04_A07=d$adult17[-c(4L,7L),,drop=FALSE],
                  archival18_mixed_age=d$adult18)
  models <- lapply(subsets,fit_adult)
  comparison <- list(); child_values <- list(); dimensionality <- list()
  for (nm in names(models)) {
    model <- models[[nm]]
    x <- subsets[[nm]]
    child <- project_fixed(model,d$pediatric_snv,k_main)
    adult_loo <- if (identical(nm,"adult17_main")) loo_main else adult_loo_stats(x,k_main)
    la <- adult_loo$table
    comparison[[nm]] <- data.frame(
       model=nm,n_adult=nrow(x),adult_variance_PC1_PC3_pct=sum(model$variance_pct[1:3]),
       max_principal_angle_vs_adult17_deg=max_principal_angle(main,model,k_main),
       median_child_Q=median(child$Q),
       median_child_Q_fraction=median(child$Q_fraction),
       min_child_Q_fraction=min(child$Q_fraction),
       max_child_Q_fraction=max(child$Q_fraction),
       spearman_child_Q_fraction_vs_adult17=spearman(child$Q_fraction,p_main$Q_fraction),
       median_abs_change_child_Q_fraction_vs_adult17=median(abs(child$Q_fraction-p_main$Q_fraction)),
       max_abs_change_child_Q_fraction_vs_adult17=max(abs(child$Q_fraction-p_main$Q_fraction)),
       child_Q_fraction_increased_vs_adult17=sum(child$Q_fraction>p_main$Q_fraction+1e-10),
       median_child_T2=median(child$T2),
       spearman_child_T2_vs_adult17=spearman(child$T2,p_main$T2),
       adult_LOO_median_Q_fraction=median(la$Q_fraction_LOO),
       adult_LOO_max_Q_fraction=max(la$Q_fraction_LOO),
       children_above_adult_LOO_max_Q_fraction=sum(child$Q_fraction>max(la$Q_fraction_LOO)),
       adult_LOO_median_T2=median(la$T2_LOO),
       adult_LOO_max_T2=max(la$T2_LOO),
       children_above_adult_LOO_max_T2=sum(child$T2>max(la$T2_LOO)))
    child_values[[nm]] <- data.frame(model=nm,pediatric_id=d$ids,Q=child$Q,
                                     Q_fraction=child$Q_fraction,T2=child$T2)
    dimensionality[[nm]] <- do.call(rbind,lapply(2L:5L,function(k) {
      p <- project_fixed(model,d$pediatric_snv,k)
      base <- project_fixed(main,d$pediatric_snv,k)
      data.frame(model=nm,K=k,
         adult_variance_pct=sum(model$variance_pct[seq_len(k)]),
         max_angle_vs_adult17_deg=max_principal_angle(main,model,k),
         median_child_Q_fraction=median(p$Q_fraction),
         median_child_T2=median(p$T2),
         rho_Q_fraction_vs_adult17_at_same_K=spearman(p$Q_fraction,base$Q_fraction))
    }))
  }
  comparison <- do.call(rbind,comparison)
  child_values <- do.call(rbind,child_values)
  dimensionality <- do.call(rbind,dimensionality)
  rownames(comparison) <- rownames(child_values) <- rownames(dimensionality) <- NULL
  write_csv(comparison,file.path(out_dir,"model_comparison_K3.csv"))
  write_csv(child_values,file.path(out_dir,"pediatric_per_model_K3.csv"))
  write_csv(dimensionality,file.path(out_dir,"model_sensitivity_K2_K5.csv"))

  summary_tbl <- data.frame(metric=c("adult_n","pediatric_n","K_primary",
      "variance_explained_PC1_PC3_n17", "median_child_Q_fraction_n17",
      "median_adult_LOO_Q_fraction_n17", "max_subspace_angle_deg_17_vs_18",
      "child_Q_fraction_spearman_17_vs_18", "max_subspace_angle_deg_A04",
      "max_subspace_angle_deg_A04_A07"),
    value=c(17,25,k_main,sum(main$variance_pct[1:3]),median(p_main$Q_fraction),
      median(loo_main$table$Q_fraction_LOO),max_principal_angle(main,mixed18,3L),
      spearman(p_main$Q_fraction,p_mixed$Q_fraction),
      influence$subspace_max_angle_deg[4L],
      comparison$max_principal_angle_vs_adult17_deg[
        comparison$model=="adult15_no_A04_A07"]))
  write_csv(summary_tbl,file.path(out_dir,"summary_metrics.csv"))
  # Defensive check that no stage inadvertently selects/reorders pediatric records.
  if (!identical(as.character(child_values$pediatric_id[1:25]), d$ids) ||
      !identical(as.character(pediatric_loo$pediatric_id[1:25]), d$ids)) {
    stop("Pediatric join/order integrity failure")
  }
  invisible(list(summary=summary_tbl,comparison=comparison,dimensions=dimensionality))
}
