# Demographic sensitivity on already computed per-model Q_fraction, base R only.
run_metadata_public <- function(root,run,out) {
  a <- utils::read.csv(file.path(root,"data","metadata","adult_demographics_candidate.csv"),stringsAsFactors=FALSE)
  p <- utils::read.csv(file.path(root,"data","metadata","pediatric_demographics_candidate.csv"),stringsAsFactors=FALSE)
  if (!identical(a$adult_id,sprintf("A%02d",1:18)) || a$age_years[18]!=15 ||
      !identical(as.character(a$primary_adult_member),c(rep("yes",17),"no")))
    stop("Adult demographic membership mismatch")
  if (nrow(p)!=25L || anyDuplicated(p$pediatric_id) ||
      !setequal(p$pediatric_id,1:25)) stop("Pediatric metadata IDs invalid")
  if (sum(a$sex[1:17]=="F")!=13 || sum(a$sex[1:17]=="M")!=4 ||
      sum(p$sex=="F")!=12 || sum(p$sex=="M")!=13) stop("Sex frequencies differ")
  if (sum(a$age_years[1:17] <=25)!=9 || sum(a$age_years[1:17]>=50)!=8)
    stop("Adult age-band totals differ")
  x <- utils::read.csv(file.path(run,"pediatric_per_model_K3.csv"),stringsAsFactors=FALSE)
  inf <- utils::read.csv(file.path(run,"adult_leave_one_out_influence_aggregated.csv"),stringsAsFactors=FALSE)
  models <- c("adult17_main","adult16_no_A04","adult16_no_A07",
              "adult15_no_A04_A07","archival18_mixed_age")
  if (nrow(x)!=125L || anyDuplicated(paste(x$model,x$pediatric_id)) ||
      !setequal(x$model,models) || nrow(inf)!=17L ||
      !identical(as.character(inf$reference_row_anonymous),a$adult_id[1:17]))
    stop("Verified model or adult influence coverage invalid")
  rho <- function(b,c) as.numeric(stats::cor(b,c,method="spearman"))
  entries <- lapply(models,function(nm) {
    z <- x[x$model==nm,,drop=FALSE]
    i <- match(p$pediatric_id,z$pediatric_id)
    if (anyNA(i) || nrow(z)!=25L || anyDuplicated(z$pediatric_id))
      stop("Pediatric ID alignment failed for ",nm)
    q <- z$Q_fraction[i]
    f <- stats::median(q[p$sex=="F"]); m <- stats::median(q[p$sex=="M"])
    data.frame(model=nm,n_female=12L,n_male=13L,
       median_female_q=f,median_male_q=m,delta_F_minus_M=f-m,
       spearman_age_q=rho(p$age_years,q),stringsAsFactors=FALSE)
  })
  child <- do.call(rbind,entries)
  adultage <- a$age_years[1:17]
  checks <- data.frame(metric=c("age_vs_angle","age_vs_max_child_q_change",
                                     "age_vs_angle_excluding_A04_A07"),
    observed=c(rho(adultage,inf$subspace_max_angle_deg),
               rho(adultage,inf$max_change_child_q_fraction),
               rho(adultage[-c(4,7)],inf$subspace_max_angle_deg[-c(4,7)])),
    stringsAsFactors=FALSE)
  expected_child <- utils::read.csv(file.path(root,"verification","frozen_R_20261018",
                        "AGGREGATED_K3_demographic_R.csv"),stringsAsFactors=FALSE)
  expected_adult <- utils::read.csv(file.path(root,"verification","frozen_R_20261018",
                        "adult_age_influence_correlations_R.csv"),stringsAsFactors=FALSE)
  keys <- c("median_female_q","median_male_q","delta_F_minus_M","spearman_age_q")
  if (!identical(child$model,expected_child$model) ||
      max(abs(as.matrix(child[,keys])-as.matrix(expected_child[,keys])))>1e-10)
    stop("Demographic pediatric numerical mismatch with verified R")
  if (!identical(checks$metric,expected_adult$metric) ||
      max(abs(checks$observed-expected_adult$observed))>1e-10)
    stop("Demographic adult numerical mismatch with verified R")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write_csv(child,file.path(out,"AGGREGATED_K3_demographic_R.csv"))
  write_csv(checks,file.path(out,"adult_age_influence_correlations_R.csv"))
  writeLines("PUBLIC-LEVEL METADATA DIAGNOSTICS PASS: 5 models, 25 children, 17 adults",
             file.path(out,"R_METADATA_VERIFICATION_PASS.txt"))
  invisible(out)
}
