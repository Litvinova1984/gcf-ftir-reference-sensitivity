# Core deterministic mathematics for the MDPI Data PRIVATE reproducibility workbench.
# All models are fit ONLY to adults. Pediatric observations are external projections.
# Base R only. Does not recreate the historical DPT -> OriginPro processed adult matrix.

as_num_matrix <- function(x, expected_dim = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(x) == 0L || anyNA(x) || any(!is.finite(x))) stop("Missing or nonfinite numerical matrix")
  if (!is.null(expected_dim) && !identical(dim(x), as.integer(expected_dim))) {
    stop("Unexpected matrix dimensions: got ", paste(dim(x), collapse = "x"),
         ", expected ", paste(expected_dim, collapse = "x"))
  }
  x
}

snv_rows <- function(x) {
  x <- as_num_matrix(x)
  mu <- rowMeans(x)
  sigma <- apply(x, 1L, stats::sd) # sample SD with n-1 denominator
  if (any(!is.finite(sigma) | sigma <= 0)) stop("Invalid within-spectrum SNV denominator")
  sweep(sweep(x, 1L, mu, "-"), 1L, sigma, "/")
}

read_inputs <- function(root) {
  # PUBLIC starting level: pre-normalized 18 x 793 adults and 25 x 793 children;
  # no upstream raw DPT, KBr, OriginPro, OPUS or pediatric interpolation.
  a <- utils::read.csv(file.path(root,"data","adult_SNV_18x793.csv"),
                       check.names=FALSE,stringsAsFactors=FALSE)
  p <- utils::read.csv(file.path(root,"data","pediatric_SNV_25x793.csv"),
                       check.names=FALSE,stringsAsFactors=FALSE)
  w <- utils::read.csv(file.path(root,"data","wavenumbers_exact.csv"),
                       check.names=FALSE,stringsAsFactors=FALSE)
  if (nrow(a)!=18L || ncol(a)!=794L || nrow(p)!=25L || ncol(p)!=794L ||
      nrow(w)!=793L || ncol(w)!=3L) stop("Public input dimensions invalid")
  if (!identical(as.character(a$adult_id),sprintf("A%02d",1:18)))
    stop("Adult record ordering or membership changed")
  ids <- as.character(p$pediatric_id)
  if (anyDuplicated(ids) || !setequal(ids,as.character(1:25)) ||
      !identical(ids, sort(ids))) stop("Pediatric lexical participant order changed")
  wn <- as.numeric(w$wavenumber_cm_1)
  if (length(wn)!=793L || any(!is.finite(wn)) || anyDuplicated(wn) ||
      any(diff(wn)>=0)) stop("Invalid exact adult wavenumber vector")
  expected <- sprintf("wn_%.5f",wn)
  if (!identical(as.character(w$spectral_column),expected) ||
      !identical(names(a)[-1],expected) || !identical(names(p)[-1],expected) ||
      !identical(as.integer(w$variable_index_1based),1:793))
    stop("Spectral grid/variable column mismatch")
  adult18 <- as_num_matrix(a[,-1,drop=FALSE],c(18L,793L))
  Y <- as_num_matrix(p[,-1,drop=FALSE],c(25L,793L))
  if (max(abs(rowMeans(adult18)))>1e-8 ||
      max(abs(apply(adult18,1L,stats::sd)-1))>1e-8 ||
      max(abs(rowMeans(Y)))>1e-8 ||
      max(abs(apply(Y,1L,stats::sd)-1))>1e-8)
    stop("Normalized input fails row-wise SNV checks")
  # A18 is the 15-year-old historical record; never part of the primary adults.
  list(adult18=adult18,adult17=adult18[-18L,,drop=FALSE],
       wn=wn,ids=ids,pediatric_snv=Y)
}

fit_adult <- function(x, max_k = 5L) {
  x <- as_num_matrix(x)
  if (nrow(x) < 12L || ncol(x) != 793L) stop("Invalid adult reference dimensions")
  if (max_k < 2L || max_k > min(nrow(x)-1L,ncol(x))) stop("Invalid PCA rank")
  p <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = max_k)
  V <- p$rotation[, seq_len(max_k), drop = FALSE]
  for (j in seq_len(max_k)) {
    idx <- which.max(abs(V[, j]))
    if (V[idx, j] < 0) V[, j] <- -V[, j]
  }
  ev <- as.numeric(p$sdev[seq_len(max_k)]^2)
  total <- sum(sweep(x, 2L, p$center, "-")^2) / (nrow(x)-1L)
  if (any(!is.finite(ev)) || any(ev <= 0) || !is.finite(total) || total <= 0) {
    stop("Degenerate PCA fit")
  }
  if (max(abs(crossprod(V) - diag(max_k))) > 1e-8) stop("PCA vectors not orthonormal")
  list(n = nrow(x), center = as.numeric(p$center), V = V,
       eigenvalues = ev, variance_pct = 100 * ev / total)
}

project_fixed <- function(m, x, K = 3L) {
  x <- as_num_matrix(x)
  if (ncol(x) != length(m$center) || K > ncol(m$V)) stop("Projection has incompatible grid or rank")
  z <- sweep(x, 2L, m$center, "-")
  V <- m$V[, seq_len(K), drop = FALSE]
  t_score <- z %*% V
  res <- z - t_score %*% t(V)
  Q <- rowSums(res^2)
  denominator <- rowSums(z^2)
  if (any(!is.finite(denominator) | denominator <= 0)) stop("Undefined Q_fraction denominator")
  Qf <- Q / denominator
  T2 <- rowSums(sweep(t_score^2, 2L, m$eigenvalues[seq_len(K)], "/"))
  if (any(Qf < -1e-10 | Qf > 1 + 1e-10) || any(T2 < 0)) stop("Invalid projection diagnostics")
  list(scores = t_score, Q = Q, Q_fraction = Qf, T2 = T2)
}

max_principal_angle <- function(m1, m2, K = 3L) {
  if (K > min(ncol(m1$V), ncol(m2$V))) stop("Invalid subspace comparison dimension")
  s <- svd(t(m1$V[, seq_len(K), drop = FALSE]) %*%
             m2$V[, seq_len(K), drop = FALSE], nu = 0L, nv = 0L)$d
  max(acos(pmax(0, pmin(1, s)))) * 180 / pi
}

adult_loo_stats <- function(x, K = 3L) {
  output <- data.frame(reference_row_anonymous = sprintf("A%02d", seq_len(nrow(x))),
                      Q_residual_LOO = numeric(nrow(x)),
                      Q_fraction_LOO = numeric(nrow(x)), T2_LOO = numeric(nrow(x)))
  models <- vector("list", nrow(x))
  for (i in seq_len(nrow(x))) {
    m <- fit_adult(x[-i, , drop = FALSE])
    models[[i]] <- m
    v <- project_fixed(m, x[i, , drop = FALSE], K)
    output$Q_residual_LOO[i] <- v$Q[1]
    output$Q_fraction_LOO[i] <- v$Q_fraction[1]
    output$T2_LOO[i] <- v$T2[1]
  }
  list(table = output, models = models)
}

spearman <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
