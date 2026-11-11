# Same 18-row and pediatric normalized public inputs as the main computation.
# Mathematical kernels and numerical tests follow the verified 2026-10-18 R diagnostic module.
run_diagnostics_public <- function(root,run,out) {
  source(file.path(root,"R","00_core.R"),local=TRUE)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
for (sub in c('tables','figures','logs')) dir.create(file.path(out,sub), recursive=TRUE, showWarnings=FALSE)

d <- read_inputs(root)
adult18 <- d$adult18
adult17 <- d$adult17
ids <- d$ids
Y <- d$pediatric_snv
wn <- d$wn
ped_df <- utils::read.csv(file.path(run,"pediatric_on_adult_grid_SNV.csv"),check.names=FALSE)
if (!identical(as.character(ped_df$sample_id),ids) ||
    max(abs(Y-as.matrix(ped_df[,-1,drop=FALSE])))>1e-9)
  stop("Public pediatric inputs differ from current run outputs")
wcheck <- utils::read.csv(file.path(run,"adult_wavenumbers_exact.csv"))[[1]]
if (max(abs(wn-wcheck))>1e-9) stop("Full-precision adult grid changed")

fit_sets <- list(adult17_main=adult17,
                 adult16_no_A04=adult17[-4L,,drop=FALSE],
                 adult16_no_A07=adult17[-7L,,drop=FALSE],
                 adult15_no_A04_A07=adult17[-c(4L,7L),,drop=FALSE],
                 archival18_mixed_age=adult18)
models <- lapply(fit_sets,fit_adult,max_k=5L)
# Assert that these are the SAME models used in the independently verified R run.
center_csv <- read.csv(file.path(run,'adult17_PCA_center.csv'))
load_csv <- read.csv(file.path(run,'adult17_PCA_loadings_PC1_PC5.csv'))
if (max(abs(models[[1]]$center-center_csv$adult17_SNV_mean))>1e-8 ||
    max(abs(models[[1]]$V-as.matrix(load_csv[,-1,drop=FALSE])))>1e-8)
  stop('Primary reference model does not match verified R basis')

reference <- read.csv(file.path(run,'pediatric_per_model_K3.csv'))
reference$pediatric_id <- as.character(reference$pediatric_id)
if (nrow(reference)!=125L || anyDuplicated(paste(reference$model,reference$pediatric_id)))
  stop('Stored R model x child reference table invalid')
K <- 3L
fingerprint <- wn>=870 & wn<=1800
if (sum(fingerprint)!=482L || sum(!fingerprint)!=311L)
  stop('Unexpected exact spectral window counts')
Q_of <- function(x,center,V,K=3L) {
  z <- sweep(x,2L,center,'-')
  VV <- V[,seq_len(K),drop=FALSE]
  e <- z - (z %*% VV) %*% t(VV)
  rowSums(e^2)
}

all_spectra <- list(); all_models <- list(); all_children <- list()
base <- models[[1]]
q00 <- project_fixed(base,Y,K)$Q
qf00 <- project_fixed(base,Y,K)$Q_fraction
base_mean_r2 <- NULL
verification_max <- 0
for (key in names(models)) {
  m <- models[[key]]
  proj <- project_fixed(m,Y,K)
  idx <- match(paste(key,ids),paste(reference$model,reference$pediatric_id))
  if (anyNA(idx)) stop('Missing model-ID pairing in frozen verified outputs: ',key)
  required <- reference[idx,c('Q','Q_fraction','T2')]
  diffs <- c(max(abs(proj$Q-required$Q)),
             max(abs(proj$Q_fraction-required$Q_fraction)),
             max(abs(proj$T2-required$T2)))
  verification_max <- max(verification_max, diffs)
  if (any(diffs>c(1e-7,1e-10,1e-8))) {
    stop('Reconstructed model fails frozen R output comparison: ',key,
         ': ',paste(format(diffs,scientific=TRUE),collapse=', '))
  }
  z <- sweep(Y,2L,m$center,'-')
  VV <- m$V[,seq_len(K),drop=FALSE]
  e <- z - (z %*% VV)%*%t(VV)
  sq <- e^2
  if (max(abs(rowSums(sq)-proj$Q))>1e-8) stop('Element-wise residuals fail Q summation: ',key)
  mean_r2 <- colMeans(sq)
  if (is.null(base_mean_r2)) base_mean_r2 <- mean_r2
  all_spectra[[key]] <- data.frame(
     model=key,wavenumber_cm_1=wn,
     window=ifelse(fingerprint,'fingerprint','high'),
     mean_r_squared=mean_r2,
     share_of_total_mean_Q=mean_r2/mean(proj$Q),
     delta_mean_r_squared_vs_primary=mean_r2-base_mean_r2)
  q10 <- Q_of(Y,m$center,base$V,K)
  q01 <- Q_of(Y,base$center,m$V,K)
  q11 <- proj$Q
  c_only <- q10-q00
  v_only <- q01-q00
  inter <- q11-q10-q01+q00
  if (max(abs(c_only+v_only+inter-(q11-q00)))>1e-8)
    stop('Center/subspace/interaction decomposition does not add: ',key)
  all_children[[key]] <- data.frame(model=key,pediatric_id=ids,
      Q_baseline=q00,Q_target=q11,Q_center_only=q10,Q_subspace_only=q01,
      delta_Q_total=q11-q00,delta_Q_center_only=c_only,
      delta_Q_basis_only=v_only,delta_Q_interaction=inter,
      Q_fraction_baseline=qf00,Q_fraction_target=proj$Q_fraction)
  all_models[[key]] <- data.frame(model=key,n_adult=m$n,K=K,
      mean_Q=mean(q11),median_Q=median(q11),
      median_Q_fraction=median(proj$Q_fraction),
      mean_total_delta_Q=mean(q11-q00),
      mean_center_only_delta_Q=mean(c_only),
      mean_basis_only_delta_Q=mean(v_only),
      mean_interaction_delta_Q=mean(inter),
      mean_shapley_center_contribution_Q=mean(c_only+0.5*inter),
      mean_shapley_subspace_contribution_Q=mean(v_only+0.5*inter),
      mean_spectral_Q_fingerprint=sum(mean_r2[fingerprint]),
      mean_spectral_Q_high=sum(mean_r2[!fingerprint]),
      fingerprint_fraction_of_mean_Q=sum(mean_r2[fingerprint])/mean(q11),
      max_mean_r2_wavenumber_cm_1=wn[which.max(mean_r2)],
      max_mean_delta_r2_wavenumber_cm_1=if(key=='adult17_main')NA_real_ else wn[which.max(mean_r2-base_mean_r2)],
      max_mean_delta_r2_at_wavenumber=if(key=='adult17_main')0 else max(mean_r2-base_mean_r2))
}
perw <- do.call(rbind,all_spectra)
summary <- do.call(rbind,all_models)
bychild <- do.call(rbind,all_children)
row.names(perw) <- row.names(summary) <- row.names(bychild) <- NULL
write.csv(perw,file.path(out,'tables','wavenumber_Q_contributions_K3.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'tables','Table_D1_Q_decomposition_K3.csv'),row.names=FALSE)
write.csv(bychild,file.path(out,'tables','PRIVATE_pediatric_Q_decomposition_K3.csv'),row.names=FALSE)

# Two monochrome-friendly spectral windows using the same full-precision grid.
colors <- c('#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd')
# One distinguishable line style per model (including black-and-white print).
linetypes <- c(1,2,3,4,5)
names(colors) <- names(linetypes) <- names(models)
model_labels <- c('Adult n=17 (fixed)','Without A04 (n=16)',
                  'Without A07 (n=16)','Without A04/A07 (n=15)',
                  'Historical n=18 (mixed age)')

# Figure 4 reviewer clarification: two different markers, computed from the
# EXISTING K=3 residual curves (not from q = Q_fraction).
# Circle: absolute maximum of the PRIMARY mean squared residual curve.
# Triangle + arrow: maximum POSITIVE difference (without A04 minus primary).
primary_spectral <- all_spectra[['adult17_main']]
a04_spectral <- all_spectra[['adult16_no_A04']]
if (!identical(primary_spectral$wavenumber_cm_1,
               a04_spectral$wavenumber_cm_1))
  stop('Figure 4 marker grid mismatch')
mark_base_i <- which.max(primary_spectral$mean_r_squared)
mark_change_i <- which.max(a04_spectral$mean_r_squared -
                          primary_spectral$mean_r_squared)
mark_base_wn <- primary_spectral$wavenumber_cm_1[mark_base_i]
mark_base_y <- primary_spectral$mean_r_squared[mark_base_i]
mark_change_wn <- primary_spectral$wavenumber_cm_1[mark_change_i]
mark_change_y0 <- primary_spectral$mean_r_squared[mark_change_i]
mark_change_y1 <- a04_spectral$mean_r_squared[mark_change_i]
if (abs(mark_base_wn - 991.29517) > 0.05 ||
    primary_spectral$window[mark_base_i] != 'fingerprint' ||
    abs(mark_change_wn - 2921.81358) > 0.05 ||
    primary_spectral$window[mark_change_i] != 'high' ||
    mark_change_y1 <= mark_change_y0)
  stop('Figure 4 annotation coordinates disagree with verified K=3 outputs')

plot_spectral <- function() {
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  layout(matrix(c(1,2,3),nrow=3),heights=c(3,3,0.9))
  for (wins in c('high','fingerprint')) {
    first <- TRUE
    for (ii in seq_along(models)) {
      nm <- names(models)[ii]
      ss <- all_spectra[[nm]]
      ss <- ss[ss$window==wins,]
      if (first) {
        ymax <- max(vapply(all_spectra,function(s) max(s$mean_r_squared[s$window==wins]),numeric(1)))
        panel_title <- if (wins=='high') {
          expression('3400-2800 ' * plain(cm)^{-1})
        } else {
          expression('1800-870 ' * plain(cm)^{-1})
        }
        par(mar=c(3.4,5.2,2.5,1.0))
        plot(ss$wavenumber_cm_1,ss$mean_r_squared,type='l',col=colors[nm],lty=linetypes[nm],lwd=2.0,
             xlim=range(ss$wavenumber_cm_1)[2:1],ylim=c(0,ymax*1.05),
             ylab=expression('Mean squared residual (' * plain(SNV)^2 * ')'),
             xlab='',cex.lab=0.91,
             main=panel_title)
        grid(nx=NA,ny=NULL,col='gray90')
        first <- FALSE
      } else lines(ss$wavenumber_cm_1,ss$mean_r_squared,col=colors[nm],lty=linetypes[nm],lwd=2.0)
    }
    if (wins=='high') {
      # The arrow compares TWO model curves at the wavenumber of maximum delta.
      graphics::arrows(mark_change_wn,mark_change_y0,
                       mark_change_wn,mark_change_y1,
                       code=2,length=0.10,col='black',lwd=1.25)
      graphics::points(mark_change_wn,mark_change_y1,pch=24,
                       bg='white',col='black',cex=1.22,lwd=1.25)
      graphics::text(mark_change_wn,mark_change_y1,
                     labels=bquote(.(sprintf('%.2f',mark_change_wn))~plain(cm)^{-1}),
                     pos=2,offset=0.6,cex=0.82,col='black')
    } else {
      # This circle is an ABSOLUTE primary-curve maximum, not a model difference.
      graphics::abline(v=mark_base_wn,lty=3,col='grey45',lwd=0.85)
      graphics::points(mark_base_wn,mark_base_y,pch=21,
                       bg='white',col='black',cex=1.22,lwd=1.25)
      graphics::text(mark_base_wn,mark_base_y,
                     labels=bquote(.(sprintf('%.2f',mark_base_wn))~plain(cm)^{-1}),
                     pos=2,offset=0.6,cex=0.82,col='black')
    }
    mtext(expression('Wavenumber (' * plain(cm)^{-1} * '), descending'),
          side=1,line=2.3,cex=0.91)
  }
  par(mar=c(0,0,0,0));plot.new()
  legend('center',legend=model_labels,col=colors,lty=linetypes,lwd=2.2,
         ncol=3,bty='n',cex=0.92)
}

plot_decomp <- function() {
  op <- par(no.readonly=TRUE);on.exit(par(op))
  sub <- summary[-1L,,drop=FALSE]
  names_parts <- c('Center only','Subspace only','Interaction')
  bars <- cbind(sub$mean_center_only_delta_Q,sub$mean_basis_only_delta_Q,sub$mean_interaction_delta_Q)
  compcolors <- c('#1f77b4','#ff7f0e','#2ca02c')
  par(mar=c(5.0,10.7,4.0,1.2))
  yy <- rev(seq_len(nrow(sub)))
  xlim <- c(min(-20,sum(pmin(bars,0),na.rm=TRUE)*1.2),max(400,max(rowSums(pmax(bars,0)))*1.3))
  plot(NA,xlim=xlim,ylim=c(0.5,nrow(sub)+0.6),axes=FALSE,
       xlab=expression('Change in mean Q, 25 fixed pediatric spectra (' * plain(SNV)^2 * ')'),
       ylab='')
  abline(v=0,col='gray35')
  axis(1);axis(2,at=yy,labels=model_labels[-1],las=1,cex.axis=.85)
  for (i in seq_len(nrow(sub))) {
    pos <- 0; neg <- 0
    for (j in 1:3) {
      v <- bars[i,j]
      if (v>=0) {rect(pos,yy[i]-.28,pos+v,yy[i]+.28,col=compcolors[j],border=NA);pos <- pos+v}
      else {rect(neg+v,yy[i]-.28,neg,yy[i]+.28,col=compcolors[j],border=NA);neg <- neg+v}
    }
    text(pos+6,yy[i],labels=sprintf('Total %+0.2f',sub$mean_total_delta_Q[i]),pos=4,cex=.8)
  }
  title('K=3: center vs subspace vs interaction (algebraic decomposition)',cex.main=0.95)
  legend('topright',legend=names_parts,fill=compcolors,bty='n',horiz=TRUE,cex=.8)
}

png(file.path(out,'figures','Figure_D1_spectral_Q_by_wavenumber_R.png'),width=1700,height=1060,res=155)
plot_spectral();dev.off()
pdf(file.path(out,'figures','Figure_D1_spectral_Q_by_wavenumber_R.pdf'),width=12,height=7.5)
plot_spectral();dev.off()
png(file.path(out,'figures','Figure_D2_Q_center_subspace_decomposition_R.png'),width=1700,height=850,res=155)
plot_decomp();dev.off()
pdf(file.path(out,'figures','Figure_D2_Q_center_subspace_decomposition_R.pdf'),width=12,height=6)
plot_decomp();dev.off()

bal <- max(abs(summary$mean_total_delta_Q - summary$mean_center_only_delta_Q -
               summary$mean_basis_only_delta_Q - summary$mean_interaction_delta_Q))
if (bal>1e-8) stop('Aggregate decomposition balance failed')
writeLines(c('PASS: 125 previously verified individual Q, q, and T2 matched.',
             sprintf('Max R-vs-R per-child diagnostic discrepancy: %.12g',verification_max),
             sprintf('Max mean-Q additivity discrepancy: %.12g',bal),
             'This analysis uses only the PRE-EXISTING K=3 adult compositions and frozen adult SNV inputs.',
             'This R diagnostic module is a new separate check, not a rerun of historical DPT/OriginPro preprocessing.',
             'No clinical claims or inferential statistics.'),
           file.path(out,'logs','R_DIAGNOSTICS_PASS.txt'))
capture.output(sessionInfo(),file=file.path(out,'logs','R_diagnostics_sessionInfo.txt'))
cat('R diagnostic verification PASS, models=',length(models),', children=',length(ids),'\n')

  invisible(out)
}
