# Regenerate candidate research figures from VERIFIED R CSVs only.
# No PCA recomputation and no private clinical metadata.

render_figures_from_outputs <- function(out_dir, figure_subdirectory = "figures") {
  if (!file.exists(file.path(out_dir,"cross_language_verification.csv"))) {
    stop("Figures blocked: missing cross-language numeric verification")
  }
  fig_dir <- file.path(out_dir, figure_subdirectory)
  dir.create(fig_dir,showWarnings=FALSE,recursive=TRUE)
  read_tab <- function(f) utils::read.csv(file.path(out_dir,f),check.names=FALSE,stringsAsFactors=FALSE)
  ev <- read_tab("adult_reference_eigenvalues_17_18_PC1_PC5.csv")
  ev <- ev[ev$model=="adult17",,drop=FALSE]
  inf <- read_tab("adult_leave_one_out_influence_aggregated.csv")
  child <- read_tab("pediatric_per_model_K3.csv")

  save_pair <- function(name, draw_fn, width=8, height=5) {
    pdf_path <- file.path(fig_dir,paste0(name,".pdf"))
    png_path <- file.path(fig_dir,paste0(name,".png"))
    grDevices::pdf(pdf_path, width=width,height=height,useDingbats=FALSE)
    tryCatch(draw_fn(),finally=grDevices::dev.off())
    grDevices::png(png_path,width=round(width*175),height=round(height*175),res=175)
    tryCatch(draw_fn(),finally=grDevices::dev.off())
  }

  save_pair("Figure1_adult_reference_variance",function() {
    graphics::par(mar=c(4.5,5.4,3,4.5))
    pct <- ev$variance_pct
    x <- seq_along(pct)
    pos <- graphics::barplot(pct,names.arg=paste0("PC",x),col="grey55",border=NA,
      ylim=c(0,55),ylab="Adult reference variance (%)",
      main="Adult-only PCA: first five components",cex.main=0.95)
    graphics::text(pos,pct+2, labels=sprintf("%.1f",pct),cex=0.8)
    graphics::mtext(sprintf("PC1–PC3: %.2f%% (n = 17)",sum(pct[1:3])),side=3,line=0,cex=0.9)
  })

  save_pair("Figure2_pediatric_Qfraction_reference_sensitivity",function() {
    keep <- c("adult17_main","adult16_no_A04","adult16_no_A07","adult15_no_A04_A07")
    ids <- as.character(child$pediatric_id[child$model==keep[1L]])
    mat <- sapply(keep,function(m) {
      v <- child[child$model==m,,drop=FALSE]
      if (!setequal(as.character(v$pediatric_id),ids)) stop("Pediatric ID drift in figure data")
      v$Q_fraction[match(ids,as.character(v$pediatric_id))]
    })
    graphics::par(mar=c(6.5,4.8,2.1,0.6))
    graphics::matplot(1:4,t(mat),type="l",lty=1,col="grey80",xaxt="n",
      xlab="",ylab="Pediatric Q fraction",ylim=c(0,1),lwd=0.7)
    graphics::axis(1,at=1:4,labels=c("n=17","n=16 (-A04)","n=16 (-A07)","n=15 (-A04/-A07)"),
                   las=2,cex.axis=0.78)
    graphics::lines(1:4,apply(mat,2,median),type="b",pch=19,lwd=2.5,col="black")
    graphics::mtext("Reference composition (adult-only fits; same 25 external spectra)",
                    side=1,line=5.4,cex=0.8)
    graphics::legend("topleft",legend=c("25 pediatric records","Median"),
                     col=c("grey75","black"),lwd=c(1,2.5),pch=c(NA,19),bty="n",cex=0.85)
  },width=8.6,height=5.6)

  save_pair("Figure3_adult_reference_LOO_angles",function() {
    graphics::par(mar=c(4.7,4.9,2.1,0.7))
    x <- seq_len(nrow(inf))
    y <- inf$subspace_max_angle_deg
    graphics::plot(x,y,type="h",lwd=2,col="grey60",ylim=c(0,max(y)*1.13),
                   xaxt="n",xlab="Omitted adult reference record",
                   ylab="Maximum principal angle (degrees)",
                   main="Leave-one-adult-out changes in the K=3 subspace")
    graphics::points(x,y,pch=19,col="grey35")
    graphics::axis(1,at=x,labels=inf$reference_row_anonymous,las=2,cex.axis=0.7)
    sel <- which(inf$reference_row_anonymous %in% c("A04","A07"))
    graphics::points(x[sel],y[sel],pch=19,cex=1.25,col="black")
    graphics::text(x[sel],y[sel],labels=round(y[sel],1),pos=3,cex=0.85)
  },width=8.6,height=5.2)
  invisible(fig_dir)
}
