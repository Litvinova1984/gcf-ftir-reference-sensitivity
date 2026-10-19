# Revalidate previously R-verified K=3 diagnostic exports against independent run.
verify_extension_table <- function(path,expected,keys,compare_fields=NULL,tol=1e-7) {
 a<-utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
 b<-utils::read.csv(expected,check.names=FALSE,stringsAsFactors=FALSE)
 if (anyDuplicated(do.call(paste,c(a[keys],sep="|"))) ||
     anyDuplicated(do.call(paste,c(b[keys],sep="|"))) ||
     !setequal(names(a),names(b))) stop("Extension schema/keys mismatch")
 k<-function(d) do.call(paste,c(d[keys],sep="|"))
 if (!setequal(k(a),k(b))) stop("Extension ID-set mismatch")
 a<-a[match(k(b),k(a)),,drop=FALSE]
 if (is.null(compare_fields)) compare_fields<-setdiff(names(b),keys)
 for (col in compare_fields) {
   if (is.numeric(a[[col]]) && is.numeric(b[[col]])) {
     xx<-a[[col]];yy<-b[[col]]
     if (!identical(is.na(xx),is.na(yy)) || max(abs(xx-yy),na.rm=TRUE)>tol)
       stop("Extension numerical discrepancy: ",basename(path)," : ",col)
   } else if (!identical(as.character(a[[col]]),as.character(b[[col]])))
     stop("Extension text mismatch: ",col)
 }
 invisible(TRUE)
}
verify_diagnostics_public <- function(root,diagnostic_dir) {
 e<-file.path(root,"verification","frozen_R_20261018")
 t<-file.path(diagnostic_dir,"tables")
 verify_extension_table(file.path(t,"Table_D1_Q_decomposition_K3.csv"),
       file.path(e,"Table_D1_Q_decomposition_K3.csv"),"model",tol=1e-6)
 verify_extension_table(file.path(t,"wavenumber_Q_contributions_K3.csv"),
       file.path(e,"wavenumber_Q_contributions_K3.csv"),
       c("model","wavenumber_cm_1"),tol=1e-6)
 verify_extension_table(file.path(t,"PRIVATE_pediatric_Q_decomposition_K3.csv"),
       file.path(e,"PRIVATE_pediatric_Q_decomposition_K3.csv"),
       c("model","pediatric_id"),tol=1e-6)
 message("PASS: K=3 diagnostics match verified source R tables")
 invisible(TRUE)
}
