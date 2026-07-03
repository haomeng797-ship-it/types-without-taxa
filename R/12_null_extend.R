# 12_null_extend.R — extend the matched-null resampling of 03_copula_null.R to R = 200 per slice.
# Replicate i has always used set.seed(20260620 + i), so extending is seamless: the merged draws are
# exactly what a single R = 200 run would produce. Replicates run in parallel (one BLAS thread each;
# cap math-library threads on the command line as for 08).
suppressMessages(library(mclust)); suppressMessages(library(parallel))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
TARGET<-as.integer(Sys.getenv("R_TARGET","200")); NC<-max(1L,detectCores()-1L)
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
gen_copula_null<-function(d){ n<-nrow(d); p<-ncol(d)
  C<-cor(d); L<-tryCatch(chol(C),error=function(e)chol(C+diag(1e-6,p)))
  Z<-matrix(rnorm(n*p),n,p)%*%L
  m<-matrix(0,n,p); for(j in 1:p){ s<-sort(d[,j]); m[,j]<-s[rank(Z[,j],ties.method="first")] }; m }
fit_one<-function(x){
  perm<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(c(k=NA,bb=NA,b1=NA))
    B<-as.numeric(r$BIC[,1]); if(all(is.na(B))) return(c(k=NA,bb=NA,b1=NA))
    c(k=which.max(B), bb=max(B,na.rm=TRUE), b1=B[1]) })
  kvec<-perm["k",]
  list(kvec=kvec, median_k=median(kvec,na.rm=TRUE), mean_k=mean(kvec,na.rm=TRUE),
       dBIC=max(perm["bb",],na.rm=TRUE)-max(perm["b1",],na.rm=TRUE)) }
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
for(nm in names(slices)){
  f<-file.path(P2,"out",paste0("nullmedian_",nm,".rds")); res<-readRDS(f)
  if(res$R>=TARGET){ cat(nm,"already at R=",res$R,"\n"); next }
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); dr<-d[sample(nrow(d),min(slices[[nm]],nrow(d))),,drop=FALSE]
  idx<-(res$R+1):TARGET; t0<-Sys.time()
  new<-mclapply(idx,function(i){ set.seed(20260620+i); f1<-fit_one(gen_copula_null(dr))
    list(med=f1$median_k,mn=f1$mean_k,db=f1$dBIC,kv=f1$kvec) },mc.cores=NC)
  ok<-!sapply(new,is.null)
  res$null_median<-c(res$null_median,sapply(new[ok],`[[`,"med"))
  res$null_mean  <-c(res$null_mean,  sapply(new[ok],`[[`,"mn"))
  res$null_dBIC  <-c(res$null_dBIC,  sapply(new[ok],`[[`,"db"))
  res$kmat<-rbind(res$kmat,do.call(rbind,lapply(new[ok],`[[`,"kv")))
  res$R<-length(res$null_median)
  res$qmed<-quantile(res$null_median,c(.025,.5,.975),na.rm=TRUE)
  res$h3_med_indist<-(res$real$median_k>=res$qmed[1] && res$real$median_k<=res$qmed[3])
  saveRDS(res,f)
  cat(sprintf("%-14s extended to R=%d | null-median 95%%[%.2f,%.2f] | H1 95th dBIC=%.0f | H3=%s | %dm\n",
    nm,res$R,res$qmed[1],res$qmed[3],quantile(res$null_dBIC,.95,na.rm=TRUE),
    ifelse(res$h3_med_indist,"null-like","real>null"),as.integer(difftime(Sys.time(),t0,units="mins"))))
}
cat("=== NULL EXTENSION DONE ===\n")
