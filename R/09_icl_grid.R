# 09_icl_grid.R — the registered grid's ICL arm: the number of components ICL selects under each
# covariance parameterization, on the same held-out subsamples as the BIC grid (02_real_grid.R).
# Reported as a supplementary specification curve; the registered null-referenced tests use BIC.
suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
out<-list()
for(nm in names(slices)){
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); x<-d[sample(nrow(d),min(slices[[nm]],nrow(d))),,drop=FALSE]
  t0<-Sys.time()
  k<-sapply(models,function(m) tryCatch({
    ic<-mclustICL(x,G=1:10,modelNames=m,verbose=FALSE)
    v<-as.numeric(ic[,1]); if(all(is.na(v))) NA_integer_ else which.max(v)
  },error=function(e)NA_integer_))
  names(k)<-models; out[[nm]]<-k
  saveRDS(out,file.path(P2,"out","icl_grid.rds"))
  cat(sprintf("%-14s ICL k: %s | median %.1f | %ds\n",nm,paste(k,collapse=" "),
    median(k,na.rm=TRUE),as.integer(difftime(Sys.time(),t0,units="secs"))))
}
cat("=== ICL GRID DONE ===\n")
