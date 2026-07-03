# 13_n_sensitivity.R — how the selected count grows with subsample size, domain slices.
# Motivates matching real data and null on n: the count is n-conditional by design.
suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
gridk<-function(x) sapply(models,function(m) tryCatch(Mclust(x,G=1:10,modelNames=m,verbose=FALSE)$G,error=function(e)NA_integer_))
plan<-list(NEO120_domain=c(2000,8000,32000),IPIP50_domain=c(2000,8000,32000),HEXACO_domain=c(2000,8000,NA))
out<-list()
for(nm in names(plan)){
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  ns<-plan[[nm]]; ns[is.na(ns)]<-nrow(d); ns<-ns[ns<=nrow(d)]
  for(n in ns){
    set.seed(20260620); x<-d[sample(nrow(d),n),,drop=FALSE]
    t0<-Sys.time(); k<-gridk(x)
    out[[paste(nm,n,sep="_")]]<-list(slice=nm,n=n,kvec=k,median=median(k,na.rm=TRUE))
    saveRDS(out,file.path(P2,"out","n_sensitivity.rds"))
    cat(sprintf("%-14s n=%6d median k=%.1f | k: %s | %ds\n",nm,n,median(k,na.rm=TRUE),
      paste(k,collapse=" "),as.integer(difftime(Sys.time(),t0,units="secs"))))
  }
}
cat("=== N SENSITIVITY DONE ===\n")
