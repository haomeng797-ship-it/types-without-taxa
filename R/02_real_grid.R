suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
gridk<-function(x) sapply(models,function(m) tryCatch(Mclust(x,G=1:10,modelNames=m,verbose=FALSE)$G,error=function(e)NA_integer_))
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
out<-list()
for(nm in names(slices)){
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); x<-d[sample(nrow(d),min(slices[[nm]],nrow(d))),,drop=FALSE]
  t0<-Sys.time(); k<-gridk(x); names(k)<-models; out[[nm]]<-k
  saveRDS(out,file.path(P2,"out","real_grid_matched.rds"))
  cat(sprintf("%-14s n=%d  k: %s  | mean %.2f range %d-%d | %ds\n",nm,nrow(x),
    paste(k,collapse=" "),mean(k,na.rm=TRUE),min(k,na.rm=TRUE),max(k,na.rm=TRUE),
    as.integer(difftime(Sys.time(),t0,units="secs"))))
}
cat("=== REALGRID MATCHED DONE ===\n")
