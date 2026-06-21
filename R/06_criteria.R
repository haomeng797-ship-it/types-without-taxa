suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"/Users/menghao/Downloads/study1_repro/phase2"
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
res<-list()
for(nm in names(slices)){
  n<-slices[[nm]]
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); idx<-sample(nrow(d),min(n,nrow(d))); x<-d[idx,,drop=FALSE]
  t0<-Sys.time(); m<-Mclust(x,G=1:10,verbose=FALSE)
  conf<-if(is.null(m$z)) NA_real_ else mean(apply(m$z,1,max)>0.8)
  B<-m$BIC; env<-apply(B,1,function(r) suppressWarnings(max(r,na.rm=TRUE))); env[!is.finite(env)]<-NA
  bic1<-env[1]; bicbest<-env[m$G]; bic2<-env[2]
  fss<-if(is.na(bic1)||is.na(bicbest)||bicbest<=bic1) NA_real_ else (bic2-bic1)/(bicbest-bic1)
  res[[nm]]<-list(model=m$modelName,G=m$G,conf=conf,fss=fss)
  saveRDS(res,file.path(P2,"out","criteria23.rds"))
  cat(sprintf("%-14s sel=%s G=%d | (2) conf-assign(>.8)=%.1f%% [<50%%? %s] | (3) first-split share=%.1f%% [>75%%? %s] | %ds\n",
    nm,m$modelName,m$G,100*conf,!is.na(conf)&&conf<0.5,100*fss,!is.na(fss)&&fss>0.75,as.integer(difftime(Sys.time(),t0,units="secs"))))
}
cat("\n=== CRITERIA 2&3 DONE ===\n")
