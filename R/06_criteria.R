# 06_criteria.R — criteria 2 & 3 (confident assignment, first-split share) for the selected solution.
# Fits each covariance model separately, in the same order and with the same seed as 02_real_grid.R,
# so the fits (and mclust's initialization stream) are identical to the grid's and the selected
# (model, G) agrees with real_grid_matched.rds by construction. The earlier all-models-in-one-call
# version could disagree with the grid on near-ties because mclust draws its initialization subsample
# from the RNG stream, which advances differently in the two call patterns.
suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
slices<-list(NEO120_domain=8000,NEO120_facet=5000,IPIP50_domain=8000,HEXACO_domain=8000,HEXACO_facet=5000)
grid_ref<-tryCatch(readRDS(file.path(P2,"out","real_grid_matched.rds")),error=function(e)NULL)
res<-list()
for(nm in names(slices)){
  d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
  set.seed(20260620); x<-d[sample(nrow(d),min(slices[[nm]],nrow(d))),,drop=FALSE]
  t0<-Sys.time()
  fits<-lapply(models,function(m) tryCatch(Mclust(x,G=1:10,modelNames=m,verbose=FALSE),error=function(e)NULL))
  names(fits)<-models
  k<-sapply(fits,function(f) if(is.null(f)) NA_integer_ else f$G)
  if(!is.null(grid_ref)&&!is.null(grid_ref[[nm]])){
    same<-identical(as.integer(k),as.integer(grid_ref[[nm]]))
    cat(sprintf("%-14s grid-consistency check: %s\n",nm,if(same)"OK (matches S1)" else "MISMATCH"))
    if(!same){cat("  new:",k,"\n  old:",as.integer(grid_ref[[nm]]),"\n")}
  }
  # BIC envelope over models at each G, from the same fits
  B<-sapply(fits,function(f) if(is.null(f)) rep(NA_real_,10) else {b<-as.numeric(f$BIC[,1]); length(b)<-10; b})
  env<-apply(B,1,function(r) suppressWarnings(max(r,na.rm=TRUE))); env[!is.finite(env)]<-NA
  best_m<-models[which.max(apply(B,2,function(c) suppressWarnings(max(c,na.rm=TRUE))))]
  f<-fits[[best_m]]
  conf<-if(is.null(f$z)) NA_real_ else mean(apply(f$z,1,max)>0.8)
  bic1<-env[1]; bicbest<-env[f$G]; bic2<-env[2]
  fss<-if(is.na(bic1)||is.na(bicbest)||bicbest<=bic1) NA_real_ else (bic2-bic1)/(bicbest-bic1)
  res[[nm]]<-list(model=best_m,G=f$G,conf=conf,fss=fss)
  saveRDS(res,file.path(P2,"out","criteria23.rds"))
  cat(sprintf("%-14s sel=%s G=%d | conf(>.8)=%.1f%% | first-split=%.1f%% | %ds\n",
    nm,best_m,f$G,100*conf,100*fss,as.integer(difftime(Sys.time(),t0,units="secs"))))
}
cat("=== CRITERIA (grid-consistent) DONE ===\n")
