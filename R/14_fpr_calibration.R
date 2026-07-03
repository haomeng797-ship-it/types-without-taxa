# 14_fpr_calibration.R — false-positive calibration of the matched-null count test (H3) on typeless
# data, against the false-positive behavior of BIC selection, the gap statistic, and the bootstrap LRT.
# Two typeless scenarios, both built on the real NEO-120 domain correlation matrix:
#   gauss — jointly Gaussian, so margins and copula are both innocent
#   skew  — the same Gaussian copula with skewed margins (chi-square and beta transforms):
#           no types by construction, but non-normal, the over-extraction regime of
#           Bauer & Curran (2003). A calibrated test should stay quiet here too.
# M datasets per scenario at n = 2000; each dataset compared to R = 40 of its own copula nulls;
# the test fires when the real fourteen-model median exceeds the nulls' 97.5th percentile.
suppressMessages({library(mclust); library(cluster); library(parallel)})
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
M<-as.integer(Sys.getenv("M","100")); R<-40L; n<-2000L; NC<-max(1L,detectCores()-1L)
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
d<-as.matrix(read.table(file.path(P2,"data","conf_NEO120_domain.csv"),sep=",",header=FALSE))
set.seed(20260620); C<-cor(d[sample(nrow(d),8000),]); Lc<-chol(C); p<-ncol(C)
gen_copula_null<-function(dd){ nn<-nrow(dd); pp<-ncol(dd)
  Cn<-cor(dd); L<-tryCatch(chol(Cn),error=function(e)chol(Cn+diag(1e-6,pp)))
  Z<-matrix(rnorm(nn*pp),nn,pp)%*%L
  m<-matrix(0,nn,pp); for(j in 1:pp){ s<-sort(dd[,j]); m[,j]<-s[rank(Z[,j],ties.method="first")] }; m }
skew_margins<-function(Z){ U<-pnorm(Z)
  X<-cbind(qchisq(U[,1],df=3),qchisq(U[,2],df=5),qbeta(U[,3],2,5),qchisq(U[,4],df=4),U[,5]^2)
  scale(X) }
kmed<-function(x){ kv<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(NA_integer_)
    B<-as.numeric(r$BIC[,1]); if(all(is.na(B))) NA_integer_ else which.max(B) })
  median(kv,na.rm=TRUE) }
bic_k<-function(x){ perm<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(c(k=NA,bb=NA))
    B<-as.numeric(r$BIC[,1]); c(k=which.max(B),bb=max(B,na.rm=TRUE)) })
  unname(perm["k",which.max(perm["bb",])]) }
gap_k<-function(x){ g<-clusGap(x,FUNcluster=kmeans,K.max=10,B=50,spaceH0="scaledPCA",nstart=10,iter.max=30,verbose=FALSE)
  maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"],method="Tibs2001SEmax") }
blrt_k<-function(x,maxG=6){ o<-tryCatch(mclustBootstrapLRT(x,modelName="VVV",nboot=99,maxG=maxG,verbose=FALSE),error=function(e)NULL)
  if(is.null(o)) return(NA_integer_)
  ns<-which(o$p.value>.05); if(length(ns)) o$G[ns[1]] else maxG+1L }
run_scenario<-function(name){
  rows<-mclapply(1:M,function(m){
    set.seed(31000000+m); Z<-matrix(rnorm(n*p),n,p)%*%Lc
    X<-if(name=="gauss") Z else skew_margins(Z)
    real<-kmed(X); kb<-bic_k(X)
    nulls<-sapply(1:R,function(r){ set.seed(31000000+m*1000+r); kmed(gen_copula_null(X)) })
    thr<-quantile(nulls,.975,na.rm=TRUE)
    gp<-if(m<=M) gap_k(X) else NA
    bl<-if(m<=30) blrt_k(X) else NA_integer_
    data.frame(scenario=name,m=m,real_med=real,null975=thr,fired=real>thr,bic_k=kb,gap_k=gp,blrt_k=bl) },mc.cores=NC)
  do.call(rbind,rows) }
acc<-list()
for(sc in c("gauss","skew")){
  acc[[sc]]<-run_scenario(sc)
  saveRDS(do.call(rbind,acc),file.path(P2,"out","fpr_calibration.rds"))
  cat("checkpoint saved after scenario:",sc,"\n")
}
out<-do.call(rbind,acc)
for(s in unique(out$scenario)){ o<-out[out$scenario==s,]
  cat(sprintf("%-6s M=%d | H3 fired: %.1f%% (nominal 2.5%%) | BIC k>1: %.0f%% (mean k %.1f) | gap k>1: %.0f%% | bLRT k>1: %.0f%% of %d\n",
    s,nrow(o),100*mean(o$fired,na.rm=TRUE),100*mean(o$bic_k>1,na.rm=TRUE),mean(o$bic_k,na.rm=TRUE),
    100*mean(o$gap_k>1,na.rm=TRUE),100*mean(o$blrt_k>1,na.rm=TRUE),sum(!is.na(o$blrt_k)))) }
cat("=== FPR CALIBRATION DONE ===\n")
