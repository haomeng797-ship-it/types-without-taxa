# 11_baselines.R — competitor methods on the positive-control settings of 08_positive_control.R.
# Same generators, same seed, so the datasets are identical to Figure 3's. For each setting:
#   gap  — gap statistic with a uniform reference (Tibshirani, Walther, & Hastie, 2001), k-means base
#   bLRT — bootstrap likelihood-ratio test for the number of components (McLachlan, 1987), VVV
#   BIC  — the number selected by the BIC envelope over the fourteen covariance parameterizations
# Our test's verdicts for the same datasets are in figures/fig3_data.rds (from 08).
suppressMessages({library(mclust); library(cluster); library(parallel)})
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
SEED<-20260620L; n<-4000L; p<-5L; NC<-max(1L,detectCores()-1L)
models<-c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
gen_dep<-function(n,p,rho){ z<-sample(2,n,replace=TRUE); X<-matrix(rnorm(n*p),n,p)
  L1<-chol(matrix(c(1,rho,rho,1),2,2)); L2<-chol(matrix(c(1,-rho,-rho,1),2,2))
  i1<-z==1; i2<-z==2
  X[i1,1:2]<-X[i1,1:2]%*%L1; X[i1,3:4]<-X[i1,3:4]%*%L1
  X[i2,1:2]<-X[i2,1:2]%*%L2; X[i2,3:4]<-X[i2,3:4]%*%L2; X }
gen_types<-function(n,p,K,Delta){ M<-matrix(rnorm(K*p),K,p); M<-M/sqrt(rowSums(M^2))*(Delta/sqrt(2))
  z<-sample(K,n,replace=TRUE); M[z,,drop=FALSE]+matrix(rnorm(n*p),n,p) }
gap_k<-function(x){ g<-clusGap(x,FUNcluster=kmeans,K.max=10,B=50,spaceH0="scaledPCA",nstart=10,iter.max=30,verbose=FALSE)
  maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"],method="Tibs2001SEmax") }
blrt_k<-function(x,maxG=9){ o<-tryCatch(mclustBootstrapLRT(x,modelName="VVV",nboot=99,maxG=maxG,verbose=FALSE),error=function(e)NULL)
  if(is.null(o)) return(NA_integer_)
  ns<-which(o$p.value>.05); if(length(ns)) o$G[ns[1]] else maxG+1L }
bic_k<-function(x){ perm<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(c(k=NA,bb=NA))
    B<-as.numeric(r$BIC[,1]); c(k=which.max(B),bb=max(B,na.rm=TRUE)) })
  unname(perm["k",which.max(perm["bb",])]) }
settings<-rbind(data.frame(panel="A",signal=c(0,.3,.5,.7,.85)),data.frame(panel="B",signal=c(0,1,1.5,2,3)))
res<-mclapply(seq_len(nrow(settings)),function(i){
  pn<-settings$panel[i]; s<-settings$signal[i]
  set.seed(SEED); X<-if(pn=="A") gen_dep(n,p,s) else gen_types(n,p,2L,s)
  set.seed(SEED+i)
  data.frame(panel=pn,signal=s,truth=ifelse(s==0,1L,2L),gap=gap_k(X),blrt=blrt_k(X),bic=bic_k(X)) },mc.cores=min(NC,10L))
tab<-do.call(rbind,res)
saveRDS(tab,file.path(P2,"out","baselines.rds"))
print(tab,row.names=FALSE)
cat("=== BASELINES DONE ===\n")
