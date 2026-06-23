suppressMessages(library(mclust))
P2<-Sys.getenv("P2"); if(P2=="")P2<-"."
a<-commandArgs(trailingOnly=TRUE); nm<-a[1]; n<-as.integer(a[2]); R<-as.integer(a[3])
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
d<-as.matrix(read.table(file.path(P2,"data",paste0("conf_",nm,".csv")),sep=",",header=FALSE))
set.seed(20260620); dr<-d[sample(nrow(d),min(n,nrow(d))),,drop=FALSE]
real<-fit_one(dr)
med<-numeric(R); mn<-numeric(R); db<-numeric(R); km<-matrix(NA_real_,R,14)
for(i in 1:R){ set.seed(20260620+i); z<-gen_copula_null(dr); f<-fit_one(z)
  med[i]<-f$median_k; mn[i]<-f$mean_k; db[i]<-f$dBIC; km[i,]<-f$kvec }
qm<-quantile(med,c(.025,.5,.975),na.rm=TRUE)
res<-list(slice=nm,n=nrow(dr),R=R,real=real,null_median=med,null_mean=mn,null_dBIC=db,kmat=km,qmed=qm,
          h3_med_indist=(real$median_k>=qm[1] && real$median_k<=qm[3]))
saveRDS(res,file.path(P2,"out",paste0("nullmedian_",nm,".rds")))
cat(sprintf("%-14s real median-k=%.1f  null-median 95%%[%.1f,%.1f]  H3(median)=%s\n",
  nm,real$median_k,qm[1],qm[3],ifelse(res$h3_med_indist,"null-like","real>null")))
