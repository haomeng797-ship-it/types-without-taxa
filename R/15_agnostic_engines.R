# agnostic_engines.R  —  B3: algorithm-agnostic demonstration.
# The matched null is defined at the level of the data, so any clustering pipeline can be tested
# against it. The main analysis demonstrates the count test with Gaussian mixtures (BIC) only.
# A reviewer will ask whether the "algorithm-agnostic" claim holds beyond that one engine. Here
# the identical test -- real count vs the copula null's count distribution, at a matched sample
# size -- is run with two further engines of different paradigm AND different count-selection rule:
#   GMM     : model-based, median selected components over the 14 covariance models (BIC)   [reference]
#   k-means : centroid/geometric, number of clusters chosen by the gap statistic
#   Ward    : agglomerative/linkage, number of clusters chosen by the average silhouette
# on the same held-out domain data. If the verdict (real count not exceeding the null) replicates
# across engines, the framework is agnostic in practice, not only in principle. Seed 20260620.

suppressMessages({library(mclust); library(cluster); library(parallel)})
P2 <- Sys.getenv("P2"); if (P2=="") P2 <- "."
n  <- as.integer(Sys.getenv("N","2000"))
R  <- as.integer(Sys.getenv("R_NULL","100"))
Bg <- as.integer(Sys.getenv("GAP_B","20"))
NC <- as.integer(Sys.getenv("NC", max(1L, detectCores()-1L)))
SEED <- 20260620L
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
datasets <- c(NEO120="conf_NEO120_domain.csv", IPIP50="conf_IPIP50_domain.csv", HEXACO="conf_HEXACO_domain.csv")

gen_copula_null <- function(dd){ nn<-nrow(dd); pp<-ncol(dd)
  Cn<-cor(dd); L<-tryCatch(chol(Cn),error=function(e)chol(Cn+diag(1e-6,pp)))
  Z<-matrix(rnorm(nn*pp),nn,pp)%*%L
  m<-matrix(0,nn,pp); for(j in 1:pp){ s<-sort(dd[,j]); m[,j]<-s[rank(Z[,j],ties.method="first")] }; m }

# --- three engines, each returns a single count on a data matrix ---
eng_gmm <- function(x){ kv<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(NA_integer_)
    B<-as.numeric(r$BIC[,1]); if(all(is.na(B))) NA_integer_ else which.max(B) })
  median(kv,na.rm=TRUE) }
eng_kmeans <- function(x){
  g<-tryCatch(clusGap(x,FUNcluster=kmeans,K.max=10,B=Bg,spaceH0="scaledPCA",nstart=10,iter.max=30,verbose=FALSE),
              error=function(e)NULL)
  if(is.null(g)) return(NA_integer_)
  maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"],method="Tibs2001SEmax") }
eng_ward <- function(x){
  d<-dist(x); hc<-hclust(d,method="ward.D2")
  sil<-sapply(2:10,function(k){ cl<-cutree(hc,k); mean(silhouette(cl,d)[,3]) })
  as.integer(which.max(sil)+1L) }   # silhouette is undefined at k=1, so this floors at 2
engines <- list(gmm=eng_gmm, kmeans=eng_kmeans, ward=eng_ward)

test_one <- function(X, fn){
  set.seed(SEED); real <- fn(X)
  nm <- unlist(mclapply(1:R, function(i){ set.seed(SEED+i); fn(gen_copula_null(X)) }, mc.cores=NC))
  q <- quantile(nm, c(.025,.5,.975), na.rm=TRUE)
  data.frame(real=real, lo=unname(q[1]), med=unname(q[2]), hi=unname(q[3]),
             p_exceed=(1+sum(nm>=real,na.rm=TRUE))/(sum(!is.na(nm))+1), within=real<=q[3]) }

OUT <- file.path(P2,"out","agnostic_engines.rds")
acc <- if (file.exists(OUT)) readRDS(OUT) else data.frame()
cat(sprintf("B3 agnostic: n=%d, R=%d, gap B=%d, NC=%d, start %s\n", n, R, Bg, NC, format(Sys.time(),"%H:%M")))
for (dn in names(datasets)){
  d <- as.matrix(read.table(file.path(P2,"data",datasets[dn]), sep=",", header=FALSE))
  set.seed(SEED); X <- d[sample(nrow(d), n), ]
  for (en in names(engines)){
    if (nrow(acc) && any(acc$dataset==dn & acc$engine==en)) next
    cat(sprintf("[%s] %-7s / %-6s ... ", format(Sys.time(),"%H:%M"), dn, en))
    res <- cbind(dataset=dn, engine=en, test_one(X, engines[[en]]))
    acc <- rbind(acc, res); saveRDS(acc, OUT)
    cat(sprintf("real=%g  null[%.1f, %.1f]  p=%.2f  %s\n", res$real, res$lo, res$hi, res$p_exceed,
                ifelse(res$within,"within (no excess)","EXCEEDS")))
  }
}
cat("\n=== B3 COMPLETE ===  (out/agnostic_engines.rds)\n")
for (dn in names(datasets)){ x<-acc[acc$dataset==dn,]
  cat(sprintf("%-7s | ", dn))
  for (en in names(engines)){ r<-x[x$engine==en,]
    cat(sprintf("%-6s %g in[%.1f,%.1f]%s  ", en, r$real, r$lo, r$hi, ifelse(r$within,"","  EXCEEDS!"))) }
  cat("\n") }
