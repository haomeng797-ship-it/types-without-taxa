# agnostic_kmeans_n8000.R  —  B3 (corrected scale): the k-means arm of the algorithm-agnostic
# demonstration, run at n = 8000 so it sits at the same sample size as the registered GMM result
# in Table 1 (the GMM arm is reused from real_grid_matched.rds, not recomputed). k-means is a
# non-model-based, geometric engine with a gap-statistic count selection (the same competitor
# reported in Table S5), tested against the identical copula matched null. Seed 20260620.

Sys.setenv(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")  # cap BLAS threads: vecLib SVD segfaults inside forked mclapply workers otherwise
suppressMessages({library(cluster); library(parallel)})
P2 <- Sys.getenv("P2"); if (P2=="") P2 <- "."
n  <- as.integer(Sys.getenv("N","8000"))
R  <- as.integer(Sys.getenv("R_NULL","40"))
Bg <- as.integer(Sys.getenv("GAP_B","10"))
NS <- as.integer(Sys.getenv("NSTART","5"))
NC <- as.integer(Sys.getenv("NC", max(1L, detectCores()-1L)))
SEED <- 20260620L
datasets <- c(NEO120="conf_NEO120_domain.csv", IPIP50="conf_IPIP50_domain.csv", HEXACO="conf_HEXACO_domain.csv")

gen_copula_null <- function(dd){ nn<-nrow(dd); pp<-ncol(dd)
  Cn<-cor(dd); L<-tryCatch(chol(Cn),error=function(e)chol(Cn+diag(1e-6,pp)))
  Z<-matrix(rnorm(nn*pp),nn,pp)%*%L
  m<-matrix(0,nn,pp); for(j in 1:pp){ s<-sort(dd[,j]); m[,j]<-s[rank(Z[,j],ties.method="first")] }; m }

eng_kmeans <- function(x){
  g<-tryCatch(clusGap(x,FUNcluster=kmeans,K.max=10,B=Bg,spaceH0="scaledPCA",nstart=NS,iter.max=30,verbose=FALSE),
              error=function(e)NULL)
  if(is.null(g)) return(NA_integer_)
  maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"],method="Tibs2001SEmax") }

test_one <- function(X){
  set.seed(SEED); real <- eng_kmeans(X)
  nm <- unlist(mclapply(1:R, function(i){ set.seed(SEED+i); eng_kmeans(gen_copula_null(X)) }, mc.cores=NC))
  q <- quantile(nm, c(.025,.5,.975), na.rm=TRUE)
  data.frame(real=real, lo=unname(q[1]), med=unname(q[2]), hi=unname(q[3]),
             p_exceed=(1+sum(nm>=real,na.rm=TRUE))/(sum(!is.na(nm))+1), within=real<=q[3]) }

OUT <- file.path(P2,"out","agnostic_kmeans_n8000.rds")
acc <- if (file.exists(OUT)) readRDS(OUT) else data.frame()
cat(sprintf("B3 k-means @ n=%d, R=%d, gap B=%d, nstart=%d, NC=%d, start %s\n", n,R,Bg,NS,NC, format(Sys.time(),"%H:%M")))
for (dn in names(datasets)){
  if (nrow(acc) && dn %in% acc$dataset) next
  d <- as.matrix(read.table(file.path(P2,"data",datasets[dn]), sep=",", header=FALSE))
  set.seed(SEED); X <- d[sample(nrow(d), n), ]
  cat(sprintf("[%s] %-7s ... ", format(Sys.time(),"%H:%M"), dn))
  res <- cbind(dataset=dn, engine="kmeans", test_one(X))
  acc <- rbind(acc, res); saveRDS(acc, OUT)
  cat(sprintf("real=%g  null[%.1f, %.1f]  p=%.2f  %s\n", res$real, res$lo, res$hi, res$p_exceed,
              ifelse(res$within,"within (no excess)","EXCEEDS")))
}
cat("\n=== B3 k-means n=8000 COMPLETE ===  (out/agnostic_kmeans_n8000.rds)\n")
print(acc[,c("dataset","real","lo","hi","p_exceed","within")], row.names=FALSE)
