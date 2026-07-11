# calib_tcopula.R  —  B2: false-positive calibration under non-Gaussian DEPENDENCE.
# Continuation of 14_fpr_calibration.R / 14b_fpr_skew.R. Those hold the copula Gaussian and
# only bend the margins (Bauer-Curran skew). This asks the sharper question: when the typeless
# continuum's DEPENDENCE is non-Gaussian (a t-copula, with tail dependence), does the count
# test H3 falsely fire? The matched null preserves margins + linear correlation and Gaussianises
# the dependence, so any excess here isolates the effect of the copula, not the margins.
# Margins are held standard normal so the copula is the only thing changed. Typeless throughout.
# Correlation matrix taken from the real NEO-120 domains, matching 14b. Seed 20260620.

suppressMessages({library(mclust); library(cluster); library(parallel)})
P2 <- Sys.getenv("P2"); if (P2=="") P2 <- "."
M  <- as.integer(Sys.getenv("M","100"))
R  <- as.integer(Sys.getenv("R_NULL","40"))
n  <- as.integer(Sys.getenv("N","2000"))
NC <- as.integer(Sys.getenv("NC", max(1L, detectCores()-1L)))
df_grid <- as.numeric(strsplit(Sys.getenv("DF_GRID","3,8"), ",")[[1]])
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")

d  <- as.matrix(read.table(file.path(P2,"data","conf_NEO120_domain.csv"), sep=",", header=FALSE))
set.seed(20260620); C <- cor(d[sample(nrow(d),8000),]); Lc <- chol(C); p <- ncol(C)

gen_copula_null <- function(dd){ nn<-nrow(dd); pp<-ncol(dd)
  Cn<-cor(dd); L<-tryCatch(chol(Cn),error=function(e)chol(Cn+diag(1e-6,pp)))
  Z<-matrix(rnorm(nn*pp),nn,pp)%*%L
  m<-matrix(0,nn,pp); for(j in 1:pp){ s<-sort(dd[,j]); m[,j]<-s[rank(Z[,j],ties.method="first")] }; m }

# typeless data with a t-copula: correlation C, df tail dependence, standard-normal margins
gen_tcop <- function(n, Lc, df){
  Z  <- matrix(rnorm(n*p), n, p) %*% Lc      # correlated normal, corr = C
  W  <- rchisq(n, df)/df                     # one scaling per row -> multivariate t
  Tm <- Z / sqrt(W)                          # rows are multivariate-t(scale C, df)
  qnorm(pt(Tm, df=df))                        # t-copula, normal margins (single population)
}

kmed <- function(x){ kv<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(NA_integer_)
    B<-as.numeric(r$BIC[,1]); if(all(is.na(B))) NA_integer_ else which.max(B) })
  median(kv,na.rm=TRUE) }
bic_k <- function(x){ perm<-sapply(models,function(mm){
    r<-tryCatch(Mclust(x,G=1:10,modelNames=mm,verbose=FALSE),error=function(e)NULL)
    if(is.null(r)||is.null(r$BIC)) return(c(k=NA,bb=NA))
    B<-as.numeric(r$BIC[,1]); c(k=which.max(B),bb=max(B,na.rm=TRUE)) })
  unname(perm["k",which.max(perm["bb",])]) }
gap_k <- function(x){ g<-clusGap(x,FUNcluster=kmeans,K.max=10,B=50,spaceH0="scaledPCA",nstart=10,iter.max=30,verbose=FALSE)
  maxSE(g$Tab[,"gap"],g$Tab[,"SE.sim"],method="Tibs2001SEmax") }
blrt_k <- function(x,maxG=6){ o<-tryCatch(mclustBootstrapLRT(x,modelName="VVV",nboot=99,maxG=maxG,verbose=FALSE),error=function(e)NULL)
  if(is.null(o)) return(NA_integer_)
  ns<-which(o$p.value>.05); if(length(ns)) o$G[ns[1]] else maxG+1L }

OUT <- file.path(P2,"out","calib_tcopula.rds")
acc <- if (file.exists(OUT)) readRDS(OUT) else data.frame()
for (df in df_grid){
  sc <- sprintf("tcop_df%g", df)
  if (nrow(acc) && sc %in% acc$scenario){ cat(sprintf("skip %s (already done)\n", sc)); next }
  cat(sprintf("scenario %s: M=%d, NC=%d, start %s\n", sc, M, NC, format(Sys.time(),"%H:%M")))
  rows <- mclapply(1:M, function(m){
    set.seed(32000000 + round(df)*100000 + m); X <- gen_tcop(n, Lc, df)
    real <- kmed(X); kb <- bic_k(X)
    nulls <- sapply(1:R, function(r){ set.seed(32000000 + round(df)*100000 + m*1000 + r); kmed(gen_copula_null(X)) })
    thr <- quantile(nulls, .975, na.rm=TRUE)
    gp <- gap_k(X); bl <- if (m<=30) blrt_k(X) else NA_integer_
    data.frame(scenario=sc, df=df, m=m, real_med=real, null975=thr, fired=real>thr, bic_k=kb, gap_k=gp, blrt_k=bl)
  }, mc.cores=NC)
  acc <- rbind(acc, do.call(rbind, rows)); saveRDS(acc, OUT)
  o <- acc[acc$scenario==sc,]
  cat(sprintf("  %s done | H3 fired %.1f%% (nominal 2.5%%) | BIC k>1 %.0f%% (mean %.1f) | gap k>1 %.0f%% | bLRT k>1 %.0f%% of %d\n",
      sc, 100*mean(o$fired,na.rm=T), 100*mean(o$bic_k>1,na.rm=T), mean(o$bic_k,na.rm=T),
      100*mean(o$gap_k>1,na.rm=T), 100*mean(o$blrt_k>1,na.rm=T), sum(!is.na(o$blrt_k))))
}
cat("\n=== B2 COMPLETE ===  (out/calib_tcopula.rds)\n")
