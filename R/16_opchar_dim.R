# opchar_dim.R  —  B1: operating characteristics vs dimension.
# Reuses the positive-control machinery of 08_positive_control.R verbatim and sweeps p.
# Question a reviewer asks: does the count test keep its power to detect dependence-only
# types as irrelevant (noise) dimensions are added? Structure lives in dims 1-4 (identical
# margins, opposite correlation orientation), dims 5..p are independent N(0,1) noise.
# Self-contained (synthetic data), checkpointed per (p, rho) cell, resumable. Seed 20260620.

Sys.setenv(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
suppressMessages(library(mclust))
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
SEED <- 20260620L
n  <- as.integer(Sys.getenv("N","4000"))
R  <- as.integer(Sys.getenv("R_NULL","40"))
NC <- as.integer(Sys.getenv("NC", max(1L, parallel::detectCores()-1L)))
p_grid   <- as.integer(strsplit(Sys.getenv("P_GRID","5,10,20"), ",")[[1]])
rho_grid <- as.numeric(strsplit(Sys.getenv("RHO_GRID","0,0.3,0.5,0.7,0.85"), ",")[[1]])
dir.create("out", showWarnings=FALSE); OUT <- "out/opchar_dim.rds"

fit_med <- function(x){
  kv <- sapply(models, function(m){
    r <- tryCatch(Mclust(x, G=1:10, modelNames=m, verbose=FALSE), error=function(e) NULL)
    if (is.null(r) || is.null(r$BIC)) return(NA_integer_)
    B <- as.numeric(r$BIC[,1]); if (all(is.na(B))) NA_integer_ else which.max(B)
  }); median(kv, na.rm=TRUE)
}
gen_copula_null <- function(d){
  n <- nrow(d); p <- ncol(d); C <- cor(d)
  L <- tryCatch(chol(C), error=function(e) chol(C + diag(1e-6, p)))
  Z <- matrix(rnorm(n*p), n, p) %*% L; m <- matrix(0, n, p)
  for (j in 1:p){ s <- sort(d[,j]); m[,j] <- s[rank(Z[,j], ties.method="first")] }; m
}
gen_dep <- function(n, p, rho){   # types only in dims 1-4; dims 5..p are noise
  z <- sample(2, n, replace=TRUE); X <- matrix(rnorm(n*p), n, p)
  L1 <- chol(matrix(c(1,rho,rho,1),2,2)); L2 <- chol(matrix(c(1,-rho,-rho,1),2,2))
  i1 <- z==1; i2 <- z==2
  X[i1,1:2] <- X[i1,1:2] %*% L1; X[i1,3:4] <- X[i1,3:4] %*% L1
  X[i2,1:2] <- X[i2,1:2] %*% L2; X[i2,3:4] <- X[i2,3:4] %*% L2; X
}
run_point <- function(X){        # real median-k vs matched-null 95% interval; nulls parallel
  set.seed(SEED); real <- fit_med(X)
  nm <- unlist(parallel::mclapply(1:R, function(i){ set.seed(SEED+i); fit_med(gen_copula_null(X)) }, mc.cores=NC))
  q <- quantile(nm, c(.025,.5,.975), na.rm=TRUE)
  data.frame(real=real, lo=unname(q[1]), med=unname(q[2]), hi=unname(q[3]), fired=real>q[3])
}

done <- if (file.exists(OUT)) readRDS(OUT) else data.frame()
cat(sprintf("B1 dimension sweep: p in {%s}, rho in {%s}, n=%d, R=%d, NC=%d, start %s\n",
            paste(p_grid,collapse=","), paste(rho_grid,collapse=","), n, R, NC, format(Sys.time(),"%H:%M")))
for (pp in p_grid) for (rr in rho_grid){
  if (nrow(done) && any(done$p==pp & done$rho==rr)) next
  set.seed(SEED); X <- gen_dep(n, pp, rr)
  res <- cbind(p=pp, rho=rr, run_point(X)); done <- rbind(done, res); saveRDS(done, OUT)
  cat(sprintf("[%s] p=%2d rho=%.2f  real=%g  null[%.1f, %.1f]  %s\n",
              format(Sys.time(),"%H:%M"), pp, rr, res$real, res$lo, res$hi,
              ifelse(res$fired,"DETECTED","quiet")))
}
cat("\n=== B1 COMPLETE ===  (out/opchar_dim.rds)\n")
for (pp in p_grid){ x <- done[done$p==pp,]; x <- x[order(x$rho),]
  cat(sprintf("p=%2d | ", pp))
  for (k in seq_len(nrow(x))) cat(sprintf("rho %.2f: %g[%.1f,%.1f]%s  ",
      x$rho[k], x$real[k], x$lo[k], x$hi[k], ifelse(x$fired[k],"*","")))
  cat("\n")
}
