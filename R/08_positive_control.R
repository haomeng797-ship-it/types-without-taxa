# 08_positive_control.R  —  power check for the matched-null count test (Figure 3).
# Self-contained: generates its own synthetic data, so it runs with no external files.
#
# Two regimes, each compared to its OWN Gaussian-copula matched null:
#   A. dependence-structure types: two components with identical N(0,1) margins and
#      opposite correlation orientation on pairs (1,2) and (3,4), so the separation
#      lives only in the joint dependence. The test SHOULD detect these.
#   B. well-separated types: K components separated by Delta SD, so the separation
#      enters the univariate margins, which the margin-preserving null then reproduces.
#      The test SHOULD stay quiet here -- its honest, stated blind spot.
#
# Writes figures/fig3_positive_control.png. Seed of record: 20260620.
# Note: compute-heavy (mclust fit across 14 models x (1 real + R nulls) x grid points).

suppressMessages({library(mclust); library(ggplot2)})
FIG <- "figures"; FONT <- Sys.getenv("FIG_FONT", "sans"); dir.create(FIG, showWarnings = FALSE)
models <- c("EII","VII","EEI","VEI","EVI","VVI","EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")
SEED <- 20260620L; n <- 4000L; R <- as.integer(Sys.getenv("FIG_R", "40")); p <- 5L

# median selected k across the 14 covariance models (same statistic as the main analysis)
fit_med <- function(x){
  kv <- sapply(models, function(m){
    r <- tryCatch(Mclust(x, G = 1:10, modelNames = m, verbose = FALSE), error = function(e) NULL)
    if (is.null(r) || is.null(r$BIC)) return(NA_integer_)
    B <- as.numeric(r$BIC[, 1]); if (all(is.na(B))) NA_integer_ else which.max(B)
  }); median(kv, na.rm = TRUE)
}

# Gaussian-copula matched null: fix every margin and the full covariance, free the rest
gen_copula_null <- function(d){
  n <- nrow(d); p <- ncol(d); C <- cor(d)
  L <- tryCatch(chol(C), error = function(e) chol(C + diag(1e-6, p)))
  Z <- matrix(rnorm(n * p), n, p) %*% L; m <- matrix(0, n, p)
  for (j in 1:p){ s <- sort(d[, j]); m[, j] <- s[rank(Z[, j], ties.method = "first")] }; m
}

# A. types that live only in the dependence (identical margins, opposite correlation sign)
gen_dep <- function(n, p, rho){
  z <- sample(2, n, replace = TRUE); X <- matrix(rnorm(n * p), n, p)
  L1 <- chol(matrix(c(1, rho, rho, 1), 2, 2)); L2 <- chol(matrix(c(1, -rho, -rho, 1), 2, 2))
  i1 <- z == 1; i2 <- z == 2
  X[i1, 1:2] <- X[i1, 1:2] %*% L1; X[i1, 3:4] <- X[i1, 3:4] %*% L1
  X[i2, 1:2] <- X[i2, 1:2] %*% L2; X[i2, 3:4] <- X[i2, 3:4] %*% L2; X
}

# B. well-separated types (separation enters the margins)
gen_types <- function(n, p, K, Delta){
  M <- matrix(rnorm(K * p), K, p); M <- M / sqrt(rowSums(M^2)) * (Delta / sqrt(2))
  z <- sample(K, n, replace = TRUE); M[z, , drop = FALSE] + matrix(rnorm(n * p), n, p)
}

# one signal level: real median-k vs the matched-null median-k 95% interval
run_point <- function(X){
  set.seed(SEED); real <- fit_med(X)
  nm <- numeric(R); for (i in 1:R){ set.seed(SEED + i); nm[i] <- fit_med(gen_copula_null(X)) }
  q <- quantile(nm, c(.025, .5, .975), na.rm = TRUE)
  data.frame(real = real, lo = q[1], hi = q[3], fired = real > q[3])
}

rho_grid   <- c(0, 0.3, 0.5, 0.7, 0.85)   # panel A: within-component correlation
delta_grid <- c(0, 1, 1.5, 2, 3)          # panel B: centroid separation, in SD

A <- do.call(rbind, lapply(rho_grid, function(rho){
  set.seed(SEED); X <- gen_dep(n, p, rho)
  cbind(panel = "A. dependence-structure types (margins identical)", signal = rho, run_point(X)) }))
B <- do.call(rbind, lapply(delta_grid, function(D){
  set.seed(SEED); X <- gen_types(n, p, 2L, D)
  cbind(panel = "B. well-separated types (separation enters margins)", signal = D, run_point(X)) }))
df <- rbind(A, B)
df$status <- ifelse(df$fired, "detected (real > null)", "null-like")
saveRDS(df, file.path(FIG, "fig3_data.rds"))  # cache: re-plot (font/labels) without re-running the sim

p3 <- ggplot(df, aes(signal, real)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#9aa0a6", alpha = .30) +
  geom_line(color = "grey55", linewidth = .4) +
  geom_point(aes(color = status), size = 3) +
  facet_wrap(~panel, scales = "free_x") +
  scale_color_manual(values = c("null-like" = "#9aa0a6", "detected (real > null)" = "#D1495B"), name = NULL) +
  scale_y_continuous(breaks = c(1, 2, 4, 6, 8, 10), limits = c(.5, 10.5)) +
  labs(title = "Positive controls: when the test fires, and when it stays quiet",
       subtitle = "Synthetic data with known types. Grey: matched-null 95% interval; red: real median-k exceeds it.",
       x = "type-signal strength  (A: within-component correlation;  B: centroid separation in SD)",
       y = "selected number of types  (median k)") +
  theme_minimal(base_size = 12, base_family = FONT) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9.5, color = "grey25"))
ggsave(file.path(FIG, "fig3_positive_control.png"), p3, width = 9, height = 4.6, dpi = 140, device = ragg::agg_png)
cat("wrote figures/fig3_positive_control.png\n")
