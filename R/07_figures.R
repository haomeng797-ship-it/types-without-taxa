# Figures 1 and 2: reads the saved model outputs in out/, writes the figures to figures/.
# Run after 02_real_grid.R, 03_copula_null.R (per slice), 04_ccfi.R, 05_ccfi_extract.R.
suppressMessages({library(ggplot2); library(dplyr)})
OUT <- Sys.getenv("FIG_OUT", "out"); FIG <- "figures"; FONT <- Sys.getenv("FIG_FONT", "sans")
dir.create(FIG, showWarnings = FALSE)
ord  <- c("NEO120_domain","NEO120_facet","IPIP50_domain","HEXACO_domain","HEXACO_facet")
labs <- c(NEO120_domain="NEO-120 / domain (5)", NEO120_facet="NEO-120 / facet (30)",
          IPIP50_domain="IPIP-50 / domain (5)", HEXACO_domain="HEXACO / domain (6)",
          HEXACO_facet="HEXACO / facet (24)")
fam  <- function(m) if (m %in% c("EII","VII")) "spherical" else
                    if (m %in% c("EEI","VEI","EVI","VVI")) "diagonal" else "ellipsoidal"
cols <- c(spherical="#185FA5", diagonal="#EF9F27", ellipsoidal="#1D9E75")

## ---- Figure 1: spec curve vs typeless matched null (pre-registered median) ----
rg <- readRDS(file.path(OUT, "real_grid_matched.rds"))
df <- do.call(rbind, lapply(ord, function(s){
  k <- rg[[s]]; data.frame(slice=s, k=as.integer(k), family=sapply(names(k), fam)) }))
df <- df %>% group_by(slice) %>% arrange(k, .by_group=TRUE) %>% mutate(x=row_number()) %>% ungroup()
band <- do.call(rbind, lapply(ord, function(s){
  r <- readRDS(file.path(OUT, paste0("nullmedian_", s, ".rds"))); q <- r$qmed
  data.frame(slice=s, lo=q[1], hi=q[3], realmed=r$real$median_k, h3=r$h3_med_indist) }))
band$anno <- ifelse(band$h3, "H3 (median): null-like (PASS)", "H3 (median): real > null (FAIL)")
df$slice   <- factor(df$slice,   levels=ord, labels=labs[ord])
band$slice <- factor(band$slice, levels=ord, labels=labs[ord])
df$family  <- factor(df$family,  levels=c("spherical","diagonal","ellipsoidal"))
p1 <- ggplot(df, aes(x, k)) +
  geom_rect(data=band, inherit.aes=FALSE, aes(xmin=-Inf, xmax=Inf, ymin=lo, ymax=hi), fill="#9aa0a6", alpha=.30) +
  geom_hline(data=band, inherit.aes=FALSE, aes(yintercept=realmed), linetype="dashed", color="#333333", linewidth=.45) +
  geom_step(aes(group=slice), color="grey55", linewidth=.35) + geom_point(aes(color=family), size=2.5) +
  geom_text(data=band, inherit.aes=FALSE, aes(x=.7, y=10.6, label=anno), hjust=0, size=2.9, color="#222222") +
  facet_wrap(~slice, ncol=2) + scale_color_manual(values=cols) +
  scale_y_continuous(breaks=c(1,2,4,6,8,10), limits=c(.5,11)) + scale_x_continuous(breaks=NULL) +
  labs(title="The typical number of 'types' matches a typeless null (pre-registered median statistic)",
       x="14 covariance specifications, sorted by selected k", y="selected number of types  (k)", color=NULL) +
  theme_minimal(base_size=12, base_family=FONT) +
  theme(legend.position="top", panel.grid.minor=element_blank(), strip.text=element_text(face="bold"))
ggsave(file.path(FIG, "fig1_spec_curve_null_median.png"), p1, width=9, height=10.5, dpi=140, device=ragg::agg_png)

## ---- Figure 2: CCFI by trait (no trait reaches the taxonic threshold) ----
traits <- c("NEO_Neuroticism","NEO_Extraversion","NEO_Openness","NEO_Agreeableness","NEO_Conscientiousness",
            "HEX_HonestyHumility","HEX_Emotionality","HEX_eXtraversion","HEX_Agreeableness","HEX_Conscientiousness","HEX_Openness")
tlab   <- c("Neuroticism","Extraversion","Openness","Agreeableness","Conscientiousness",
            "Honesty-Humility","Emotionality","eXtraversion","Agreeableness","Conscientiousness","Openness")
d2 <- do.call(rbind, Map(function(t,l){
  r <- readRDS(file.path(OUT, paste0("ccfi_", t, ".rds"))); v <- c(r$CCFI.MAMBAC, r$CCFI.MAXEIG)
  data.frame(trait=t, lab=l, mean=mean(v), lo=min(v), hi=max(v),
             inst=ifelse(grepl("^NEO", t), "NEO-120", "HEXACO")) }, traits, tlab))
d2$trait <- factor(d2$trait, levels=rev(traits)); d2$inst <- factor(d2$inst, levels=c("NEO-120","HEXACO"))
p2 <- ggplot(d2, aes(mean, trait, color=inst)) +
  annotate("rect", xmin=-Inf, xmax=.45, ymin=-Inf, ymax=Inf, fill="#1D9E75", alpha=.07) +
  annotate("rect", xmin=.45, xmax=.55, ymin=-Inf, ymax=Inf, fill="#9aa0a6", alpha=.14) +
  annotate("rect", xmin=.55, xmax=Inf, ymin=-Inf, ymax=Inf, fill="#D1495B", alpha=.10) +
  geom_vline(xintercept=c(.45,.55), linetype="dotted", color="grey50") +
  geom_segment(aes(x=lo, xend=hi, yend=trait), linewidth=.7, alpha=.45) + geom_point(size=3.1) +
  facet_grid(inst ~ ., scales="free_y", space="free") +
  scale_y_discrete(labels=setNames(d2$lab, as.character(d2$trait))) +
  scale_color_manual(values=c("NEO-120"="#185FA5","HEXACO"="#EF9F27"), guide="none") +
  scale_x_continuous(limits=c(.15,.66), breaks=c(.2,.3,.45,.55)) +
  labs(title="No trait reaches the taxonic threshold",
       subtitle="Point = mean of the two registered procedures (MAMBAC, MAXEIG); segment = their range.\nThresholds: below .45 dimensional, .45-.55 ambiguous, above .55 taxonic.",
       x="Comparison Curve Fit Index (CCFI)", y=NULL) +
  theme_minimal(base_size=12, base_family=FONT) +
  theme(panel.grid.major.y=element_blank(), strip.text.y=element_text(angle=0, face="bold"),
        plot.subtitle=element_text(size=9.5, color="grey25"))
ggsave(file.path(FIG, "fig2_ccfi_by_trait.png"), p2, width=8.6, height=5.6, dpi=140, device=ragg::agg_png)
cat("wrote figures/fig1_spec_curve_null_median.png and figures/fig2_ccfi_by_trait.png\n")
