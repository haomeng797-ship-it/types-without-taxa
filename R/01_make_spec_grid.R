# 01_make_spec_grid.R  —  Day 2: freeze the specification grid.
# One row = one runnable analysis ("parallel universe"). The confirmatory engine
# reads spec_grid.csv and runs every row, recording the selected number of types k.

DEST <- "prereg"

cov <- c("EII","VII",                         # spherical
         "EEI","VEI","EVI","VVI",             # diagonal
         "EEE","EVE","VEE","VVE","EEV","VEV","EVV","VVV")  # ellipsoidal

g <- expand.grid(
  instrument       = c("NEO120","IPIP50","HEXACO"),
  resolution       = c("domain","facet"),
  data             = c("real","null"),        # null = Gaussian-copula matched null
  covariance_model = cov,
  criterion        = c("BIC","ICL"),
  stringsAsFactors = FALSE)

g$role <- ifelse(g$instrument == "HEXACO", "robustness", "confirmatory")  # HEXACO hold-out not naive
g <- g[order(g$instrument, g$resolution, g$data, g$criterion, match(g$covariance_model, cov)), ]
g$spec_id <- sprintf("S%03d", seq_len(nrow(g)))
g <- g[, c("spec_id","instrument","role","resolution","data","covariance_model","criterion")]

write.csv(g, file.path(DEST, "spec_grid.csv"), row.names = FALSE)

writeLines(c(
 "SPEC_GRID — column meanings (plain):",
 "  spec_id          : id for this one run",
 "  instrument       : which questionnaire (NEO120 / IPIP50 / HEXACO) — run separately, never pooled",
 "  role             : confirmatory (NEO120, IPIP50) or robustness (HEXACO; its hold-out is not naive)",
 "  resolution       : zoom level — domain (5-6 broad traits) or facet (24-30 fine facets)",
 "  data             : real (held-out people) or null (matched 'no-types' twin data)",
 "  covariance_model : the assumed CLUSTER SHAPE (mclust's 14 codes) — we sweep all 14 on purpose",
 "  criterion        : the rule for how-many-clusters (BIC or ICL); both tried because they can disagree",
 "",
 "Per-row output (added when the engine runs): selected k = the number of types that run picked.",
 paste("Total rows:", nrow(g))),
 file.path(DEST, "SPEC_GRID_README.txt"))

cat("wrote spec_grid.csv with", nrow(g), "rows\n")
