SPEC_GRID — column meanings (plain):
  spec_id          : id for this one run
  instrument       : which questionnaire (NEO120 / IPIP50 / HEXACO) — run separately, never pooled
  role             : confirmatory (NEO120, IPIP50) or robustness (HEXACO; its hold-out is not naive)
  resolution       : zoom level — domain (5-6 broad traits) or facet (24-30 fine facets)
  data             : real (held-out people) or null (matched 'no-types' twin data)
  covariance_model : the assumed CLUSTER SHAPE (mclust's 14 codes) — we sweep all 14 on purpose
  criterion        : the rule for how-many-clusters (BIC or ICL); both tried because they can disagree

Per-row output (added when the engine runs): selected k = the number of types that run picked.
Total rows: 336
