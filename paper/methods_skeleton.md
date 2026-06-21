# Methods — skeleton draft (Types Without Taxa)

*Methods-paper framing; conservative wording. Pairs with results_skeleton.md.*

## Transparency and openness
All analyses were pre-registered before the held-out data were examined (OSF registration
2ekcg; the decision rule and specification grid were locked, seed 20260620). We report all
datasets, all instruments, every covariance parameterization, and the full pre-registered
decision rule; nothing was added or dropped after the hold-out was opened. Analyses were run
in R [version]; data, code, and the registration are openly available [OSF link].

## Data
We reanalyzed three large, public personality datasets, each scored to the trait (and, where
available, facet) level:

- the Johnson IPIP-NEO-120 (N = 410,376), scored to five domains and thirty facets;
- the open-source IPIP "Big Five" 50-item set (N = 603,322), scored to five domains (this
  instrument defines no facet scores);
- the IPIP-HEXACO (N = 22,734), scored to six domains and twenty-four facets.

Using three instruments lets a verdict be required to replicate across measurement systems
rather than rest on one.

**Held-out confirmation.** Each dataset was split 50/50 (seed 20260620) into an exploratory and
a confirmatory partition, sealed before any model reported here was fit. For NEO-120 and
IPIP-50 the held-out respondents are almost entirely unseen (~95–97%); for HEXACO, whose
smaller size forced overlap, ~88% of held-out respondents also appear in the exploratory split,
so HEXACO is treated as a robustness check, not a naive confirmation. All models below were fit
to random subsamples of the held-out partition (n = 8,000 at domain resolution, n = 5,000 at
facet resolution; seed 20260620). Because the number of components selected by BIC tends to
grow with sample size in finite mixtures, the real data and its matched null are always
compared at the same n.

## The covariance multiverse
We model each dataset as a Gaussian mixture and treat "types" as the components selected by
BIC. mclust [cite] supplies fourteen covariance parameterizations that constrain a component's
covariance in volume, shape, and orientation, grouped into spherical (EII, VII), diagonal (EEI,
VEI, EVI, VVI), and ellipsoidal (EEE, EVE, VEE, VVE, EEV, VEV, EVV, VVV) families. For each
parameterization we fit G = 1 to 10 components and record the BIC-selected G. The set of
fourteen selected counts for a dataset is its specification curve: the apparent number of types
as a function of a covariance assumption that carries no substantive meaning.

## The Gaussian-copula matched null
To ask whether an apparent clustering is more than its own margins and covariance produce, we
compare each dataset to a matched null that fixes both and randomizes only the higher-order
dependence. For a real n × p subsample we (i) compute the sample correlation matrix R and its
Cholesky factor L; (ii) draw an n × p matrix Z of independent standard-normal variates and
correlate it as ZL; and (iii) for each variable, replace the correlated Gaussian column with
the sorted real values of that variable placed at the ranks of the Gaussian column. The result
reproduces every univariate marginal exactly and the full covariance to within sampling error
(maximum absolute correlation difference ≈ 0.009 across datasets), while imposing Gaussian
(linear) dependence and no latent clustering. This is the precise null for the question at
hand: everything a single-cluster description captures — each margin and the whole covariance —
is held fixed, so any excess clustering in the real data reflects structure the null lacks. We
drew 80 null replicates per domain slice and 50 per facet slice.

## Taxometric analysis (CCFI)
For a direct categorical-versus-continuous test we submitted each trait's facet indicators (six
per NEO-120 trait, four per HEXACO trait) to three taxometric procedures — MAMBAC, MAXEIG, and
L-Mode — using RTaxometrics [cite]. The Comparison Curve Fit Index (CCFI) contrasts the
observed curves against curves simulated from explicitly dimensional and explicitly categorical
comparison data (10,000-case comparison populations), returning a value in [0, 1] scored on the
standard band: below .45 dimensional, .45–.55 ambiguous, above .55 taxonic [Ruscio cite]. We
report each trait's CCFI as the mean across the three procedures and retain the range across
procedures as a transparency check (Figure 2). Because the IPIP-50 instrument has no facet
scores, taxometric analysis was possible only for NEO-120 and HEXACO.

## Pre-registered decision rule
The decision rule was fixed before the hold-out was opened (full text in the registration). In
brief:

- **Real structure (Verdict 1)** is claimed for an instrument only if its real-data clustering
  gain — the BIC improvement of the best mixture over a single component, across the fourteen
  parameterizations — exceeds the 95th percentile of the matched-null gain.
- The structure is judged **continuous (no separated categories)** only if all four hold:
  (1) mean within-trait CCFI below .45 with no trait at or above .55, in both instrument
  families; (2) fewer than 50% of respondents confidently classifiable (maximum posterior
  probability above .8) under the BIC-selected solution; (3) the first split — the BIC gain from
  one to two components — accounting for more than 75% of the total gain; and (4) the real
  distribution of selected counts across the fourteen parameterizations statistically
  indistinguishable from the matched-null distribution (the real mean count lying within the
  null's 95% interval).
- A **categorical** verdict requires the converse pre-specified signatures: a CCFI at or above
  .55 replicated across instruments, a clear BIC minimum with density gaps between component
  centers, a majority of respondents confidently classifiable, and a stable selected count
  across datasets.
- NEO-120 and IPIP-50 are confirmatory and must converge; HEXACO is robustness only.

**Statistics.** ΔBIC = BIC(best G) − BIC(G = 1) on the best-of-fourteen envelope; the mean
selected count is averaged over the fourteen parameterizations; the first-split share is
[BIC(2) − BIC(1)] / [BIC(best) − BIC(1)] on that envelope; the confident-assignment share is
the proportion of respondents whose maximum posterior component probability exceeds .8 under
the BIC-selected solution.

## Software
R [version]; mclust [version] (mixtures), RTaxometrics [version] (CCFI); seed 20260620
throughout. Data, scripts, and registration: [OSF].

---

### Open drafting notes (not for the manuscript)
- Criterion (4): the registration says "real *median* selected-k within the null interval"; the
  confirmatory run used the *mean* across the fourteen models. Either switch the reported
  statistic to the median or justify the mean explicitly — do not leave them mismatched.
- Fill versions: R, mclust, RTaxometrics; add the Ruscio CCFI primary citation and the mclust
  (Scrucca et al.) citation.
- One line on why BIC (not ICL) for the confirmatory selected-count, since the registered grid
  named both.
- State n explicitly and note the n-dependence of selected-k (real vs null always matched on n)
  — this is itself a small forking-path point worth one sentence.
