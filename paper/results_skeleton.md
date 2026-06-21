# Results — skeleton draft (Types Without Taxa)

*Methods-paper framing; conservative wording. All analyses are confirmatory and
pre-registered (OSF registration 2ekcg, decision rule locked before the hold-out was
opened; seed 20260620). Five held-out slices: NEO-120 and HEXACO at domain and facet
resolution, IPIP-50 at domain resolution. "Types" are mixture components selected by BIC
across 14 Gaussian covariance parameterizations in mclust.*

## 1. The number of "types" is set by the covariance assumption

We first ask how many types each dataset contains and find that the question has no single
answer. Across the 14 covariance parameterizations, the BIC-selected number of components
ranged from one to ten (Figure 1). At facet resolution the count spanned the entire
admissible grid in both NEO-120 and HEXACO; even at domain resolution it ran from three-to-ten
(HEXACO) and five-to-ten (NEO-120). The parameterization is a modeling convenience with no
substantive interpretation, yet it moves the apparent number of personality types across
nearly the whole range a reader might report. Richer (ellipsoidal) covariance models
consistently returned fewer components than simpler (spherical, diagonal) ones, because a
flexible covariance can absorb the continuous shape of the data that a constrained one fits
only by adding components. The number of types is not a property of respondents; it is
co-determined by a choice the analyst makes before seeing a single result.

## 2. There is real structure, but weak relative to the data's own margins and covariance

A multiverse of counts does not by itself show the structure is illusory. We compare each
dataset to a Gaussian-copula matched null that holds every univariate margin and the full
covariance matrix fixed and randomizes only the higher-order dependence. In all five slices
the real-data clustering gain (BIC improvement over a single component) exceeded the 95th
percentile of the matched-null distribution (NEO-120 domain 717 vs 485; IPIP-50 domain 1369
vs 913; HEXACO domain 816 vs 213; NEO-120 facet 4013 vs 1073; HEXACO facet 2701 vs 337). Real
personality data is therefore not perfectly reproduced by its margins and covariance alone:
some genuine higher-order structure is present everywhere. This is a statement about departure
from a Gaussian dependence model, not about categories.

## 3. The typical count is no different from what a typeless null produces

We then ask whether that real structure inflates the *number* of components beyond the matched
null. The null contains no latent types by construction, yet it produces multi-component BIC
solutions of its own; the question is whether the real count exceeds the null count. On the
pre-registered statistic — the median selected count across the fourteen parameterizations — the
real count fell inside the matched-null 95% interval in all five slices (NEO-120 domain 5 in
[4.0, 7.0]; NEO-120 facet 4 in [3.0, 6.4]; IPIP-50 domain 8 in [6.0, 8.8]; HEXACO domain 7 in
[3.5, 7.5]; HEXACO facet 4 in [1.0, 4.9]). The grey band in Figure 1 marks the null interval in
each panel; where the real median (dashed line) sits inside it, the typical number of types is
statistically indistinguishable from a distribution with no types at all.

For transparency we also computed the mean across parameterizations, which exceeds the null in
IPIP-50 and HEXACO. The two statistics diverge for a reason worth stating: the per-model count is
right-skewed, because a few covariance models drive the count to the ceiling of ten, and that
upper tail lifts the mean but not the robust median. The excess a mean detects therefore lives in
the tail of extreme parameterizations, not in the typical analysis — and the pre-registered
statistic was the median. That the verdict can turn on the choice between two summaries of the
same numbers is, fittingly, one more instance of the paper's own point.

## 4. No trait shows a taxonic signature

For a direct categorical-versus-continuous test we computed the Comparison Curve Fit Index
(CCFI) per trait, scored on the standard band (< .45 dimensional, .45–.55 ambiguous, > .55
taxonic) and averaged over the MAMBAC, MAXEIG, and L-Mode procedures. No trait reached the
taxonic threshold in either instrument family (Figure 2): across 11 traits the maximum was
HEXACO Honesty-Humility at .530. Mean CCFI was .377 for the NEO-120 facets and .453 for the
HEXACO facets. Several individual traits fell in the ambiguous band, more so in HEXACO; none
crossed into taxonic territory.

## 5. Most respondents cannot be confidently typed, and one split dominates

Two further pre-registered checks point the same way. Under the BIC-selected solution, the
share of respondents assignable to a single component with posterior probability above .8 was
below 50% in every slice (14.6% and 48.0% for NEO-120 domain and facet; 25.1% for IPIP-50;
30.0% and 46.0% for HEXACO). And the first split — the BIC gain from one to two components —
accounted for more than 75% of the total gain in four of five slices (88.7% and 89.6% for
NEO-120; 82.8% and 98.9% for HEXACO; the exception was IPIP-50 domain at 61.8%). Where types
were "found," respondents lay on a continuum across them rather than in confidently separated
groups.

## 6. Convergent verdict

Read against the pre-registered decision rule, the categorical signatures are absent everywhere:
no replicated CCFI ≥ .55, no clear BIC minimum, an unstable selected-k. On the pre-registered
statistic the typical number of types is indistinguishable from a typeless null in all five
slices. NEO-120 — the instrument with a naive hold-out — satisfies all four continuity criteria
and is unambiguously continuous; IPIP-50 and HEXACO meet the continuity criteria as well, with
only minor, non-categorical shortfalls (IPIP-50's BIC gain is a little less concentrated in its
first split; HEXACO's mean CCFI sits in the ambiguous band). We therefore find no evidence for
separated personality types in any of the three instruments. The apparent cluster solutions are
reproducible by a covariance-matched null; there is a mild real departure from that null
(Section 2), but it does not raise the type-count above the null or take any categorical form.
HEXACO is reported as a robustness check only, because its hold-out is not naive.

---

### Figure captions

**Figure 1.** The BIC-selected number of mixture components ("types") across 14 covariance
parameterizations, on held-out data (matched n), for five instrument/resolution slices. Each
point is one parameterization, coloured by covariance family and sorted by selected count. The
grey band marks the 95% interval of the matched-null median count; the dashed line is the real
median count. A real line inside the band indicates a typical count indistinguishable from data
with no latent types. *(file: phase2/out/spec_curve_null_median.png)*

**Figure 2.** CCFI by trait, averaged over MAMBAC, MAXEIG, and L-Mode (point = mean; segment =
range across the three procedures). Below .45 indicates dimensional structure, .45–.55
ambiguous, above .55 taxonic. No trait reaches the taxonic threshold.
*(file: phase2/out/ccfi_by_trait.png)*

---

### Open drafting notes (not for the manuscript)
- §3 reconciliation to write carefully: H1 (real > null on ΔBIC) and H3 (count ≈ null on NEO)
  are not in tension — there is real higher-order structure, it just does not take the form of
  extra *types*. Keep this explicit so a reader does not read §2 as evidence of categories.
- Avoid "types do not exist" anywhere. Use "no evidence for separated types beyond the matched
  null."
- Scope every claim to: these three instruments, 14 covariance models, this null, BIC/ICL.
- Possible §3 add: report criterion-2 (confident assignment) and the spec-curve together as the
  "even where the model splits, people are not separated" point.
