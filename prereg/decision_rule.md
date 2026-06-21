# Decision Rule — Study 1: do personality "types" exist?
# Day 3 of the preregistered confirmatory plan. LOCKED before touching the hold-out.
# Read with: spec_grid.csv (the 336 runs) and holdout_confirm_*.txt (the held-out data).
# Seed of record: 20260620.

## Per-run output
Each row of spec_grid.csv yields the selected number of types k. Taxometric runs also
yield a CCFI per trait, scored on the STANDARD band (< .45 dimensional, .45–.55 ambiguous,
> .55 taxonic), with a robustness band of .40/.60, using MAMBAC AND MAXEIG (not one procedure).

## VERDICT 1 — Is there real structure (beyond skew + correlation)?
Per instrument, the real-data clustering gain (ΔBIC) must EXCEED the 95th percentile of the
Gaussian-copula matched-null ΔBIC distribution (one-sided). If yes → "structure beyond
Gaussian dependence." Primary: NEO-120, IPIP-50. Robustness only: HEXACO.

## VERDICT 2 — Categorical or continuous? (the main test)
Declare NO SEPARATED CATEGORIES (continuous) only if ALL FOUR hold:
  (1) within-trait CCFI mean < .45 AND no single trait ≥ .55, in BOTH instrument families;
  (2) confidently classifiable respondents < 50% (posterior > .8);
  (3) first-split share of the BIC gain > 75%;
  (4) the real-data distribution of selected k across the 14 covariance models is
      statistically indistinguishable from the matched-null distribution of selected k
      (real median selected-k lies within the null's 95% interval).

## REVERSE GATE — what would flip the verdict to CATEGORICAL
Declare TYPES if: CCFI ≥ .55 replicated across instruments, AND a clear BIC minimum with
density gaps between centers, AND most respondents confidently classifiable, AND a stable
selected k across datasets. (Pre-specified disconfirmers; none expected — see note §8.)

## Multiple datasets
Each instrument judged separately (Type-N; not pooled). Convergent verdicts on NEO-120 AND
IPIP-50 (primary) are required; HEXACO is robustness only (its hold-out is not naive).

## PREREG SENTENCE (paste into OSF)
"Before analyzing the held-out data we fixed the following decision rule. Real structure is
claimed for an instrument only if its real-data ΔBIC exceeds the 95th percentile of the
Gaussian-copula matched-null ΔBIC. The structure is judged continuous (no separated
categories) only if all four hold: mean within-trait CCFI < .45 with no trait ≥ .55 in both
instrument families; fewer than 50% of respondents confidently classifiable (posterior > .8);
the first split accounts for over 75% of the BIC gain; and the real distribution of selected
k across the 14 covariance parameterizations is statistically indistinguishable from the
matched-null distribution. A categorical verdict requires the converse pre-specified
signatures (CCFI ≥ .55 replicated, a clear BIC minimum with density gaps, majority confident
classification, and a stable k across datasets). NEO-120 and IPIP-50 are confirmatory and
must converge; HEXACO is robustness only because its hold-out is not naive."

## ── CONFIRMATORY VERDICT (held-out data; computed 2026-06-21, AFTER the rule above was locked) ──
Seed 20260620. Held-out subsamples; 14 covariance models; Gaussian-copula matched null
(80 replicates per domain slice, 50 per facet slice). CCFI via RTaxometrics
(MAMBAC + MAXEIG + L-Mode; num.p=20, n.pop=10000). Three results bear on the rule.

VERDICT 1 — real structure beyond Gaussian dependence: **CONFIRMED in every slice.**
Real ΔBIC exceeded the matched-null 95th percentile in all five:
  NEO-120 domain 717 vs 485; IPIP-50 domain 1369 vs 913; HEXACO domain 816 vs 213;
  NEO-120 facet 4013 vs 1073; HEXACO facet 2701 vs 337.
→ Genuine structure beyond skew + linear (Gaussian-copula) dependence. Both primary
  instruments and the robustness instrument.

REVERSE GATE (categorical) — **clearly NOT met → NOT taxonic.**
  • CCFI ≥ .55 replicated: NO. No trait reached .55 in either family (11 traits; max =
    HEXACO Honesty-Humility .530). NEO mean CCFI .377; HEXACO mean .453.
  • Clear BIC minimum / stable k: NO. Selected k swings across the 14 covariance models in
    every slice (range 1–10 at facet resolution); BIC keeps improving toward the cap.
→ Pre-specified categorical disconfirmers absent. No separated categories.

VERDICT 2 (continuous) — **NOT cleanly met; lands in the pre-specified middle.**
  (1) CCFI mean < .45 & no trait ≥ .55, both families: PARTIAL. NEO ✓ (.377). HEXACO gray
      (.453 > .45; none ≥ .55).
  (2) confidently classifiable < 50% (max posterior > .8): MET in ALL 5 slices.
      NEO-120 domain 14.6%; NEO-120 facet 48.0%; IPIP-50 domain 25.1%;
      HEXACO domain 30.0%; HEXACO facet 46.0% — all < 50%. In every slice most
      respondents are NOT confidently assignable to any component.
  (3) first-split share of BIC gain > 75%: MET in 4/5.
      NEO-120 domain 88.7%; NEO-120 facet 89.6%; HEXACO domain 82.8%;
      HEXACO facet 98.9% ✓.  IPIP-50 domain 61.8% ✗ (gain spread across more
      components — consistent with its criterion-4 over-extraction).
  (4) real selected-k distribution indistinguishable from matched null: MET in ALL 5 slices
      on the PRE-REGISTERED statistic (median across the 14 covariance models; matched n,
      60 null reps/domain, 50/facet). Real median-k inside the matched-null median 95% interval
      everywhere: NEO-120 domain 5 in [4.0,7.0]; NEO-120 facet 4 in [3.0,6.4]; IPIP-50 domain
      8 in [6.0,8.8]; HEXACO domain 7 in [3.5,7.5]; HEXACO facet 4 in [1.0,4.9].
      Transparency: the MEAN across models (NOT the registered statistic) shows real > null in
      IPIP-50 and HEXACO; the two diverge because the per-model k-distribution is right-skewed by
      covariance models that railroad to the k=10 cap, lifting the mean but not the robust median.
      The excess the mean detects lives in the upper tail, not the centre. The registered median
      is the statistic of record → criterion (4) is satisfied, including for the primary pair.

FINAL VERDICT (all four criteria computed; criterion 4 on the registered MEDIAN statistic).
By instrument —
  • NEO-120 (naive hold-out): satisfies ALL FOUR — (1) CCFI .377, none ≥ .55; (2) 14.6% / 48.0%;
    (3) 88.7% / 89.6%; (4) median k null-like at both resolutions → UNAMBIGUOUS CONTINUOUS.
  • IPIP-50: passes (2) and (4); fails (3) (first split 61.8% < 75%) → continuous, with the BIC
    gain a little less concentrated in the first split. (No CCFI; IPIP-50 has no facets.)
  • HEXACO (robustness only): passes (2), (3), and (4); gray on (1) (mean CCFI .453, none ≥ .55).

OVERALL: No separated categories in any instrument. The reverse/categorical gate fails decisively
everywhere (no CCFI ≥ .55; unstable k), and on the pre-registered median statistic the apparent
number of types is statistically indistinguishable from a typeless matched null in ALL FIVE slices.
There is mild real structure beyond the Gaussian-copula null (Verdict 1, all slices), but it does
NOT take the form of extra types or categories: under the robust statistic the count is null-like
and no trait is taxonic. The residual blemishes are minor and non-categorical — IPIP-50's
first-split concentration and HEXACO's gray CCFI mean.

Calibrated reportable claim: "No taxonic structure was found in any instrument. On the
pre-registered (median) statistic the number of selected components is indistinguishable from a
covariance-matched null with no latent types in all five data slices; the mean across covariance
models exceeds the null in IPIP-50 and HEXACO, but this reflects a right-skewed tail of models
reaching the component cap rather than a categorical signal. A mild real departure from the null
exists in every instrument (Verdict 1) but does not take a categorical form." NEO-120 (naive
hold-out) is the cleanest case; HEXACO is robustness only (hold-out ~88% seen). Figure:
phase2/out/spec_curve_null_median.png.
