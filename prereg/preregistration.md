# PREREGISTRATION — Do personality "types" exist?
# A matched-null multiverse test of categorical vs. continuous personality structure
# Template: OSF "Secondary Data Preregistration" (van den Akker et al., 2021).
# Author: Miura Meng (UPenn, SMART).  Seed of record: 20260620.
# Attached, frozen files: spec_grid.csv, decision_rule.md, knowledge_of_data.md,
#   holdout_confirm_{NEO120,IPIP50,HEXACO}.txt, 00_seal_holdout.R, 01_make_spec_grid.R.

## 1. STUDY INFORMATION
Research questions
 RQ1. Does the apparent clustering in large personality datasets exceed a null that holds every
      marginal AND the full covariance fixed (i.e., is there structure beyond Gaussian dependence)?
 RQ2. If so, is that structure categorical (separated types) or continuous?
 RQ3. Does the selected number of "types" k depend on covariance-parameterization and criterion
      choices that should not matter (i.e., is the type count an analyst artifact)?
Hypotheses (directional; the author's expectation is disclosed in Knowledge of Data)
 H1. Real ΔBIC exceeds the matched-null 95th percentile (structure exists).
 H2. The structure is continuous: all four criteria in decision_rule.md (Verdict 2) indicate
     "no separated categories."
 H3. Selected k varies widely across the 14 covariance models and is statistically
     indistinguishable from the matched-null k distribution.
 Disconfirmers are pre-specified in decision_rule.md (Reverse Gate); the data can overturn H2/H3.

## 2. DATA DESCRIPTION
Three public datasets (see knowledge_of_data.md §1 for sources/keys):
 - Johnson IPIP-NEO-120, N = 410,376;  - openpsychometrics IPIP-50, N = 603,322;
 - openpsychometrics IPIP-HEXACO, N = 22,734.
Obtained and scored before this registration [exact access date: __________].

## 3. VARIABLES
Indicators: standardized trait scores. Resolutions analyzed: domain (5 Big-Five / 6 HEXACO) and
facet (30 Big-Five / 24 HEXACO). No manipulated variables (observational).
Inclusion: complete responses only (already applied to the N above; IPIP-50 de-duplicated by
network submission). Missing data: complete-case. 
Hold-out: each dataset partitioned with seed 20260620 into a 50% exploration set and a 50% held-out
confirmation set (indices in holdout_confirm_*.txt). Confirmatory analyses use ONLY the held-out set.
NEO-120 and IPIP-50 are confirmatory; HEXACO is robustness only (its hold-out is not naive, ~88% seen).

## 4. KNOWLEDGE OF DATA
See knowledge_of_data.md (full disclosure of the completed exploratory pass, the results already
known to the author, additional analyses since, prior use by others, and bias safeguards).

## 5. ANALYSES
Engine: for each row of spec_grid.csv (336 runs = instrument x role x resolution x data x
14 covariance models x {BIC, ICL}), fit a Gaussian mixture (mclust) on subsamples of the held-out
set and record the selected number of types k; run the IDENTICAL grid on Gaussian-copula
matched-null draws (Cholesky covariance + rank-mapped marginals); compute within-trait taxometric
CCFI (RTaxometrics; MAMBAC and MAXEIG; standard band .45/.55, robustness .40/.60) at facet level;
compute separation diagnostics (first-split share, posterior>.8 rate, density gaps).
Null resampling: 1,000 replicates (raised from the exploratory 200-500).
Decision: exactly as in decision_rule.md — Verdict 1 (real ΔBIC vs null 95th pct), Verdict 2 (four
criteria for "continuous"), and the Reverse Gate (signatures that would force a "categorical" verdict).
Instruments judged separately (not pooled); NEO-120 and IPIP-50 must converge.
Software: R (mclust, RTaxometrics); all seeds recorded.

## 6. STATEMENT OF INTEGRITY
The analyses, specification grid, and decision rule were fixed and registered before any
confirmatory analysis of the held-out data. Any deviation from this plan will be reported and
justified in the resulting manuscript.
