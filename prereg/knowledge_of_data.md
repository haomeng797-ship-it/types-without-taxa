# Knowledge of Data — Study 1: do personality "types" exist?
# Day 4 of the preregistered confirmatory plan.
# The "Knowledge of Data" section of the secondary-data preregistration: an honest
# disclosure of prior contact with the data, which is what makes such a prereg credible.

## 1. Data provenance and access
All three datasets are public and were obtained from their original sources before this
preregistration:
- Johnson IPIP-NEO-120 (N = 410,376) — OSF osf.io/tbmh5 (scoring key osf.io/ycvdk);
- openpsychometrics IPIP-50 Big-Five markers (N = 603,322, after removing duplicate network submissions);
- openpsychometrics IPIP-HEXACO (N = 22,734).
[Obtained and scored prior to this registration; exact download dates in the project log.]

## 2. Prior analyses by the author (exploratory pass — already completed)
Before this preregistration I ran a complete exploratory analysis on all three datasets:
Gaussian mixture models (full covariance, BIC over k = 1..10, on 20,000-respondent subsamples);
a Gaussian-copula matched null (Cholesky covariance + rank-mapped marginals); a pre-calibrated
taxometric CCFI (MAMBAC with comparison data); and separation diagnostics (first-split share,
confident-assignment rate, density gaps). Run at domain and facet level, replicated across the
three datasets and two software implementations (Python/scikit-learn and base R/webR).

## 3. Prior RESULTS known to me (full disclosure)
I am not blind to the exploratory results. They were: the real data exceeded the matched null
(delta-BIC; p ~ .002-.007 across datasets); but the structure read as continuous, not categorical
-- ~75-89% of the BIC gain came from the first split, only ~16-23% of respondents were confidently
classifiable, domain-level CCFI was ~.48-.49 and within-trait facet-level CCFI averaged ~.42 (two
traits, Big-Five Openness ~.52 and HEXACO Honesty-Humility ~.53, just inside the ambiguous band,
none >= .55), and the selected number of components swung between 2 and 10 with covariance choices.
The HEXACO mixture reproduced Gerlach et al.'s k = 4. My confirmatory expectation is therefore
"continuous"; the decision rule and its reverse gate were written so the data could still overturn it.

## 4. Additional analyses since the working note
After the note I also ran, on the same (exploratory) data: mclust's 14 parsimonious covariance
models with BIC/ICL (selected k = 7 / 5 / 10 for NEO-120 / HEXACO / IPIP-50), and a standard
RTaxometrics CCFI on Openness facets (~.30, dimensional, on a reduced n = 1,500 run). Disclosed as
exploratory; their confirmatory versions will be run only on the hold-out.

## 5. Prior use by others
These datasets are well studied. The Johnson IPIP-NEO-120 and related large IPIP samples were
analyzed by Gerlach et al. (2018) and in the subsequent exchange (Freudenstein et al., 2019;
Katahira et al., 2020). I take the published characterizations of these data as background.

## 6. Bias risk and safeguards
Because I have already seen the exploratory results, plain preregistration is weaker than usual.
Three structural safeguards address this:
(a) Held-out confirmation set: after exploration, each dataset was partitioned with a recorded seed
    (20260620) into a 50% exploration set and a 50% held-out confirmation set (indices frozen in
    holdout_confirm_*.txt); confirmatory specifications run only on the held-out set. The exploratory
    pass used ~20,000-row subsamples drawn before this split, so for NEO-120 and IPIP-50 the held-out
    set is ~95-97% unseen, but for HEXACO ~88% was already seen -- HEXACO is therefore robustness, not
    naive confirmation.
(b) Self-calibrating contrast: the confirmatory claim is the REAL-vs-MATCHED-NULL contrast, not a raw
    value; prior knowledge of the real data cannot bias that contrast, since the same pipeline is
    applied to null draws sharing the real data's marginals and covariance.
(c) The specification grid (spec_grid.csv) and decision rule (decision_rule.md) were fixed before the
    confirmatory analysis of the held-out data.
