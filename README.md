# Types Without Taxa

### A Covariance-Matched-Null Multiverse Test of Categorical versus Continuous Personality Structure
*Do the "personality types" that clustering algorithms report reflect real categories of people, or are they artifacts of the analysis?*

A preregistered reanalysis, in R, of three large public personality datasets (~1.0M respondents).
When you run a clustering algorithm on personality data it almost always returns some number of
"types," whether or not real types exist. This study asks how to tell the difference, and answers
it with a matched null: a synthetic population built to share the real data's trait distributions
and its full correlation structure, but with no latent types by construction. If a typeless null
produces as many "types" as the real data, the types were never there.

The point is methodological. Change one routine setting, the covariance assumption, and the number
of types found slides from one to ten. Held against the matched null and a set of preregistered
criteria, none of the three instruments shows separated categories. Personality is just where I try
the test out. The test is the thing I'd want others to pick up and use.

**Author:** Miura Meng

## Read it
- Preregistration (OSF): https://osf.io/2ekcg
- Paper: [`paper/paper.pdf`](paper/paper.pdf) (Quarto source `paper.qmd`, also `paper.html`).

## What's here

| Path | Contents |
|---|---|
| `R/` | analysis scripts `00`–`07`, run in order |
| `prereg/` | the preregistered decision rule, specification grid, knowledge-of-data, and the sealed hold-out indices |
| `figures/` | the figures |
| `paper/` | the manuscript: Quarto source, plus rendered PDF and HTML |
| `data/` | `SOURCES.md`, where to obtain the public datasets (the data are not redistributed here) |

## Method, in one paragraph
Each dataset is modeled as a Gaussian mixture; "types" are the components selected by BIC across
mclust's fourteen covariance parameterizations. To ask whether an apparent clustering exceeds what
the data's own margins and covariance produce, each dataset is compared to a Gaussian-copula matched
null that fixes every univariate marginal and the full covariance and randomizes only the
higher-order dependence. A taxometric analysis (CCFI) provides an independent categorical-versus-
continuous test. A decision rule, fixed before the held-out data were opened, declares the structure
continuous only if all four preregistered criteria hold. See [`prereg/decision_rule.md`](prereg/decision_rule.md).

## Reproduce
Requires R (≥ 4.4) with `mclust` and `RTaxometrics`. The datasets are public but large and are not
redistributed here; obtain them via [`data/SOURCES.md`](data/SOURCES.md), score them to the
trait/facet level, then run, in order:

```bash
Rscript R/00_seal_holdout.R                       # 50/50 split, seed 20260620
Rscript R/01_make_spec_grid.R                     # the specification grid
Rscript R/02_real_grid.R                          # selected k across 14 covariance models, matched n
Rscript R/03_copula_null.R <slice> <n> <reps>     # matched-null distribution, per slice
Rscript R/04_ccfi.R && Rscript R/05_ccfi_extract.R  # taxometric CCFI per trait
Rscript R/06_criteria.R                           # confident-assignment + first-split shares
Rscript R/07_figures.R                            # Figures 1 and 2
```

The scripts are research code: paths and slice subsamples reflect the author's setup and may need
adjusting. Seed of record: **20260620**. The exact held-out respondents are listed in
`prereg/holdout_confirm_*.txt`.

## Data sources
All public. Each dataset belongs to its providers; cite them as the data source. See
[`data/SOURCES.md`](data/SOURCES.md).
- IPIP-NEO-120 (Johnson): 5 domains, 30 facets
- IPIP Big-Five 50-item (Open-Source Psychometrics): 5 domains
- IPIP-HEXACO: 6 domains, 24 facets

## Preregistration
Registered on OSF before the held-out data were analyzed (https://osf.io/2ekcg). The decision rule
and specification grid in `prereg/` are the registered versions; the confirmatory verdict appended
to `prereg/decision_rule.md` was computed after the rule was locked.

## Citation
If you use this material, please cite the registration:

> Miura Meng. (2026). *Types Without Taxa: A Covariance-Matched-Null Multiverse Test of Categorical
> versus Continuous Personality Structure.* OSF. https://osf.io/2ekcg

(Manuscript in preparation.)

## License
Code is released under the MIT License (see [`LICENSE`](LICENSE)). The data are not covered by that
license and remain subject to the terms of their original providers.
