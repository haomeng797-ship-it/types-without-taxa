# Data sources

All three datasets are public and were used as scored trait/facet matrices. They are large and are
**not redistributed in this repository**; obtain them from the providers below and cite them as the
data source.

## IPIP-NEO-120 (N = 410,376): 5 domains, 30 facets
Johnson, J. A. (2014). Measuring thirty facets of the Five Factor Model with a 120-item public
domain inventory. *Journal of Research in Personality, 51,* 78–89.
Source: Open Science Framework (J. A. Johnson's IPIP-NEO data).

## IPIP Big-Five 50-item (N = 603,322): 5 domains
Goldberg's IPIP Big-Five Factor Markers (50 items). Scored to five domains; this instrument defines
no facet scores.
Source: Open-Source Psychometrics Project (https://openpsychometrics.org/_rawdata/)

## IPIP-HEXACO (N = 22,734): 6 domains, 24 facets
IPIP-HEXACO item pool (Ashton, Lee, & Goldberg). Scored to six domains and twenty-four facets.
Source: Open-Source Psychometrics Project (https://openpsychometrics.org/_rawdata/).

---
Seed of record for the 50/50 hold-out split: **20260620** (see `R/00_seal_holdout.R`; the exact
held-out row indices are in `prereg/holdout_confirm_*.txt`). NEO-120 and IPIP-50 hold-outs are
nearly fully unseen (~95–97%); the HEXACO hold-out overlaps the exploratory split (~88% seen) and is
therefore used for robustness only.
