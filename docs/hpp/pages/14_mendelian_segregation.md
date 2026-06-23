# Bloc 14 — Per-family × per-LRR Mendelian segregation analysis

| | |
|---|---|
| **module** | `src/hpp/mendelian_segregation.py` |
| **framework steps** | 7–8 of the family-based-validation framework |
| **purpose** | for each (family, candidate LRR), classify the parental cross type and test offspring segregation against Mendelian expectation |

## Cross types (collapsed to canonical ordering)

| Cross | Expected | Informative for |
|---|---|---|
| HOM_REF × HOM_REF | 100% HOM_REF | — (uninformative) |
| HOM_REF × HET | 50% HOM_REF / 50% HET | segregation |
| HOM_REF × HOM_INV | 100% HET | inheritance validation only |
| HET × HET | 25 / 50 / 25 | **distortion / pseudo-overdominance** |
| HET × HOM_INV | 50% HET / 50% HOM_INV | segregation |
| HOM_INV × HOM_INV | 100% HOM_INV | — |

## Tests

- **2-class** (50/50) → exact two-sided binomial against p = 0.5
- **3-class** (25/50/25) → chi-square goodness-of-fit, df = 2,
  analytic p-value `P(X ≥ x) = exp(-x/2)` (stdlib, no scipy)
- **Fixed** (single class) → all-match check; any off-class offspring
  is a `fixed_violation` (pedigree-error candidate)

## Interpretation labels

| Category | Meaning |
|---|---|
| `fixed_inheritance_validation` | HOM_REF × HOM_INV cross → all offspring HET (the cleanest validation case) |
| `fixed_consistent` | any other fixed cross, all offspring match expected |
| `fixed_violation` | fixed cross with one or more off-class offspring (flag for pedigree error) |
| `segregation_consistent` | informative cross, p ≥ alpha |
| `segregation_distorted` | informative cross, p < alpha, no pseudo-overdom pattern |
| **`homozygote_depletion`** | HET × HET cross, both homozygotes under-represented (pseudo-overdominance signature) |
| `small_n` | n_offspring below `min_offspring_for_test` (default 4) |
| `ambiguous` | parental arrangement uncalled |

## Cohort-level aggregation

`summarise_lrr(records)` aggregates across all HET × HET families for
one LRR, building a single 3-class chi-square test on the pooled
observed counts. The pooled test surfaces a cohort-wide
`cohort_homozygote_depletion` signature even when each individual
family is too small to detect distortion on its own.

## What this bloc does NOT do

- It does **not** classify the inversion mechanism. A
  `homozygote_depletion` signature is consistent with
  pseudo-overdominance but also with technical artifacts
  (genotyping-quality dropout in homozygotes, viability filtering
  during hatchery rearing, etc.). Mechanism inference needs additional
  evidence layers (deleterious-burden, fitness, etc.).
- It does **not** correct for multiple testing across LRRs. The caller
  is responsible for applying FDR / Bonferroni across the candidate
  LRR set.
