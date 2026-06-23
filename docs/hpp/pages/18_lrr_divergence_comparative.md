# Bloc 18 — Cross-LRR comparative test (situation 1 → genome-wide)

| | |
|---|---|
| **module** | `src/hpp/lrr_divergence_comparative.py` |
| **purpose** | per-LRR metrics + cross-LRR regression of homokaryotype depletion / FIS / heterokaryotype enrichment against between-arrangement divergence |
| **hypothesis** | more divergent LRRs show stronger homokaryotype depletion and heterokaryotype enrichment — consistent with progressive pseudo-overdominance |

## The user-stated hypothesis

> Across candidate LRRs, do older / more divergent arrangements show
> stronger homokaryotype depletion and heterokaryotype enrichment?

The bloc-17 "situation 1" pattern (a hemizygous-only DEL whose host
homokaryotype is depleted) was the anecdotal version. This bloc lifts
it to a cohort-wide regression: does the magnitude of homokaryotype
depletion scale with between-arrangement divergence?

## Per-LRR metrics

For each LRR, compute:

| Quantity | Definition |
|---|---|
| `n_hom0`, `n_het`, `n_hom1` | sample counts in each arrangement class |
| `freq_*` | class frequencies (out of called total) |
| `p_allele_1` | arrangement-1 allele frequency `(2·n_hom1 + n_het) / (2·n_total)` |
| `H_obs` | observed heterozygosity `n_het / n_total` |
| `H_exp` | expected heterozygosity `2 p (1 - p)` |
| **`FIS`** | `1 - H_obs / H_exp` — negative = het excess |
| **`heterokaryotype_enrichment`** | `H_obs / H_exp` — > 1 = het excess |
| `hom0_het_ratio`, `hom1_het_ratio` | `(n_hom* + 0.5) / (n_het + 0.5)` — Haldane-corrected |
| **`min_hom_het_ratio`** | `min(hom0_het_ratio, hom1_het_ratio)` |
| `missing_hom_class` | `"HOM_REF"` / `"HOM_INV"` / `"BOTH"` / `None` |
| **`dxy_between_arrangements`** | mean over markers inside the LRR of `p_arr0(1-p_arr1) + p_arr1(1-p_arr0)`, computed from DEL genotypes in arr-0 vs arr-1 homozygous samples |

`dxy_between_arrangements` is treated as a **relative divergence proxy**
— it is influenced by ancestral diversity, mutation rate, Ne, and
selection. It is NOT an absolute age estimate. The thesis wording
(in the regression output) makes this explicit.

## Cross-LRR regression

Three response variables tested, each against `dxy_between_arrangements`:

| y | Expected sign under pseudo-overdominance |
|---|---|
| `log(min_hom_het_ratio)` | **negative** (more divergent → stronger depletion) |
| `FIS` | **negative** (more divergent → more het excess) |
| `heterokaryotype_enrichment` | **positive** |

Fit by stdlib OLS:

```
slope     = Σ(x - x̄)(y - ȳ) / Σ(x - x̄)²
intercept = ȳ - slope · x̄
r         = Σ(x - x̄)(y - ȳ) / √(Σ(x - x̄)² · Σ(y - ȳ)²)
t_stat    = r · √(n - 2) / √(1 - r²)
p_value   = erfc(|t| / √2)   # asymptotic normal approximation
```

The p-value is normal-approximation; for `n_LRRs < 30` an asymptotic
caveat is attached so consumers know to fall back to exact Student t
externally.

## Output

`lrr_comparative_summary.tsv` — one row per LRR with every metric above.

`lrr_comparative_regression.tsv` — one row per regression model:

```
y_name                       x_name                       n_lrrs slope     r        p_asym
log_min_hom_het_ratio        dxy_between_arrangements     47     -3.241    -0.84    1.2e-12
fis                          dxy_between_arrangements     47     -2.105    -0.79    3.8e-10
heterokaryotype_enrichment   dxy_between_arrangements     47      2.480     0.81    5.4e-11
```

## Thesis-safe wording (built into the module docstring)

> We tested whether candidate LRR divergence predicted genotype-class
> imbalance by regressing homokaryotype-to-heterokaryotype ratios and
> FIS against between-arrangement dXY. This analysis evaluates
> whether heterokaryotype enrichment strengthens as recombination-
> suppressed arrangements diverge, as expected under progressive
> pseudo-overdominance. Because dXY is influenced by both age and
> ancestral diversity, we interpret it as a relative divergence proxy
> rather than a direct estimate of age.

## If the trend holds (worked phrasing for the chapter)

> More divergent LRRs showed stronger homokaryotype depletion and
> heterokaryotype enrichment, consistent with progressive accumulation
> of complementary deleterious variation under recombination
> suppression.

## If the trend does not hold

> Homokaryotype depletion was observed in some candidate LRRs but did
> not scale clearly with between-arrangement divergence, suggesting
> that pseudo-overdominance intensity may depend on local gene content,
> arrangement-specific deleterious variants, or recent demographic
> structure rather than divergence alone.

Either outcome is publishable.

## What this bloc does NOT do

- It does **not** claim that age **causes** depletion — only that the
  pattern is consistent with the evolutionary model.
- It does **not** correct for multiple testing across response
  variables; with three tested y's, the caller should apply
  Bonferroni / FDR if treating them as independent.
- It does **not** add length / sample-size / marker-density covariates
  — those are obvious next-step additions (the user proposed
  `lm(log(min_hom_het_ratio) ~ dXY + LRR_length + n_samples + SNP_density)`);
  this bloc fits the simple univariate model, leaves covariate-extended
  models to downstream R / statsmodels callers.
