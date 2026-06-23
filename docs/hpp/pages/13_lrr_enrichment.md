# Bloc 13 — Family-based OR enrichment for candidate LRRs

| | |
|---|---|
| **module** | `src/hpp/lrr_enrichment.py` |
| **CLI** | `scripts/09_lrr_enrichment.py` |
| **inputs** | merged DEL catalogue (bloc 12), LRR list TSV (`--list_of_LRR`), confirmed triads / dyads |
| **purpose** | quantify whether each candidate LRR is enriched for transmission-compatible block inheritance, vs matched background |

## The correct claim

> Family-based enrichment supports an inversion-like recombination-suppressed
> haplotype block — it does **not** alone prove the molecular inversion
> breakpoint.

## 2 × 2 table per LRR

| | block-compatible | not |
|---|---|---|
| inside candidate LRR | **a** | **b** |
| matched background | **c** | **d** |

`block-compatible` = within the region, ≥ `min_markers` informative
parent-HET markers, dominant transmitted allele covers ≥
`dominance_threshold` (default 0.8) of those markers, **and** zero
Mendelian contradictions.

OR = (a/b)/(c/d). 95% CI via Woolf log-OR SE with Haldane–Anscombe
correction when any cell is zero.

Three ORs reported:
- **combined_OR** — all relationships pooled
- **triad_OR** — triads only (stronger evidence per region)
- **dyad_OR** — dyads only (weaker; informative direction-aware only)

## Input contract

LRR list TSV (`--list-of-lrr` / `--list_of_LRR`):

```
lrr_id	chrom	start	end
LRR_001	Chr1	100000	500000
LRR_002	Chr1	1200000	1800000
```

Triads / dyads as TSV (or directly from bloc 12's `candidate_PO_pairs.tsv`).

## Background sampling

For each LRR, sample `--n-background` non-overlapping windows of the
same length on the same chromosome, avoiding any LRR interval.
Deterministic via `--seed`.

## Interpretation

- **OR ≈ 1** — LRR is not special; behaves like random regions
- **OR > 1** — enriched for inherited block behavior
- **OR >> 1 with tight CI** — strong family-based support
- **OR < 1** — candidate is noisy or wrong

## What this does NOT prove

The molecular inversion breakpoint. For that, you need split-read /
assembly / inversion-specific SV evidence — separate analysis,
separate bloc.
