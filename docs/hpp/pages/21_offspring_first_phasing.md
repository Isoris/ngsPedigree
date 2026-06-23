# Bloc 21 — Offspring-first arrangement phasing

| | |
|---|---|
| **module** | `src/hpp/offspring_phasing.py` |
| **tests** | `tests/test_offspring_phasing.py` |
| **scope** | reverse-segregation marker discovery from variant dosage in sibship cohorts. **Not SNP phasing.** Offspring are split into segregation classes using *any* segregating variant (DEL, DUP, INV, BND, indel, depth, soft probabilities) and every variant is then scored by how well it separates the classes. |
| **outputs** | `offspring_phasing.tsv` (per-variant ranked markers), `offspring_class_assignments.tsv` (per-sample class label) |
| **input** | a children × variants dosage matrix per (family × chromosome × interval) |

## The premise

> You are not phasing SNPs; you are phasing inherited parental
> arrangement *blocks* using any variant that segregates among
> offspring.

For SVs the marker is not 0/1 at one SNP. It can be:

- DEL present / absent
- DUP present / absent
- INV-support present / absent
- BND-support present / absent
- indel present / absent
- normalised coverage drop / normal
- discordant-read or split-read signal
- local haplotype cluster

So the input is a **variant-by-child matrix**, not a SNP matrix. The
dosage encoding is uniform across types:

| value | meaning |
|---|---|
| `0.00` | absent / `0/0` |
| `0.50` | heterozygous / `0/1` |
| `1.00` | homozygous present / `1/1` |
| `None` | missing |

Soft probabilities in [0, 1] (from likelihoods or normalised depth) are
accepted directly. Raw 0/1/2 integers are normalised by /2.

## Method

For each (family × chromosome × interval):

1. **Build** the children × variants dosage matrix on the interval.
2. **Filter** variants by either of two modes (see below).
3. **Cluster** the remaining offspring into two classes A / B using
   k-means-like iteration on the filtered variants. The seed pair is
   the two most distant offspring (L1 distance over jointly observed
   positions). Iteration stops when class membership is stable.
4. **Score** every kept variant by class separation:

   ```
   marker_score = | mean(dosage | class A) - mean(dosage | class B) |
   ```

5. **Rank** by `marker_score`. The top markers are the
   arrangement-tagging candidates for this interval.

The polarity of class A versus class B is arbitrary (band-relabel
symmetry, same caveat as bloc 06). Only the within-family marker score
is meaningful.

## Two filter modes, run side by side

The two modes can be invoked independently and compared.

### `segregating` (default)

Keep variants whose carrier frequency in the sibship is in
`[min_freq, max_freq]` (default 0.20–0.80). This is the classic
"variant splits children into inheritance groups" filter and catches
the common patterns:

| pattern | implied parental cross | useful as marker? |
|---|---|---|
| ~50% present, ~50% absent | `AA × AB` (one parent HET) | yes (high score) |
| ~25% absent, ~50% het, ~25% hom | `AB × AB` (both parents HET) | yes |
| 100% present | `BB × BB` | no (no offspring split) |
| singleton | de novo / error | no |

### `hemizygous_only`

Keep variants where the **majority** of informative offspring are
heterozygous (intermediate dosage). This is the situation-1
"hom-depleted" pattern: both homokaryotypes are absent or extremely
rare, so every observed offspring carries one copy.

These markers are explicitly **not** reported as evidence of lethality.
The pattern is `het × absent` *or* `hom_depleted` — biological
depletion and technical hom-DEL miscall are not distinguishable from
the sibship alone. See bloc 17 for the longer disclaimer.

Both modes feed the same downstream scoring, so the two output tables
are directly comparable. Run the segregating mode first to see where
the offspring split into classes, then run the hemizygous-only mode to
catch the arrangement-tagging markers that the segregating filter
discards because they sit at carrier-frequency 1.0.

## Honesty caveats

This module:

- does **not** phase SNPs (no LD model, no IBD inference)
- does **not** estimate arrangement age or dXY
- does **not** make any fitness, balancing-selection, or
  overdominance claim
- does **not** call lethality from hom-depleted patterns

It outputs arrangement-tagging marker *candidates* only. False
positives are caught downstream by the family-based enrichment (bloc
13) and Mendelian segregation analysis (bloc 14).

## Output schema

`offspring_phasing.tsv`:

| column | meaning |
|---|---|
| `family_id`, `chrom`, `interval_start`, `interval_end` | scope |
| `variant_id`, `variant_type`, `pos` | marker identity (DEL / DUP / INV / BND / indel / depth / SNP) |
| `n_class_a`, `n_class_b` | informative offspring per class |
| `mean_dosage_class_a`, `mean_dosage_class_b` | per-class mean dosage |
| `marker_score` | `| mean_a − mean_b |`; rank key |
| `n_informative`, `missingness` | per-variant call counts |
| `filter_mode` | `segregating` or `hemizygous_only` |
| `segregation_pattern` | `1:1` / `1:1_inverted` / `1:2:1` / `all_het` / `fixed_present` / `fixed_absent` / `intermediate` |
| `inferred_parental_state` | rough hint: `AA_x_AB` / `BB_x_AB` / `AB_x_AB` / `het_x_absent_or_hom_depleted` / `BB_x_BB` / `AA_x_AA` |
| `notes` | free text |

`offspring_class_assignments.tsv`: `family_id`, `chrom`,
`interval_start`, `interval_end`, `sample_id`, `class_label` (`A` /
`B` / `U`), `filter_mode`.

## Minimal first usage

```python
from hpp.offspring_phasing import (
    VariantMarker,
    phase_offspring_interval,
    write_phasing_tsv,
)

markers = [
    VariantMarker("DEL_20491", "LG12", 10_200_000, 10_200_500, "DEL"),
    VariantMarker("DEL_20492", "LG12",  9_500_000,  9_500_300, "DEL"),
    # ...
]

# offspring × variant dosage matrix; values can be GT strings or floats.
matrix = {
    "DEL_20491": {"FAM01_O1": "1/1", "FAM01_O2": "0/0", ...},
    # ...
}

classes_seg, recs_seg = phase_offspring_interval(
    family_id="FAM01", chrom="LG12",
    interval_start=8_100_000, interval_end=14_600_000,
    markers=markers, dosage_matrix=matrix,
    offspring=["FAM01_O1", "FAM01_O2", ...],
    mode="segregating",
)

classes_hemi, recs_hemi = phase_offspring_interval(
    family_id="FAM01", chrom="LG12",
    interval_start=8_100_000, interval_end=14_600_000,
    markers=markers, dosage_matrix=matrix,
    offspring=["FAM01_O1", "FAM01_O2", ...],
    mode="hemizygous_only",
)

write_phasing_tsv("offspring_phasing.segregating.tsv", recs_seg)
write_phasing_tsv("offspring_phasing.hemizygous_only.tsv", recs_hemi)
```

The two output tables share the same schema (with `filter_mode`
distinguishing them), so they can be `cat`-ed and rank-compared
directly.

## Where this sits in the pipeline

The offspring-first phasing is **complementary** to:

- **bloc 13** (LRR enrichment) — answers "is this LRR enriched for
  some segregation pattern across families?"
- **bloc 14** (Mendelian segregation) — answers "within this LRR,
  do families show Mendelian-consistent transmission?"
- **bloc 17** (situation-1 classifier) — answers "is this DEL
  arrangement-1-linked, possibly hom-depleted?"
- **bloc 19** (chromosome inheritance map) — answers "did this
  candidate parent contribute the whole chromosome?"

Bloc 21 answers a different question: **given offspring genotypes
alone, with no parents, can we recover arrangement-tagging markers?**
The output is fed into bloc 17's situation-1 classifier as candidate
markers when no curated LRR list is available.

## Status

MVP 1. Stdlib only, no pandas / numpy. 19 tests cover encoding, both
filter modes, two-class clustering, variant scoring, the interval
orchestrator, and TSV round-trips.
