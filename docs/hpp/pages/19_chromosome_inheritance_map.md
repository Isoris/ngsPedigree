# Bloc 19 — Chromosome inheritance map (independent-assortment level)

| | |
|---|---|
| **module** | `src/hpp/chromosome_inheritance.py` |
| **scope** | whole-chromosome compatibility scoring per (offspring × candidate parent × chromosome) |
| **sibling** | `del_inheritance.py` is the **LRR inheritance map** — segment-level transmission *within* a recombination-suppressed haplotype block. Two valid maps, two questions. |

## Two inheritance maps, one bigger picture

| | Chromosome inheritance (bloc 19) | LRR inheritance (bloc 12 `del_inheritance.py`) |
|---|---|---|
| **Level** | whole chromosome | segment within an LRR |
| **Unit** | independent assortment (Mendel's second law) | transmission within a contiguous haplotype block |
| **Pearson at** | individual / one-pair level | one-pair level inside the LRR |
| **Marker resolution** | every DEL on the chromosome | every DEL inside the LRR |
| **Score** | per-chromosome compatibility + Pearson r + parent-HET transmission rate | per-segment transmitted-allele trace + recombination event flags |
| **Question** | "did this parent contribute to this whole chromosome of the offspring?" | "within this LRR, which arrangement haplotype did the offspring inherit at each segment?" |

The chromosome map is the **outer** level. The LRR map is the
**inner** level (zoomed in inside a recombination-suppressed region
where the trace is cleaner because no recombination is breaking it up).

## What the chromosome map scores

For each (offspring, candidate parent, chromosome), three numbers + a
bucket:

| Quantity | Meaning |
|---|---|
| `compatibility_rate` | fraction of DEL markers with no opposite-homozygote contradiction (HOM_REF parent vs HOM_DEL offspring or vice versa). Should be ≈ 1.0 for a real parent. |
| `pearson_r` | Pearson correlation between parent and offspring DEL **dosage** on this chromosome. Positive r = inheritance signal; ≈ 0 for a stranger. |
| `het_inheritance_rate` | of all parent-HET DELs on this chromosome, fraction the offspring carries at least one copy of. **Should be ≈ 0.5** for a true parent (HET parent transmits the DEL 50% of the time). |
| `parent_hom_dels_inherited / parent_hom_dels_total` | every parent-HOM-DEL must be transmitted; offspring should carry the DEL at every such marker. Used as a Mendelian sanity check. |

Bucket (`inheritance_support`):

| Label | Rule |
|---|---|
| `rejected` | ≥ 1 opposite-homozygote contradiction (Mendelian-impossible) |
| `ambiguous` | < `min_markers` DELs called for both samples on this chromosome |
| `compatible` | Mendelian-clean, but only one of (Pearson ≥ 0.3, HET rate in [0.35, 0.65]) holds |
| `strong` | Mendelian-clean **and both** signals hold |

## Why this works (the user's intuition, formalised)

A diploid offspring receives one homolog of each chromosome from each
parent under independent assortment. So:

- Where the candidate parent is **HOM-REF (0/0)**: the homolog they
  contributed carries 0. Offspring must carry **≥ 1 zero allele** at
  that marker — so offspring is 0/0 or 0/1, never 1/1. A 1/1
  offspring at a parent-0/0 marker is a **strict Mendelian
  contradiction** → exclude.
- Where the candidate parent is **HOM-DEL (1/1)**: symmetric.
  Offspring must carry ≥ 1 DEL allele.
- Where the candidate parent is **HET (0/1)**: not informative for
  hard exclusion at a single marker, but the parent's HET DELs are a
  **per-chromosome fingerprint**. Under independent assortment of
  parental alleles, ~50% of parent's HET DELs are transmitted to the
  offspring's contributed homolog. The other homolog (from the other
  parent) is independent. So at the chromosome level, on average ~50%
  of parent's HET DELs should appear in the offspring's genotype as
  0/1 or 1/1.

The Pearson r captures both signals at once: positive dosage
correlation between parent and offspring on this chromosome is the
inheritance signature.

## What "we cannot phase it all" means here

The user's exact phrasing was: "the individual got chr 1 chr 2 chr 3
for sure from catfish A and from catfish B probably too but we need
to have score because we cannot phase it all."

The score *is* the un-phaseable answer. Without phasing, we cannot
say "this exact homolog came from A" — but we **can** say with high
confidence: "A is compatible as a contributor to this whole chromosome
(no opposite-hom contradictions, strong dosage correlation, ~50% HET
transmission), and B is also compatible (same)." That is enough to
confirm both parents of a triad on a per-chromosome basis. And when
one of them is NOT compatible — opposite-hom contradictions on chr X
but not on the others — we have evidence of either a half-sibship,
a sample swap, or a chromosomal anomaly.

## Use cases

- **Triad sanity-check.** Each of the two confirmed parents should
  score `strong` or `compatible` on every chromosome. A `rejected`
  per-chromosome score flags a problem.
- **Distinguishing PO from FS** when θ + IBS0 are borderline. A true
  PO is compatible on every chromosome; a full-sib has ~25% expected
  opposite-hom rate across chromosomes from the unshared parental
  alleles.
- **Cryptic chromosomal events.** Uniparental disomy or a sibship
  swap leaves one chromosome with the wrong parent profile while
  others look normal — visible per-chromosome but invisible to
  genome-wide θ.

## Pearson scale: population vs individual (the user's distinction)

| Bloc | Pearson computed across | Why |
|---|---|---|
| **16 (de novo LRR discovery)** | **DEL pairs across the cohort** (one r per DEL pair) | finds DELs that co-segregate at the cohort level → arrangement haplotype block |
| **19 (chromosome inheritance)** | **one offspring vs one candidate parent** on this chromosome | tests whether this specific pair has inheritance compatibility |

Both Pearson r — different denominators. Bloc 16's signal is
population-level (does this group of DELs co-segregate across many
samples?). Bloc 19's signal is individual-level (does this one
parent's dosage match this one offspring's dosage on this one
chromosome?).

## Output

`chromosome_inheritance_scores.tsv`:

| col | meaning |
|---|---|
| `offspring_sample_id`, `candidate_parent_sample_id`, `chrom` | identity |
| `n_markers_called` | DELs where both samples have a call |
| `n_compatible`, `n_excluding` | Mendelian-compatible vs opposite-hom |
| `compatibility_rate` | n_compatible / n_called |
| `pearson_r` | dosage correlation |
| `n_parent_het_markers`, `n_parent_het_inherited`, `het_inheritance_rate` | HET-DEL transmission stats |
| `parent_hom_dels_inherited`, `parent_hom_dels_total` | HOM-DEL passage check |
| `inheritance_support` | `rejected` / `ambiguous` / `compatible` / `strong` |

`chromosome_best_parent.tsv`: per (offspring, chrom), best non-rejected
candidate parent and their Pearson r.

## What this bloc does NOT do

- It does **not** phase the chromosomes. We do not claim which
  specific homolog came from which parent — only that each parent is
  compatible as a contributor.
- It does **not** call chromosome-level recombination. Switches in
  inheritance source mid-chromosome are LRR-bloc territory.
- It does **not** prove uniparental disomy or trisomy on its own.
  A `rejected` per-chromosome score is a flag, not a diagnosis.

## How it fits the manual

Section 4.7 of `docs/MANUAL.md` covers the LRR inheritance map.
The chromosome inheritance map slots in alongside as a **parallel**
analysis: run it after the candidate-PO + triad-detection step,
before the LRR-level analyses. On the master pipeline output, it
becomes the per-chromosome score sheet for each offspring against
every candidate first-degree partner.
