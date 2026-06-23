# Blocs 16–17 — One-command end-to-end + de novo LRR discovery + arrangement-linkage classifier

| | |
|---|---|
| **CLI** | `scripts/10_run_full_pipeline.py` |
| **modules** | `lrr_discovery.py`, `del_arrangement_linkage.py` |
| **inputs** | one Delly VCF + one Manta VCF (uncompressed or `.gz`); optional `--list-of-lrr` OR `--discover-lrrs`; optional `--karyotype-catalogue` |
| **outputs** | merged catalogue, pairwise table, candidate PO + triads, chromosome inheritance map, discovered LRRs, enrichment table, arrangement-linkage table |

## Run the whole stack with one command

```bash
python ngsPedigree_v0.5.0/scripts/10_run_full_pipeline.py \
    --delly delly.vcf.gz \
    --manta manta.vcf.gz \
    --discover-lrrs \
    [--karyotype-catalogue cat.json] \
    --outdir results/
```

That single command runs blocs 12 → 16 (discovery) → 13 → 17 in order
and emits everything to `--outdir`. Skip `--discover-lrrs` and pass
`--list-of-lrr lrr.tsv` instead if you already have a curated LRR list.

## Bloc 16 — De novo LRR discovery from DEL correlation

When no curated LRR list is supplied, scan each chromosome in sliding
windows and emit candidate LRRs where DEL genotypes are
**strongly haplotype-linked**:

1. For each window of length `--window-size` (default 1 Mb), collect
   all DEL markers whose midpoint falls inside.
2. If the window has ≥ `--min-markers-per-window` markers, compute the
   mean **absolute pairwise Pearson correlation** of DEL genotype
   dosages (0/1/2) across samples for every pair of markers in the
   window.
3. Windows with mean correlation ≥ `--correlation-threshold`
   (default 0.5) are candidate LRRs.
4. Adjacent (or near-adjacent, within `max_merge_gap`) candidate
   windows are merged into one LRR interval.

The detector is intentionally conservative — output is a candidate LRR
list with the **same TSV shape that `--list-of-lrr` accepts**, plus
QC columns (`n_markers`, `mean_pairwise_correlation`). False positives
are caught downstream by the family-based enrichment (bloc 13).

## Bloc 17 — Per-(DEL, LRR) arrangement-linkage classifier ("situation 1")

When a `--karyotype-catalogue` is supplied (per-(sample, LRR)
arrangement classes from PC1 bands, schema
`ngspedigree_karyotype_catalogue_v1`), classify every DEL inside every
LRR. For each DEL × LRR pair compute DEL allele frequency stratified
by host arrangement:

|                       | DEL freq |
|---|---|
| HOM_REF (band 0/0) samples | f_HR |
| HET (band 0/1) samples     | f_H  |
| HOM_INV (band 1/1) samples | f_HI |

A real arrangement-linked DEL produces the expected gradient
`0.0 → 0.5 → 1.0` (arrangement-1 linked) or its mirror
(arrangement-0 linked). Interpretation labels:

| Label | Pattern |
|---|---|
| `arrangement_0_marker` | f_HR ≥ 0.85, f_HI ≤ 0.10 |
| `arrangement_1_marker` | f_HR ≤ 0.10, f_HI ≥ 0.85 |
| **`arrangement_1_marker_hom_depleted`** | f_HR ≤ 0.10, f_H ≈ 0.5, **no/few HOM_INV samples** — the "situation 1" pattern: arrangement-linked DEL only ever seen hemizygous because the host homokaryotype is absent or depleted |
| **`arrangement_0_marker_hom_depleted`** | symmetric for arrangement-0 |
| `unlinked` | all three classes intermediate |
| `ambiguous` | < 2 classes with at least 2 samples |

## Honesty pass for "situation 1"

When the homokaryotype is missing, the `notes` field carries the
exact disclaimer the user asked for:

> no/few HOM_INV samples observed; arrangement-linked DEL only seen
> hemizygous; cannot distinguish biological depletion from low
> frequency or technical hom_DEL miscall — interpret as depletion
> candidate, **not lethality**

No claim that the homozygous arrangement is lethal. The pattern is
flagged as **depletion candidate**; mechanism (low frequency vs
sampling vs viability filtering vs technical miscall) is left to
external analysis.

## Output table — bloc 17

`del_arrangement_linkage.tsv`:

| Column | Meaning |
|---|---|
| `del_id`, `lrr_id`, `chrom`, `del_start`, `del_end` | identity |
| `n_hom_ref` / `n_het` / `n_hom_inv` | sample counts per class |
| `del_freq_hom_ref` / `_het` / `_hom_inv` | DEL allele frequency per class |
| `interpretation` | one of the six labels above |
| `notes` | disclaimer text (especially for depletion cases) |

## Bloc 12 + 13 + 14 + 16 + 17 — what the master CLI actually runs

| Step | Bloc | Triggered when |
|---|---|---|
| read & merge Delly + Manta VCFs | 12 | always |
| KING-robust pairwise relatedness | 12 | always |
| Mendelian-exclusion + forced-offspring triad detection | 12 | always |
| chromosome inheritance map (per triad) | 12 | always |
| de novo LRR discovery | **16** | `--discover-lrrs` |
| family-based OR enrichment | 13 | LRR set is non-empty |
| per-DEL × per-LRR arrangement linkage | **17** | `--karyotype-catalogue` supplied |

300 unit tests cover the path.

## What this still does NOT do

- It does **not** call recombination tracts inside the LRR — that is
  ngsTracts' job (downstream).
- It does **not** prove molecular inversion breakpoints — that needs
  independent SV-caller / assembly evidence.
- It does **not** prove that hemizygous-only DELs reflect lethal
  homokaryotypes — the explicit disclaimer in bloc 17's `notes` field
  rules out that overclaim.
- It does **not** replace local-PCA discovery — bloc 16 is the broke-
  grad-student stand-in for population-genomic LRR discovery and is
  expected to be less sensitive than proper local-PCA / Cramér's V.
