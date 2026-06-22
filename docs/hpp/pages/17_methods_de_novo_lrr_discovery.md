# Methods — how de novo LRR discovery works from DEL VCFs alone

This page explains in detail how bloc 16 detects candidate LRR
(low-recombination region / inversion-like haplotype block) intervals
**without** any user-supplied dosage table, karyotype call, or
local-PCA output. The only inputs are Delly + Manta DEL VCFs.

The intuition behind why this works is small but specific. Let me
unpack it step by step.

## Step 0 — the only thing we have

Two structural-variant VCFs:

```
delly.vcf:
#CHROM  POS    ID            REF ALT    QUAL FILTER INFO                       FORMAT P_F P_M C  UNREL
Chr1    1000   DEL00000001   N   <DEL>  100  PASS   SVTYPE=DEL;END=1500;...   GT     0/1 0/0 0/0 0/0
Chr1    5000   DEL00000002   N   <DEL>  100  PASS   SVTYPE=DEL;END=5500;...   GT     0/0 0/1 0/0 0/0
...

manta.vcf:
... similar structure ...
```

Each row carries a DEL call with a per-sample genotype: `0/0` (no
deletion), `0/1` (heterozygous / hemizygous), `1/1` (homozygous DEL).

**That's it.** No PCA bands. No karyotype labels. No HOM/HET annotation
provided by you. Just per-sample 0/0, 0/1, 1/1 genotypes for each DEL
call.

## Step 1 — collapse the genotype to a numeric dosage

For Mendelian work, the genotype is conceptually a count of DEL
alleles:

| Genotype | Dosage |
|---|---|
| 0/0 (no DEL) | 0 |
| 0/1 (hemizygous) | 1 |
| 1/1 (homozygous DEL) | 2 |
| ./. (missing) | drop |

So each DEL marker becomes a **vector of integers**, one per sample:

```
DEL_001 dosage vector across [P_F, P_M, C, UNREL] = [1, 0, 0, 0]
DEL_002 dosage vector                              = [0, 1, 0, 0]
DEL_003 dosage vector                              = [1, 1, 2, 0]
DEL_004 dosage vector                              = [2, 0, 1, 0]
...
```

That's the per-sample dosage you said you "didn't give me". It's
already in the VCF — the genotype encoding **is** the dosage. We
don't need you to tell us anything beyond the VCF.

## Step 2 — the core signal: DELs on the same haplotype co-segregate

Here is the central observation.

An LRR is, by definition, a region where recombination is suppressed
— typically because there are two structurally distinct arrangement
haplotypes (REF vs INV) that don't pair up cleanly at meiosis. In such
a region, **every variant that arose on one arrangement haplotype gets
carried along with that arrangement haplotype** in every meiosis. They
don't get separated by crossing-over.

Including DELs.

So if DEL A and DEL B both happen to lie on the INV arrangement of an
LRR, every sample that carries that INV arrangement once carries both
DELs once; every sample that carries it twice carries both DELs twice;
every sample that carries it zero times carries neither.

Concretely, for two DELs both on arrangement-1 of the same LRR:

| Sample arrangement | DEL A dosage | DEL B dosage |
|---|---|---|
| arrangement 0/0 | 0 | 0 |
| arrangement 0/1 | 1 | 1 |
| arrangement 1/1 | 2 | 2 |

So **DEL A's dosage vector and DEL B's dosage vector are identical
across samples** (up to genotyping noise). Their Pearson correlation is
≈ 1.

Outside an LRR, recombination cuts variants apart. Two DELs that happen
to be near each other but on different haplotypes — or that just have
no special linkage — will show **uncorrelated dosage vectors across the
cohort**. Their Pearson correlation hovers around 0 (with some scatter
depending on sample size).

So the signal is:

> **DELs whose dosage vectors are highly correlated across samples
> are sitting on the same haplotype block. A genomic window with many
> DELs whose pairwise dosage correlations are all high is a candidate
> LRR.**

You don't need to know which sample is HOM or HET. You just need to
notice that **DEL A's pattern matches DEL B's pattern matches DEL C's
pattern**, which is the LRR fingerprint.

## Step 3 — quantify "correlated" with Pearson r

For two dosage vectors `x = (x_1, x_2, …, x_n)` and `y = (y_1, …, y_n)`
across `n` samples:

```
            sum_i (x_i - mean(x)) * (y_i - mean(y))
r(x, y) = ─────────────────────────────────────────────
          sqrt( sum_i (x_i - mean(x))^2 * sum_j (y_j - mean(y))^2 )
```

- `r = +1` : perfect positive co-segregation (DEL A on the same
  arrangement as DEL B)
- `r = -1` : perfect anti-correlation (DEL A on arrangement 1, DEL B on
  arrangement 0 — also strong LRR signal, just flipped)
- `r ≈ 0` : independent / unlinked

In bloc 16 we use **|r|** (absolute value) so that anti-correlated
markers count as LRR-supporting too. Both arrangements' markers
contribute to the signal.

The implementing function is
`hpp.lrr_discovery.pearson_corr(xs, ys)` and `marker_pair_correlation`
above it.

## Step 4 — average the pairwise correlations inside a window

A single highly-correlated pair could be a fluke (two DELs nearby that
just happen to look alike). The LRR signature is that **many DELs in
a region all pairwise-correlate**. So we average:

For a window with markers `M_1, …, M_k`:

```
window_score = (1 / (k choose 2)) * sum over all pairs (i, j) of |r(M_i, M_j)|
```

That number lives in [0, 1]. The closer to 1, the more uniformly the
DELs in the window are linked.

For the synthetic fixture:
- a Chr1 window with 16 DELs whose dosage patterns mostly co-segregate
  scores `window_score ≈ 0.64`
- a Chr2 window with 7 DELs whose dosage patterns partly co-segregate
  scores `≈ 0.50`
- a window of 6 randomly-genotyped DELs scores `~ 0.2`

## Step 5 — sliding-window scan + adjacency merging

Walk each chromosome in steps of `--window-size / 2` (default step =
500 kb when window = 1 Mb). For each window:

1. Collect all DEL markers whose midpoint falls inside the window.
2. Skip windows with fewer than `--min-markers-per-window` markers
   (default 4) — too few pairs to compute a stable mean.
3. Compute the window's mean |r|.
4. If it is ≥ `--correlation-threshold` (default 0.5), emit the
   window as a candidate LRR.

Then merge: any two emitted windows that overlap or sit within
`max_merge_gap` (default 500 kb) are combined into one LRR interval.
The merged span runs from the earliest start to the latest end of any
constituent window.

The output is a TSV in **exactly the same shape that `--list-of-lrr`
consumes**, plus QC columns:

```
lrr_id     chrom  start    end       n_markers  mean_pairwise_correlation  notes
cLRR_0001  Chr1   1269     101269    16         0.6428                      merged
cLRR_0002  Chr2   1249     51249     7          0.5045
```

These flow straight into bloc 13's family-based OR enrichment — same
shape, same downstream.

## Step 6 — why this works without you telling us "this is HOM, this is HET"

The information you thought you needed to provide is **already implicit
in the DEL genotype patterns themselves**:

- "Which samples are HOM_REF for this arrangement?" → the ones whose
  dosage is 0 at all DELs sitting on the INV arrangement.
- "Which samples are HET?" → those whose dosage is 1 at all DELs
  sitting on the INV arrangement.
- "Which samples are HOM_INV?" → those whose dosage is 2 across the
  DEL pattern.

We don't need you to label them, because the pattern itself recovers
the labelling — the three groups fall out as natural clusters of
the dosage vectors. (We don't actually cluster the samples in bloc 16
— we just measure that the DELs co-segregate, which is the LRR
signature. Sample-level arrangement calling is what bloc 17 does, IF
you give us a karyotype catalogue.)

## Step 7 — what happens with no karyotype catalogue at all

If you do not provide `--karyotype-catalogue`, the pipeline still
works through bloc 16 (discovery) and bloc 13 (family-based enrichment).
What it cannot do without a catalogue is bloc 17 (per-DEL
arrangement-linkage labelling — the "situation 1" table), because that
specifically requires knowing whether each sample is HOM_REF / HET /
HOM_INV at each LRR.

The pipeline degrades gracefully:

| With | We can also produce |
|---|---|
| just VCFs | merged catalogue, pairwise relatedness, dyads/triads, inheritance map |
| + `--discover-lrrs` | candidate LRR list, OR enrichment |
| + `--list-of-lrr` curated | same, but using the curated set |
| + `--karyotype-catalogue` | per-(DEL, LRR) arrangement-linkage table (situation 1) |

## A worked example — synthetic_svvcfs fixture

8 PASS DEL records on Chr1 + 7 on Chr2 from Delly, 5 on Chr1 + 4 on
Chr2 from Manta. Merge by 500 bp tolerance + 0.5 reciprocal overlap →
18 unified markers (6 found by both callers, 12 by one).

Cohort: 4 samples (P_F, P_M, C, UNREL). For each DEL marker, build the
dosage vector across those four samples:

```
DEL @ Chr1:1000-1500    dosage [1, 0, 0, 0]
DEL @ Chr1:5000-5500    dosage [0, 1, 0, 0]
DEL @ Chr1:10000-10500  dosage [1, 1, 2, 0]
DEL @ Chr1:20000-21000  dosage [2, 0, 1, 0]
DEL @ Chr1:30000-30500  dosage [0, 2, 1, 2]
DEL @ Chr1:40000-40500  dosage [1, 0, 0, 0]
DEL @ Chr1:50000-50800  dosage [1, 0, 1, 0]
...
```

Slide a 50 kb window along Chr1. The first window (1–50 kb) contains
~8 markers. Compute all C(8, 2) = 28 pairwise |r|. Some of these pairs
co-segregate (e.g. DEL @ 1000 and DEL @ 40000 share `[1, 0, 0, 0]`,
r = 1.0). Others don't. The mean |r| works out to 0.64 — above the
0.4 threshold — so this window is emitted as a candidate LRR.

Run again on Chr2 → second candidate emerges.

Adjacent candidate windows merge. Final output:
`cLRR_0001 Chr1 1269-101269` and `cLRR_0002 Chr2 1249-51249`. Each
gets fed downstream into the family-based enrichment.

## Limitations and honest caveats

- **Sample size matters.** With < ~10 samples, Pearson correlations are
  noisy. The synthetic fixture has 4 samples, which is why the discovery
  threshold has to be lowered to 0.4. On real cohort data (≥ 30
  samples) the default 0.5 is realistic.
- **DEL density matters.** A real LRR with only one or two DEL markers
  inside will not pass `min_markers_per_window`. Increase
  `--window-size` for sparse regions; decrease `--min-markers-per-window`
  if you want sensitivity at the cost of more false positives.
- **No replacement for proper population-genomic discovery.** This is
  the broke-grad-student stand-in for local-PCA / Cramér's V / long-
  range votes. Real local-PCA on SNPs is more sensitive. Bloc 16's job
  is to give you something to feed the family-based test when the
  PCA-based discovery isn't available.
- **It does not prove the inversion.** A high-correlation cluster is
  consistent with any haplotype-block process: real inversion, real
  LRR around a centromere, a misassembly that locks DEL calls
  together, etc. Family-based OR enrichment (bloc 13) is the test
  that separates real inherited blocks from these artefacts.

## TL;DR

The dosage table you thought you needed to provide — "this sample is
HOM, this sample is HET" — **already exists in the Delly / Manta VCF
genotypes**: the per-sample `0/0` / `0/1` / `1/1` calls are the
dosage. We exploit the fact that DELs on the same haplotype
co-segregate across samples, measure that co-segregation as average
pairwise correlation in a sliding window, and call high-correlation
windows as candidate LRRs. No external dosage call, no PCA input, no
karyotype label is required.

If you do provide a karyotype catalogue, we additionally produce the
per-DEL × per-LRR arrangement-linkage table (bloc 17) that surfaces
the "situation 1" hemizygous-only pattern with the not-lethality
disclaimer.
