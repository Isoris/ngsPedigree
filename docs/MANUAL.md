# ngsPedigree — full pipeline manual (start to finish)

This document is the single read-this-first guide. It covers the whole
chain end-to-end: what the pipeline does, why it works, what to feed
it, how to run it with one command, what each output means, and how
to interpret results for the thesis. Every step that has its own bloc
page is referenced by number; cross-reference if you want the deep
math.

---

## 0. What this pipeline actually does

You start with two Structural-Variant VCFs that you have already
produced (Delly2 + Manta). You finish with:

1. **A reconstructed pedigree** — first-degree parent-offspring pairs,
   detected triads, mtDNA-validated maternal lines when supplied.
2. **A chromosome inheritance map** per confirmed parent → offspring
   relationship — which arrangement was transmitted, segment by
   segment, with candidate recombination break points.
3. **A list of candidate LRRs** (low-recombination regions /
   inversion-like haplotype blocks), either supplied by you or
   discovered de novo from the DEL data.
4. **Family-based evidence per LRR** — odds-ratio enrichment for
   transmission-compatible block inheritance, Mendelian segregation
   per family-cross-type, the situation-1 hemizygous-only depletion
   pattern, and a genome-wide regression of depletion against
   between-arrangement divergence.
5. **A drop-in handoff to ngsTracts** — the OUT JSON that gives
   ngsTracts the phase polarity it needs to call CO/NCO/DCO inside
   the inversion intervals.

**No HPC. No ngsRelate. No BEAGLE. stdlib Python only.**
Real PhD-defensible analysis, free.

---

## 1. The biology, in one paragraph

A real inversion (or other recombination-suppressed LRR) creates two
distinct chromosomal arrangements that don't recombine cleanly at
meiosis. Whatever variants happen to lie on one arrangement haplotype
**stay there across generations** — including deletions detected by
Delly + Manta. So:

- **DELs on the same arrangement co-segregate across samples** →
  this is the de novo discovery signal (bloc 16).
- **Parent → offspring transmission of these DELs builds an inheritance
  map** at DEL-marker resolution (bloc 12).
- **DEL Mendelian exclusion** (a parent who is 0/0 cannot give a child
  1/1) lets us identify which candidate is the real parent (bloc 12).
- **Per (family, LRR) Mendelian segregation** lets us test arrangement
  inheritance the same way classical genetics tests SNP segregation
  (bloc 14).
- **Family-based odds ratios** quantify whether the candidate LRR
  behaves as an inherited block compared with matched random regions
  (bloc 13).
- **Situation 1** — when a DEL is only ever seen hemizygous because
  the homokaryotype is depleted — points to pseudo-overdominance,
  testable across LRRs by regressing depletion against divergence
  (blocs 17 + 18).

---

## 2. What you need to provide

Minimum to run anything:

| File | Why | Required? |
|---|---|---|
| `delly.vcf` or `.vcf.gz` | DEL marker calls + per-sample GT | **yes** |
| `manta.vcf` or `.vcf.gz` | DEL marker calls + per-sample GT | **yes** |

Optional, unlocks more analyses:

| File | Unlocks |
|---|---|
| `--list-of-lrr lrr.tsv` (curated LRR intervals) | family-based OR enrichment (bloc 13) and Mendelian segregation per LRR (bloc 14) on a known LRR set |
| `--discover-lrrs` (no file needed) | de novo LRR discovery (bloc 16) — emits its own LRR list, then runs blocs 13/14 on it |
| `--karyotype-catalogue catalogue.json` | per-DEL arrangement-linkage classifier (bloc 17 — situation 1) and the cross-LRR comparative regression (bloc 18) |
| `--mtdna mtdna.tsv` | maternal-lineage pre-flight before polarization (bloc 07) |
| `--triads triads.tsv` / `--dyads dyads.tsv` | use an externally-built pedigree instead of the SV-derived one |

That's it. No PCA output, no dosage table, no SNP VCF. You can also
supply karyotype calls if you have them — but you don't have to.

---

## 3. The one-command run

```bash
python ngsPedigree_v0.5.0/scripts/10_run_full_pipeline.py \
    --delly delly.vcf.gz \
    --manta manta.vcf.gz \
    --discover-lrrs \
    [--karyotype-catalogue cat.json] \
    --outdir results/
```

That's the whole thing. Replace `--discover-lrrs` with
`--list-of-lrr lrr.tsv` if you have a curated LRR list (e.g. from
local-PCA work).

When it finishes, `results/` contains:

```
merged_del_catalogue.json
pairwise_relationship_classification.tsv      ← Stage-1-shape table
candidate_PO_pairs.tsv
detected_triads.tsv
inheritance_segments.tsv                       ← the chromosome map
discovered_lrrs.tsv                            ← (if --discover-lrrs)
lrr_enrichment.tsv                             ← (if any LRRs)
del_arrangement_linkage.tsv                    ← (if --karyotype-catalogue)
```

Section 6 below explains every output. First let me walk through
what's actually happening inside the script.

---

## 4. Methodology, step by step

### 4.1 Read both VCFs (`hpp.vcf_sv`, bloc 12)

For each VCF, walk the records. Keep only `SVTYPE=DEL` (or `ALT=<DEL>`)
with `FILTER=PASS`. Extract per-sample `GT` field. Detect which caller
produced the VCF from the `##source` header or the filename.

What this produces internally: a list of `SvRecord(chrom, pos, end,
sv_type, qual, caller, genotypes)`, one per kept call.

### 4.2 Merge Delly + Manta into a unified catalogue (`hpp.catalogue_merge`, bloc 12)

Two records merge when **all three** are true:
- same chromosome
- left breakpoint within ±500 bp (`--bp-tolerance`)
- right breakpoint within ±500 bp
- reciprocal overlap ≥ 0.5 (`--reciprocal-overlap`)

For matched pairs, take consensus midpoint coordinates and reconcile
genotypes — agreement is kept, disagreement becomes `./.` (uncertain).
Unmatched Delly and Manta records pass through as solo entries.

What this produces: `MergedMarker(marker_id, chrom, start, end,
callers, n_callers, qual_max, genotypes)`. The `marker_id` is a stable
ID like `DEL_Chr1_100000_105000` that the rest of the pipeline uses
to identify markers.

A genotype matrix `{marker_id: {sample_id: "0/0"|"0/1"|"1/1"}}` falls
out of this.

### 4.3 Pairwise relatedness from DEL markers (`hpp.del_relatedness`, bloc 12)

For every pair of samples, compute the KING-robust kinship estimator
on the DEL genotypes:

```
theta = (N_both_het - 2 * N_opposite_hom) / (N_a_het + N_b_het)
IBS0  = N_opposite_hom / N_informative
```

where N_both_het is the count of markers where both samples are 0/1,
N_opposite_hom is the count where one is 0/0 and the other is 1/1,
and N_a_het / N_b_het are the per-sample HET counts.

This is exactly the same estimator ngsRelate / KING use on SNPs —
it works on biallelic markers regardless of whether they're SNPs or
DELs. The standard thresholds (Manichaikul 2010) still apply:

| Relationship | Expected theta | Expected IBS0 |
|---|---|---|
| duplicate / MZ twin | 0.50 | 0.000 |
| parent–offspring (PO) | 0.25 | ≈ 0 (always share ≥ 1 allele) |
| full sibling | 0.25 | ≈ 0.0125 |
| half-sib / second-degree | 0.125 | ≈ 0.04 |
| third-degree | 0.0625 | ≈ 0.08 |
| unrelated | 0 | high |

### 4.4 Classify each edge (Stage 1 logic, `hpp.relatedness_sim.classify_edge_stdlib`)

First-match-wins decision tree, exactly the same as Stage 1's classifier:

```
1. duplicate_or_clone        theta ≥ 0.45  AND IBS0 < 0.001
2. parent_offspring          theta ≥ 0.177 AND IBS0 < 0.005
3. full_sibling              theta ≥ 0.177 AND IBS0 ≥ 0.005
4. ambiguous_first_degree    theta ≥ 0.177 AND IBS0 missing
5. second_degree             theta ≥ 0.0884
6. third_degree              theta ≥ 0.0442
7. unrelated                 theta < 0.0442
```

The output is a Stage-1-shape TSV that any downstream consumer
(bloc 08, `STEP_PED_01`, etc.) can read directly.

### 4.5 Mendelian-exclusion parent identification (`hpp.del_inheritance.exclude_as_parent`, bloc 12)

For each candidate PO pair, test both directions. At any marker where:
- candidate parent is 0/0 and child is 1/1 → impossible
- candidate parent is 1/1 and child is 0/0 → impossible

A single high-confidence opposite-homozygote excludes the candidate.
With short-read data and `--min-excluding 1` (the default), this is
strict.

**Dyad direction is symmetric** when there are no opposite-hom markers
in either direction. In that case the script reports `assigned_parent
= ambiguous`. This is the documented `po_dyad_only` situation — a
known limitation of pedigree work without grandparents, sex info, or
the next step.

### 4.6 Triad detection (forced-offspring rule, `scripts/10_run_full_pipeline.py`)

A sample that has two candidate PO partners who are themselves
**unrelated** (or 2nd/3rd-degree) is the offspring; the two partners
are the parents. This is the natural complement to Stage 1's
forced-parent rule.

This recovers triads even when individual dyad directions were
ambiguous, because the triad's structure resolves it.

### 4.7 Two inheritance maps — chromosome and LRR (blocs 19 + 12)

**There are TWO inheritance maps. They answer different questions.**

| Map | Level | Module | Pearson at | Question |
|---|---|---|---|---|
| Chromosome (bloc 19) | whole chromosome | `chromosome_inheritance.py` | individual / one-pair | "did this parent contribute to this whole chromosome of the offspring?" |
| LRR (bloc 12) | segment within an LRR | `del_inheritance.py` | individual inside the LRR | "within this recombination-suppressed region, which arrangement haplotype was transmitted at each segment?" |

The chromosome map is the **outer** level (independent chromosomal
assortment under Mendel's second law). The LRR map is the **inner**
level (zoomed in inside a recombination-suppressed block where the
trace is clean because no recombination is breaking it up).

#### 4.7a Chromosome inheritance map (bloc 19)

For each (offspring × candidate parent × chromosome), compute:

- **`compatibility_rate`** — fraction of DEL markers with no
  opposite-homozygote contradiction. ≈ 1.0 for a real parent;
  < 1.0 for a non-parent.
- **`pearson_r`** — Pearson correlation of parent vs offspring DEL
  dosage on this chromosome. Positive r is the inheritance signal.
- **`het_inheritance_rate`** — fraction of parent's HET DELs the
  offspring carries at least one copy of. ≈ 0.5 for a true parent
  (HET parent transmits the DEL 50% of the time).
- **`inheritance_support`** ∈ {`rejected`, `ambiguous`, `compatible`,
  `strong`} — bucket of the above three numbers.

`rejected` means at least one opposite-homozygote contradiction at a
marker — Mendelian-impossible for a true parent. `strong` means
Mendelian-clean AND positive Pearson r AND HET transmission near 0.5.

Without phasing, the per-chromosome compatibility score is the
strongest statement we can make about "this parent contributed to
this chromosome of this offspring." For a confirmed triad, **both
parents should score `strong` or `compatible` on every chromosome**.

#### 4.7b LRR inheritance map (bloc 12, what `del_inheritance.py` does)

Walks every chromosome through all DEL markers in order. At each
marker, determines which parental allele was transmitted (REF or DEL).
Run-length encodes into segments. Switches between segments are
candidate recombination break points. Output is segment-by-segment,
with confidence Gold / Silver / Bronze based on marker density.

Within a recombination-suppressed LRR, the trace is clean (long runs
of same-allele transmissions). Outside an LRR, the trace switches
with every recombination event, which is exactly what reveals the
recombination map.

#### When to use which

- **Use the chromosome map** as a parent-of-origin sanity-check on
  every confirmed PO pair. It catches uniparental disomy, sample
  swaps, and half-sibship masquerading as PO.
- **Use the LRR map** to seed ngsTracts with phase polarity per
  (parent, offspring, LRR), and to identify candidate recombination
  break points within recombination-suppressed regions.

#### Walking the LRR map (the `del_inheritance.py` algorithm)

For each triad (or unambiguous dyad), walk every chromosome through
all the DEL markers in coordinate order. At each marker:

| Parent | Offspring | Transmitted allele |
|---|---|---|
| 0/0 | anything | REF (always) |
| 1/1 | anything | DEL (always) |
| 0/1 | 0/0 | REF (offspring received parent's REF) |
| 0/1 | 1/1 | DEL |
| 0/1 | 0/1 | ambiguous (could be either) |

In triad mode, the co-parent's homozygous state lets us disambiguate
parent-HET + offspring-HET sites that would be ambiguous in a dyad.

Run-length encode the transmitted-allele trace per chromosome.
Adjacent same-allele markers cluster into segments. Switches from REF
to DEL (or vice versa) are candidate recombination break points.
Short isolated minority runs (length < `min_run_length`) flanked by
the same majority neighbour get smoothed away as likely noise.

Output is a TSV with one row per (relationship, chromosome, segment),
carrying segment coordinates, transmitted allele, marker count,
confidence (Gold ≥8 markers, Silver 3–7, Bronze 1–2), and the
recombination-event flags.

This is the inheritance map you wanted: "this catfish inherited this
arrangement from this parent in this chromosome segment, and that
arrangement from that parent in this other segment."

### 4.8 LRR set: curated or discovered (`hpp.lrr_discovery`, bloc 16)

If you supply `--list-of-lrr`, the script reads your TSV
(columns: `lrr_id`, `chrom`, `start`, `end`).

If you supply `--discover-lrrs`, the script scans each chromosome in
sliding windows. **The signal** (see bloc 17m for the full deep-dive):
DELs sitting on the same arrangement haplotype have identical dosage
vectors across samples, so their Pearson correlation is ≈ 1. Random
DELs are uncorrelated.

For each window:
- collect DEL markers whose midpoint falls inside,
- if ≥ `--min-markers-per-window` markers, compute the mean absolute
  pairwise Pearson `r` across all C(k, 2) marker pairs,
- if the mean exceeds `--correlation-threshold` (default 0.5),
  the window is a candidate.

Adjacent or near-adjacent candidate windows merge into one LRR interval.
The output TSV has the same shape as the `--list-of-lrr` input, plus
QC columns (`n_markers`, `mean_pairwise_correlation`).

This is where your "200,000+ DELs gives enough power" intuition matters
— with a dense catalogue, even windows with 20+ DELs have C(20, 2) =
190 pairwise correlations averaging in, which gives a low-variance
window score.

### 4.9 Family-based odds-ratio enrichment (`hpp.lrr_enrichment`, bloc 13)

For every LRR, classify each (relationship, region) as
**block-compatible** when:
- ≥ `min_markers` informative parent-HET markers in the region,
- dominant transmitted allele covers ≥ `dominance_threshold` of those
  markers (default 0.8),
- zero Mendelian contradictions.

For each LRR, sample `--n-background` matched windows on the same
chromosome avoiding any LRR, classify those too, and build the 2×2
table:

|   | block-compat | not |
|---|---|---|
| inside LRR | a | b |
| matched background | c | d |

OR = (a/b)/(c/d). Woolf log-OR 95% CI; Haldane-Anscombe (+0.5)
correction when any cell is zero.

Three ORs reported per LRR: combined, triad-only, dyad-only.

**Interpretation:** OR > 1 means the LRR is enriched for inherited
block behaviour — consistent with a real inversion / LRR. **It does
not prove the molecular breakpoint** — that needs split-read /
assembly evidence (separate analysis).

### 4.10 Per-family × per-LRR Mendelian segregation (`hpp.mendelian_segregation`, bloc 14)

For every (family, LRR) pair where the LRR's arrangement-class
assignments are known per sample (from the karyotype catalogue), or
inferred from the inheritance map, classify the parental cross type:

| Cross | Expected offspring | Informative for |
|---|---|---|
| HOM_REF × HOM_REF | 100% HOM_REF | — |
| HOM_REF × HET | 50/50 | segregation |
| HOM_REF × HOM_INV | 100% HET | inheritance validation only |
| **HET × HET** | **25/50/25** | **distortion / pseudo-overdominance** |
| HET × HOM_INV | 50/50 | segregation |
| HOM_INV × HOM_INV | 100% HOM_INV | — |

Test:
- 2-class crosses → exact two-sided binomial against p = 0.5
- 3-class HET × HET → chi-square df=2, analytic
  `P(X ≥ x) = exp(-x/2)`
- fixed crosses → all-match check

Interpretation labels include `fixed_inheritance_validation`,
`segregation_consistent`, `segregation_distorted`,
**`homozygote_depletion`** (HET × HET cross with both homs
under-represented — pseudo-overdominance signature),
`fixed_violation` (pedigree-error candidate), `small_n`, `ambiguous`.

Cohort aggregation pools every HET × HET family at each LRR into one
chi-square; this surfaces cohort-wide
`cohort_homozygote_depletion` even when each sibship is too small to
detect distortion on its own.

### 4.11 Situation 1 — per-DEL arrangement-linkage (`hpp.del_arrangement_linkage`, bloc 17)

When `--karyotype-catalogue` is supplied (per-sample HOM_REF / HET /
HOM_INV labels per LRR), classify every DEL inside every LRR.

For each (DEL, LRR), compute DEL allele frequency stratified by the
host's arrangement class. A real arrangement-1-linked DEL should give:

| Arrangement class | DEL frequency |
|---|---|
| HOM_REF samples | ≈ 0.0 |
| HET samples | ≈ 0.5 |
| HOM_INV samples | ≈ 1.0 |

Interpretation labels:

- `arrangement_0_marker` / `arrangement_1_marker` — clean linkage
- **`arrangement_1_marker_hom_depleted`** — DEL tracks arrangement 1,
  HET samples are hemizygous, but **no/few HOM_INV samples observed**.
  This is the "situation 1" pattern.
- `arrangement_0_marker_hom_depleted` — symmetric
- `unlinked` — DEL frequency similar across all classes
- `ambiguous` — insufficient samples in one or more classes

The honesty note attached to every `hom_depleted` record:

> no/few HOM_INV samples observed; arrangement-linked DEL only seen
> hemizygous; cannot distinguish biological depletion from low
> frequency or technical hom_DEL miscall — interpret as depletion
> candidate, **not lethality**

### 4.12 Cross-LRR comparative test (`hpp.lrr_divergence_comparative`, bloc 18)

The situation-1 pattern, lifted to a genome-wide hypothesis:

**Does the magnitude of homokaryotype depletion scale with
between-arrangement divergence?**

For every LRR, compute:
- `n_hom0`, `n_het`, `n_hom1` (genotype-class counts)
- `FIS` = 1 − H_obs / H_exp
- `heterokaryotype_enrichment` = H_obs / H_exp
- `min_hom_het_ratio` = `min(hom0_het, hom1_het)` with Haldane (+0.5)
- `missing_hom_class` flag
- `dxy_between_arrangements` — mean across DEL markers in the LRR of
  `p_arr0·(1-p_arr1) + p_arr1·(1-p_arr0)`, where the two groups are
  the arr-0 and arr-1 homozygous samples

Then fit, across all LRRs:
- `log(min_hom_het_ratio)` ~ dXY → expect **negative slope**
- `FIS` ~ dXY → expect **negative slope** (more het excess)
- `heterokaryotype_enrichment` ~ dXY → expect **positive slope**

p-value is asymptotic normal-approximation t-test; flagged when
n_LRRs < 30.

**dXY is interpreted as a "relative divergence proxy", not an
absolute age.** Thesis-safe wording is in the module docstring.

### 4.13 ngsTracts hand-off (bloc 06 — when LRR-side analyses are run with `--karyotype-catalogue` + polarity hint)

Bloc 06 (`scripts/05_polarize_inversion.py`) produces the actual
`polarized_transmissions.out.json` that ngsTracts consumes. The
master pipeline doesn't run it automatically — it's a separate step
because ngsTracts needs an external polarity hint per LRR (which
arrangement is REF vs INV — see bloc 06 doc page on the symmetry
caveat). When you have polarity hints, the chain is:

```bash
python ngsPedigree_v0.5.0/scripts/05_polarize_inversion.py \
    --in   polarization_in.json   # assembled by 06_build_polarization_input.py
    --mtdna mtdna.tsv             # optional bloc 07 pre-flight
    --out  ngstracts_input.json   # bloc 06 OUT JSON, ngsTracts consumes
```

---

## 5. End-to-end command, fully spelled out

The realistic catfish workflow:

```bash
# Inputs you have:
#   delly.vcf.gz
#   manta.vcf.gz
#   (optional) lrr_list.tsv   from your local-PCA / curated work
#   (optional) catalogue.json  per-(sample, LRR) HOM_REF/HET/HOM_INV
#   (optional) mtdna.tsv       sample_id, mtdna_haplotype

# One-command run, de novo LRR discovery + arrangement linkage:
python ngsPedigree_v0.5.0/scripts/10_run_full_pipeline.py \
    --delly delly.vcf.gz \
    --manta manta.vcf.gz \
    --discover-lrrs \
    --karyotype-catalogue catalogue.json \
    --outdir results/

# Or with your own LRR list:
python ngsPedigree_v0.5.0/scripts/10_run_full_pipeline.py \
    --delly delly.vcf.gz \
    --manta manta.vcf.gz \
    --list-of-lrr lrr_list.tsv \
    --karyotype-catalogue catalogue.json \
    --outdir results/

# Optional: ngsTracts hand-off (per LRR, separate step)
python ngsPedigree_v0.5.0/scripts/06_build_polarization_input.py \
    --stage1-edges  results/pairwise_relationship_classification.tsv \
    --stage1-roster results/roster.tsv \
    --karyotype     results/karyotype_calls.json \
    --inversion-id  cLRR_0007 \
    --polarity-hint band_0_is_REF \
    --out           polarization_in.json
python ngsPedigree_v0.5.0/scripts/05_polarize_inversion.py \
    --in    polarization_in.json \
    --mtdna mtdna.tsv \
    --out   ngstracts_input.json
```

---

## 6. Outputs explained — what each file means

### 6.1 `merged_del_catalogue.json`

The unified DEL marker set. One entry per merged marker with start,
end, contributing callers, and per-sample genotypes. **Use this as
the canonical DEL marker set for downstream analysis.**

### 6.2 `pairwise_relationship_classification.tsv`

Stage-1-shape table. One row per pair: theta, IBS0, edge_class,
confidence. Filter on `edge_class == 'parent_offspring'` to get
candidate PO pairs. Filter on `edge_class == 'full_sibling'` to get
candidate sibling pairs. This table is drop-in compatible with
`scripts/06_build_polarization_input.py`.

### 6.3 `candidate_PO_pairs.tsv`

The PO subset of the pairwise table with Mendelian-exclusion
exclusion-count columns. Use it to identify your dyads.

### 6.4 `detected_triads.tsv`

Triads found by the forced-offspring rule. **One row per triad** with
paternal, maternal, offspring sample IDs. Note: until you supply sex
info, paternal/maternal is by convention only — flip if necessary.

### 6.5 `inheritance_segments.tsv`

**This is the chromosome inheritance map.** One row per (relationship,
chromosome, segment). Columns:

| Column | Meaning |
|---|---|
| `relationship_id` | the triad/dyad |
| `relationship_type` | `triad_paternal` / `triad_maternal` / `dyad` |
| `parent_sample_id` | who contributed |
| `offspring_sample_id` | who received |
| `chrom`, `seg_start`, `seg_end` | segment span |
| `transmitted_allele` | `REF` / `DEL` / `ambiguous` |
| `n_informative_markers` | parent-HET markers in this segment |
| `confidence` | `Gold` (≥8) / `Silver` (3–7) / `Bronze` (1–2) |
| `recomb_event_left/right` | candidate recombination at this boundary |
| `parental_hap_inherited` | `1` (REF-transmitting hap) / `2` (DEL) / `ambiguous` |

Interpret a row as: "from `parent_sample_id` to `offspring_sample_id`,
the chromosome interval `chrom:seg_start-seg_end` carries DEL markers
that all came from the same parental haplotype (the one that
transmits REF, or DEL, as labelled)."

A change between adjacent rows on the same chromosome is a candidate
crossover.

### 6.6 `discovered_lrrs.tsv` (when `--discover-lrrs`)

The de novo LRR list. Same shape as `--list-of-lrr` consumes, plus
QC columns `n_markers` and `mean_pairwise_correlation`. Higher
correlation = stronger LRR signal.

### 6.7 `lrr_enrichment.tsv` (when an LRR set exists)

One row per LRR with the family-based OR enrichment:

| Column | Meaning |
|---|---|
| `lrr_id`, `chrom`, `start`, `end` | identity |
| `combined_OR`, `combined_CI_low`, `combined_CI_high` | combined OR with 95% Woolf CI |
| `combined_a`, `_b`, `_c`, `_d` | the 2×2 cells |
| `combined_haldane_corrected` | True if any cell was 0 → Haldane (+0.5) applied |
| `triad_OR`, `triad_CI_low`, `triad_CI_high` | triad-only OR |
| `dyad_OR`, `dyad_CI_low`, `dyad_CI_high` | dyad-only OR |

**OR > 1 with tight CI = strong family-based evidence the LRR behaves
as an inherited block.**

### 6.8 `del_arrangement_linkage.tsv` (when `--karyotype-catalogue`)

The situation-1 table. One row per (DEL, LRR):

| Column | Meaning |
|---|---|
| `n_hom_ref`, `n_het`, `n_hom_inv` | sample counts per class |
| `del_freq_hom_ref/_het/_hom_inv` | DEL allele frequency per class |
| `interpretation` | one of `arrangement_0_marker`, `arrangement_1_marker`, `arrangement_0_marker_hom_depleted`, `arrangement_1_marker_hom_depleted`, `unlinked`, `ambiguous` |
| `notes` | the not-lethality disclaimer when hom_depleted |

### 6.9 `lrr_comparative_summary.tsv` + `lrr_comparative_regression.tsv` (bloc 18 output)

The cross-LRR comparative test. Summary table has one row per LRR
with `FIS`, `min_hom_het_ratio`, `dxy_between_arrangements`, etc.
Regression table has one row per fitted model (3 rows: log(min ratio),
FIS, het-excess all regressed against dXY) with slope, r, t-stat,
asymptotic p-value.

---

## 7. How to interpret it for the thesis

The defensible chain of inference, in order:

1. **The pedigree exists.** From `candidate_PO_pairs.tsv` and
   `detected_triads.tsv`, you have first-degree relationships and at
   least some directional triads — recovered from DEL VCFs alone.

2. **Inheritance is Mendelian.** From `inheritance_segments.tsv`,
   offspring chromosomes decompose into segments inherited from each
   parent. Recombination events appear as transmission-allele
   switches.

3. **Candidate LRRs behave as inherited blocks** (cite bloc 13's OR
   enrichment). From `lrr_enrichment.tsv`, the OR > 1 with tight CI
   tells you the LRR transmits as a block, distinct from matched
   background regions.

4. **The cross-type analysis is Mendelian-consistent or surfaces
   distortion** (cite bloc 14). From the per-family per-LRR table,
   crosses behave per the Punnett expectation, with deviations
   flagged as `segregation_distorted` or `homozygote_depletion`.

5. **Per-DEL arrangement linkage reveals situation 1 cases**
   (cite bloc 17). The `arrangement_X_marker_hom_depleted` labels in
   `del_arrangement_linkage.tsv` mark candidate pseudo-overdominance
   regions — **with the explicit not-lethality disclaimer**.

6. **The signal scales with divergence** (cite bloc 18). From
   `lrr_comparative_regression.tsv`, the negative slope of
   `log(min_hom_het_ratio)` against `dxy_between_arrangements`
   supports progressive pseudo-overdominance as a model. If the slope
   isn't significant, the negative result is still publishable —
   alternative phrasing in bloc 18 doc.

**The thesis claim that this chain supports:**

> Candidate LRRs in the hatchery cohort behaved as inherited
> recombination-suppressed haplotype blocks: they were transmitted
> across triads in a Mendelian-compatible manner, were enriched for
> block-compatible inheritance relative to matched background
> regions, and the magnitude of homokaryotype depletion scaled with
> between-arrangement divergence — a pattern consistent with
> progressive pseudo-overdominance under recombination suppression.
> Molecular breakpoint validation is treated as a separate analysis
> stream.

**The thesis claim that this chain does NOT support:**

- Specific molecular breakpoints (needs split-read / assembly).
- Lethality of any specific homokaryotype class.
- Causal age effects on depletion (dXY is a divergence proxy, not
  an age estimate).

---

## 8. Worked synthetic-fixture example

The repo's `synthetic_svvcfs` fixture is a 4-sample synthetic cohort
(P_F, P_M, C, UNREL). After the one-command run:

```
[full-pipeline] merged catalogue: 18 markers (6 two-caller, 12 single)
[full-pipeline] cohort: 4 samples
[full-pipeline] candidate PO edges: 2  detected triads: 1
[full-pipeline] de novo discovered LRRs: 2
[full-pipeline] wrote lrr_enrichment.tsv (2 LRRs)
```

Triad `P_F + P_M → C` detected without being told. Two de novo LRRs
discovered (Chr1 cLRR_0001 with mean |r| 0.64; Chr2 cLRR_0002 with
0.50). Inheritance map shows a candidate Chr1 recombination event
between markers at ~40 kb and ~50 kb in the paternal trace.

On real catfish data with ~200 samples and ~200,000+ DEL markers,
expect:
- thousands of candidate PO pairs at first-degree θ
- tens to hundreds of triads via forced-offspring rule
- ~30–200 candidate LRRs (depending on local recombination
  landscape)
- per-LRR ORs above 1 in the LRRs that are real inversions; ORs
  near 1 in noise calls

---

## 9. Tuning knobs (defaults are conservative)

| Flag | Default | When to change |
|---|---|---|
| `--bp-tolerance` | 500 | longer for noisier breakpoint calls |
| `--reciprocal-overlap` | 0.5 | raise to 0.7+ for tighter merging |
| `--min-informative-pair` | 50 | raise on dense data, lower on sparse |
| `--window-size` (discovery) | 1 Mb | shrink for fine-grained LRRs |
| `--min-markers-per-window` | 4 | lower → more sensitive, more FPs |
| `--correlation-threshold` | 0.5 | raise → tighter LRR calls |
| `--n-background` (enrichment) | 10 | raise for more stable ORs |
| `--min-excluding` (exclusion) | 1 | raise to tolerate 1/1 short-read miscalls |

---

## 10. Limitations & honest caveats

This is a no-budget pipeline. It is not a substitute for proper
local-PCA discovery, BEAGLE-based inheritance mapping, or split-read
breakpoint validation. The trade-offs:

- **DEL marker density:** sparser than SNPs. PO classification still
  robust at 200k+ DELs. LRR discovery may miss low-DEL-density
  inversions.
- **Dyad direction:** undecidable from pure DEL genotypes when both
  directions are Mendelian-compatible. The triad forced-offspring
  rule resolves it when applicable; otherwise the direction is
  reported as `ambiguous`.
- **Inheritance-map resolution:** DEL-marker resolution (typically
  10 kb–10 Mb segments), not SNP resolution. Adequate to seed
  ngsTracts; not a replacement for SNP-level mapping.
- **De novo LRR discovery:** correlation-based, less sensitive than
  local-PCA on SNPs. Use both when possible.
- **Situation 1:** the `hom_depleted` label is a candidate flag, not
  proof of overdominance. The explicit not-lethality disclaimer is
  enforced by unit test.
- **Cross-LRR regression p-values:** asymptotic when n_LRRs < 30;
  compute exact t-tests in R/statsmodels for the final manuscript
  table.

---

## 11. Where to find each piece

```
docs/hpp/pages/
  00_index.md                              ← navigation
  06_inversion_polarization.md             ← polarization layer
  07_mtdna_maternal_validation.md          ← --mtdna pre-flight
  11_hemizygous_markers.md                 ← fake-trio direction
  12_sv_only_pedigree.md                   ← bloc 12 in detail
  13_lrr_enrichment.md                     ← OR enrichment
  14_mendelian_segregation.md              ← per-family per-LRR
  15_framework.md                          ← user's 12-step framework
  16_full_pipeline_and_discovery.md        ← --discover-lrrs
  17_mtdna_methods_de_novo_lrr_discovery.md  ← how de novo works
  18_lrr_divergence_comparative.md         ← cross-LRR test

ngsPedigree_v0.5.0/
  scripts/
    08_pedigree_from_sv_vcfs.py            ← bloc 12 only
    09_lrr_enrichment.py                   ← bloc 13 only
    10_run_full_pipeline.py                ← *** EVERYTHING ***
  src/hpp/
    vcf_sv.py                              ← VCF parsing
    catalogue_merge.py                     ← bloc 12 merging
    del_relatedness.py                     ← KING-robust
    del_inheritance.py                     ← inheritance map
    lrr_discovery.py                       ← bloc 16
    lrr_enrichment.py                      ← bloc 13
    mendelian_segregation.py               ← bloc 14
    del_arrangement_linkage.py             ← bloc 17 (situation 1)
    lrr_divergence_comparative.py          ← bloc 18
    relatedness_sim.py                     ← classifier shadow
```

319 unit tests cover every step. Run them with:

```
python ngsPedigree_v0.5.0/tests/run_tests.py
```

---

## 12. TL;DR

**You have:** Delly + Manta VCFs.
**You run:** one command (section 3).
**You get:** pedigree, inheritance map, candidate LRRs (curated or
de novo), family-based evidence at each LRR, situation-1 depletion
classifier, cross-LRR comparative test.
**You write the thesis with:** section 7's chain of inference and
the honest-but-publishable wording in bloc 18's doc.
**You don't need:** ngsRelate, BEAGLE, HPC, SNP genotype tables, PCA
output, dosage labels, or money.

That's the whole pipeline.
