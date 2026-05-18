# Haplotype-aware deleterious-burden specs — umbrella README

**Status:** SPEC ONLY — awaiting audit before implementation.
**Version:** v0.8 (KBC: LoF-first three-tier damaging-variant
definition, splice-validation six-check gate, headline = Tier 1 only)
**Project:** `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort.

---

## What this is

A set of three independent specifications that together cover
haplotype-aware deleterious-burden analysis for the inversions
manuscript. They are **sibling methods, not tiers of one method**.
They share input infrastructure (variant-master table, joint VCF,
reference, GFF) and output schemas where it makes sense, but they
answer different questions, take different inputs, and produce
different outputs.

| | **KBC** | **HPP** | **HAPS** |
|---|---|---|---|
| **File** | `KBC_SPEC.md` | `HPP_SPEC.md` | `HAPS_SPEC.md` |
| **Full name** | Karyotype Burden Contrast | Haplotype Projection from Pedigree | Haplotype-Aware Protein Scoring |
| **Unit of analysis** | inversion × karyotype class × POD segment | offspring × segment × dyad (or triad) | sample × transcript × haplotype |
| **Question** | Does the AA / AB / BB contrast show the masking → exposure pattern that's POD-compatible? | What variants does this specific offspring carry on which specific chromosome copy at this segment, inferred from a parental haplotype projection? | For a missense variant on this specific haplotype, what does VESM say when we feed it the actual protein sequence the haplotype is translating (not the reference)? |
| **Required inputs** | variant-master, PCAngsd K=3 karyotypes, inversion intervals | variant-master, joint VCF, ngsPedigree Stage 3 inheritance maps, parental VCFs | Clair3 phased VCFs, GFF3, reference FASTA, candidate region set, ESM1b/VESM weights |
| **Family structure as input** | not used — corrected for as relatedness confound | only via explicit dyads / triads from ngsPedigree | not relevant |
| **Phasing required** | no — karyotype call provides haplotype background for AA / BB | yes, via inheritance map | yes — Clair3 read-backed (within phase blocks) |
| **Manuscript role** | **headline POD result** — feeds v20 directly | **mechanism / case studies / supplementary** — where Stage 3 supports it | **quality refinement** — supplementary paragraph + table listing variants where haplotype-aware VESM differs from reference-frame |
| **Blocks on** | nothing — runs on inputs already in hand | ngsPedigree Stage 3 readiness | nothing in principle; **deferred until after manuscript by author decision** |
| **Ship order** | first, in time for v20 | after KBC, once Stage 3 lands | after manuscript |
| **Critical path for v20** | yes | no | no |

## What they don't do (in common)

All three specs explicitly do not:

- re-run SnpEff or SIFT4G (HAPS does add a new VESM column, but does not replace the existing one);
- depend on GERP / phastCons / phyloP / orthology-tier / CAFE / family-evolution scores (the EGO pipeline is decoupled);
- infer balancing selection, true overdominance, or fitness — all three report structural / mutational distributions only;
- conflate the three catfish cohorts — all three operate exclusively on the 226-sample pure *C. gariepinus* hatchery cohort;
- use HWE deviation as evidence for anything;
- claim to solve every problem in this field — there are many open issues (in-frame indel calibration, isoform-aware scoring, statistical phase quality at 9×, etc.) that these specs deliberately bound around rather than solve.

## How the three fit together

```
                Clair3 phased VCFs ─────────┐
                                            │
                MODULE_CONSERVATION ────────┼────▶ variant_master_scored.tsv
                (SnpEff, SIFT4G, VESM)      │            │
                                            │            │ adds new column
                                            ▼            │ vesm_llr_haplotype_aware
                                          HAPS ──────────┤ (when built, post-manuscript)
                                                         │
                                                         ▼
                                          ┌────────────────────────────┐
                                          │                            │
                                          ▼                            ▼
   PCAngsd K=3 karyotypes ─────────▶    KBC                          HPP ◀─── ngsPedigree
   POD segments                          (headline)                   (mechanism /        Stage 3 maps
                                                                     supplementary)
```

KBC consumes `variant_master_scored.tsv` whether or not HAPS has run.
If HAPS later upgrades the VESM column for a subset of variants, KBC
reruns and burden numbers shift by a small amount — the
POD-compatible headline pattern is robust to that magnitude of
change.

## Repo home (proposed)

```
MODULE_CONSERVATION/                # existing
  haps/                             # new submodule — ships post-manuscript
    SPEC_HAPS.md
    ...

catfish-variant-analysis/Modules/   # existing
  NN_kbc/                           # ships first — v20
    SPEC_KBC.md
    ...
  NN+1_hpp/                         # ships after Stage 3 lands
    SPEC_HPP.md
    ...
  shared_io/                        # shared loaders
    ...
```

HAPS sits inside MODULE_CONSERVATION because it modifies the
canonical variant table (`variant_master_scored.tsv`). KBC and HPP
sit inside `catfish-variant-analysis` because they consume that
table.

## Cross-references

- `MODULE_CONSERVATION` produces `variant_master_scored.tsv`. KBC and HPP consume it. HAPS extends it with a new column.
- `ngsPedigree` Stage 1 / 2 / 3 produces dyad lists and inheritance maps for HPP. Stage 3 is not yet ready — HPP's spec contains a placeholder schema so the module can be implemented before Stage 3 lands.
- PCAngsd K=3 + Hungarian arrangement-label matching produces the karyotype calls KBC consumes.
- The inversion atlas defines the candidate intervals and the L / M / R POD partition.
- The splice module's `SPLICE_SUBCLASS` codes are read by all three specs where relevant.
- bcftools/csq is the haplotype-aware consequence caller that HAPS wraps; it is also referenced indirectly by KBC and HPP through MODULE_CONSERVATION's existing CSQ cross-check column.

## Changelog

- **v0.8** — KBC §1.8 added: three-tier damaging-variant definition
  with high-confidence LoF (stop-gained, frameshift, start-lost,
  canonical splice donor/acceptor) as **Tier 1 and the manuscript
  headline**; Tier 2 adds validated splice Class A from the splice
  module (gated by §1.9 splice validation checklist); Tier 3 adds
  SIFT4G-deleterious and VESM-strong missense as exploratory only.
  KBC §1.9 added: six-check splice validation gate (canonical sites,
  GFF coordinates, strand handling, boundary classification, SnpEff
  concordance, isoform stability). Table E restructured to carry
  per-tier signal classes and a `tier_concordance` diagnostic.
  Manuscript paragraph rewritten so the headline always uses Tier 1
  only.
- **v0.7** — KBC §1.6 added: each inversion is one independent test
  unit; across-inversion summaries are descriptive (signal-class
  counts) or meta-analytic (inverse-variance weighted), never
  sum-of-raw-counts. KBC §1.7 added: synonymous variants are a
  parallel negative control, not a dN/dS-style denominator; the
  headline statistic is not a ratio. New output table F:
  `kbc_across_inversion_summary.tsv` with per-class counts and a
  random-effects meta-analytic sheet. §12 (out-of-scope) updated to
  exclude pooled raw-count tests and ratio-form headline statistics
  explicitly.
- **v0.6** — added HAPS_SPEC.md (haplotype-aware protein
  reconstruction + VESM rescoring for candidate regions). Updated
  README to three siblings.
- **v0.5** — split IHC into KBC and HPP siblings; dropped the "tiered
  confidence" framing that conflated cohort-level burden contrast
  with individual-level inheritance projection.
- **v0.4** — recessive HET-masking biology made explicit; per-gene
  het-masked count restored as carrier-load measurement; F-table
  renamed `karyotype_burden_summary.tsv`.
- **v0.3** — per-karyotype-class burden (AA / AB / BB); AB as masked
  reference; `exposure_delta_vs_AB` headline.
- **v0.2** — dropped EGO / GERP dependencies.
- **v0.1** — initial draft consolidating Pastes 1 + 2 + 3.
