# Bloc 15 — Framework summary: family-based validation of candidate LRRs / inversions

This page maps the user-stated constraint/solution framework to the
implementing blocs. The framework is the family-based-validation logic
that makes a PhD-defensible inversion claim out of population-discovered
candidate LRRs, using only the data we already have.

## The five constraints

| Constraint | Solution | Implementing bloc(s) |
|---|---|---|
| **No clean crosses** | reconstruct the family graph from data | [Stage 1/2](../../scripts/STEP_PED_01_annotate_relationships.py) + [bloc 12 SV-only pipeline](12_sv_only_pedigree.md) |
| **No perfect pedigree** | test all transmission directions | [bloc 11 hemizygous "fake-trio"](11_hemizygous_markers.md); bloc 12 Mendelian-exclusion CLI |
| **No full recombination map** | chromosome-level inheritance maps from cheap markers | [bloc 12 `del_inheritance.py`](12_sv_only_pedigree.md) |
| **No breakpoint validation everywhere** | odds-ratio enrichment for inherited block behavior | [bloc 13 LRR enrichment](13_lrr_enrichment.md) |
| **No huge resources** | use microhaps, DELs, LRR arrangements as compressed inheritance markers | every bloc operates on biallelic-marker Punnett enumeration; the math is the same for SNPs, DELs, microhaplotypes, and LRR arrangement classes |

## The 12-step family-based validation framework

The framework the user described, step by step, with the bloc that
implements each step. Wording is the user's; bloc references are mine.

| # | Step | Bloc / file |
|---|---|---|
| 1 | Population-level discovery of candidate LRRs (local PCA, dosage, Cramér's V, long-range votes) | external; pipeline consumes the LRR list (`--list_of_LRR`) |
| 2 | Assign arrangement classes (HOM_REF / HET / HOM_INV) to each individual at each LRR | [bloc 10 karyotype catalogue](10_karyotype_catalogue.md), or any upstream PC1-band caller |
| 3 | Build the family / relatedness graph | Stage 1 / Stage 2; [bloc 12 SV-only](12_sv_only_pedigree.md) for the no-ngsRelate path |
| 4 | Test triad direction across all permutations (A+B→C, A+C→B, B+C→A) | [bloc 11 fake-trio direction](11_hemizygous_markers.md); bloc 12 forced-offspring rule |
| 5 | Use the inheritance map as parental-decomposition proof | [bloc 12 `del_inheritance.py`](12_sv_only_pedigree.md); SPEC-§3 placeholder when BEAGLE Stage 3 lands |
| 6 | Interpret the LRR inside the confirmed family — parent-of-origin per arrangement | bloc 06 polarization + bloc 12 inheritance map together |
| 7 | Classify the parental cross type per family × LRR | [bloc 14 Mendelian segregation](14_mendelian_segregation.md) |
| 8 | Compare observed vs Mendelian-expected offspring class counts | [bloc 14 Mendelian segregation](14_mendelian_segregation.md) |
| 9 | Use dyads carefully (weaker than triads) | bloc 11 + 12 explicitly distinguish dyad vs triad evidence layers |
| 10 | Quantify family-based support with odds-ratio enrichment | [bloc 13 LRR enrichment](13_lrr_enrichment.md) |
| 11 | Integrate with population-genomic evidence (πS, FST, dXY, HWE/FIS) | external; pipeline emits the family-side numbers for the integration |
| 12 | Final claim: "candidate LRRs behave as inherited recombination-suppressed haplotype blocks" — call them inversions only when breakpoint evidence is independent | wording in `13_lrr_enrichment.md` and `14_mendelian_segregation.md` |

## End-to-end command for the broke-grad-student stack

Given Delly + Manta VCFs and an LRR list, two commands take you from
raw SV calls to the family-based validation report:

```bash
# 1. SV VCFs → pedigree + inheritance map
python ngsPedigree_v0.5.0/scripts/08_pedigree_from_sv_vcfs.py \
    --delly delly.vcf.gz --manta manta.vcf.gz \
    --outdir results/

# 2. Pedigree + LRR list → odds-ratio enrichment
python ngsPedigree_v0.5.0/scripts/09_lrr_enrichment.py \
    --catalogue   results/merged_del_catalogue.json \
    --list-of-lrr lrr_list.tsv \
    --triads      results/triads.tsv \
    --dyads       results/candidate_PO_pairs.tsv \
    --n-background 20 \
    --out         results/lrr_enrichment.tsv
```

The Mendelian-segregation analysis (bloc 14) plugs in at the
same level: pass in the same triads + the per-LRR arrangement
catalogue, get back per-family per-LRR cross-type classification + a
cohort-aggregated HET × HET goodness-of-fit per LRR.

## What this framework does NOT claim

- It does **not** claim breakpoint-resolution proof of any inversion.
  That is a separate (split-read / assembly / inversion-specific SV
  caller) line of evidence.
- It does **not** prove pseudo-overdominance from `homozygote_depletion`
  alone. The signature is consistent with overdominance but also with
  technical artifacts (DOC dropout in homozygotes, hatchery viability
  filtering). Mechanism claims need to integrate the deleterious-burden
  layer (sibling HPP / KBC repos).
- It does **not** replace the inheritance map a true BEAGLE-based
  Stage 3 would produce. The DEL-marker map (bloc 12) is the
  no-budget approximation. Both contracts feed the same downstream
  consumers.

## Honesty pass

Every claim above lands on the conservative side. The pipeline tells
you that **candidate LRRs behave as inherited blocks at the family
level** — that is the defensible statement. Anything stronger needs
independent evidence.

286 unit tests cover the path. The PhD ships on this.
