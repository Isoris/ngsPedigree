# HANDOFF — HPP implementation into `ngsPedigree`

**For:** Claude Code session in the `ngsPedigree` repo.
**You are receiving:** `HPP_SPEC.md` (the spec), `UMBRELLA_README.md` (where HPP sits in the four-sibling stack).
**Status:** SPEC ONLY — **awaiting audit + ngsPedigree Stage 3 readiness before implementation.**
**Manuscript role:** **mechanism / case studies / supplementary** — adds per-individual haplotype-resolved deleterious-variant assignments where Stage 3 supports them.

---

## What HPP is, in one paragraph

HPP (Haplotype Projection from Pedigree) is the individual-level companion to KBC. For each ngsPedigree-confirmed dyad (parent → offspring) or triad (P1 + P2 → O), it uses the Stage 3 inheritance map to project the parents' consequence-annotated variants onto the offspring's chromosome copies. The output is, per offspring × segment, a haplotype-resolved variant list — "this variant came from this parent's hap-1 at this segment, and on the offspring it lives on the copy inherited from that parent." It does not run on the cohort as a unit. It runs dyad-by-dyad and triad-by-triad. Family hubs are not units of analysis — only explicit dyads / triads from Stage 2.

## Where HPP lives

Proposed:

```
ngsPedigree/
  ngsPedigree_v0.3.0/         ← current (Stage 1 + Stage 2 + Relatedness Atlas pages/)
  ngsPedigree_v0.4.0/         ← NEW: Stage 3 inheritance map (already on the roadmap)
  ngsPedigree_v0.5.0/         ← NEW: Stage 4 = HPP (this spec)
    SPEC_HPP.md
    scripts/
    src/hpp/
    schemas/
    tests/
```

HPP is conceptually **Stage 4** of ngsPedigree — it consumes Stage 3's inheritance maps and produces per-offspring haplotype-resolved consequence assignments. The Relatedness Atlas browser pages (already built at `pages/Relatedness_atlas.html`) get a new HPP tab when this lands.

Results land under `${BASE}/results/ngsPedigree/04_hpp/` per the project results-layout convention.

## Blocking dependency

**HPP requires ngsPedigree Stage 3.** Stage 3 builds the chromosome inheritance map from Beagle dosages for confirmed PO dyads — that's the substrate HPP projects onto.

Stage 3 is on the ngsPedigree roadmap (v0.4.0) but not yet built. HPP's spec contains a **placeholder schema** for the Stage 3 output (§3.2 of the spec). HPP can be implemented now against the placeholder for unit tests; the real implementation switches when Stage 3 ships.

## Inputs HPP consumes

| Input | Source | Status |
|---|---|---|
| `variant_master_scored.tsv` | MODULE_CONSERVATION STEP 16 (in `catfish-variant-analysis`) | exists |
| Joint VCF | MODULE_CONSERVATION STEP 03 | exists |
| ngsPedigree Stage 2 dyad/triad table | this repo, Stage 2 already ships | exists |
| ngsPedigree Stage 3 inheritance maps | this repo, Stage 3 = v0.4.0 | **does not yet exist** |
| Reference FASTA | `fClaHyb_Gar_LG.fa` | exists |
| Sample metadata | project sample sheet | exists |

Optional inputs: KBC table B (`kbc_variant_arrangement_assignments.tsv`) — when KBC has shipped, HPP can cross-check inversion-interval variants against the karyotype-derived arrangement assignment.

## Hard rules — DO NOT violate

1. **Three-cohort rule absolute.** 226-sample pure *C. gariepinus* hatchery only.
2. **No family-hub aggregation.** ngsRelate first-degree edge clustering produces hubs that often contain multiple unrelated parental pairs. A hub is NOT a pedigree family. HPP works at the dyad / triad level only.
3. **No cohort-level statistical claim.** That's KBC's job. HPP produces per-individual structural-mutational facts, not population-level inference.
4. **No claim of fitness, balancing selection, or overdominance.** Per-individual mechanism only.
5. **Confidence labels are mandatory.** Every per-variant assignment carries one of: `exact_phased`, `pedigree_supported` (Gold), `family_supported` (Silver), `statistical_phase` (Bronze), `pseudo_phase`, `unphased`. Bronze segments propagate Bronze confidence.
6. **Stage 3 placeholder schema is a contract.** When Stage 3 ships, the placeholder swaps for the real loader without changing downstream contracts. Do not change the schema casually.
7. **No HWE-as-evidence.**

## Three open questions to settle BEFORE coding

1. **Stage 3 inheritance-map schema** — the placeholder in §3.2 of the spec needs to match what Stage 3 will actually emit. In particular, where does parent-heterozygous-site phase live in the Stage 3 output? This is the single biggest implementation risk.
2. **Bronze-segment policy** — include or exclude from headline counts? Current default: include but tag confidence. The audit chat should confirm.
3. **KBC cross-check timing** — does HPP wait for KBC, or do they develop in parallel? Cross-check is optional, so HPP can ship before KBC is fully wired.

## Staged MVP

- **MVP 1** — Stage 3 placeholder loader + parental haplotype builder. Synthetic dyad fixture with known truth.
- **MVP 2** — Dyad projection + Table A (per offspring × variant × hap-copy).
- **MVP 3** — Per-gene status + Table B (using the same gene-status enum as KBC: het_masked, hom_exposed_same_variant, compound_het_trans, compound_het_cis, partially_resolved, unresolved).
- **MVP 4** — Triad projection + Mendelian check + Tables C, D.
- **MVP 5** — KBC cross-check + Table E (once KBC has shipped).
- **MVP 6** — Switch placeholder → real Stage 3 (once Stage 3 ships).

## Integration with the Relatedness Atlas

The existing `pages/Relatedness_atlas.html` has a Mendelian sub-tab and a Family/Individual Evidence Hub. HPP outputs should feed a **new tab** in that atlas:

- per-offspring view: which dyad/triad they belong to, segment-by-segment inheritance map (Gold/Silver/Bronze), per-gene status colour-coded;
- KBC cross-check view: concordance between HPP's haplotype assignment and KBC's arrangement assignment, per offspring per inversion;
- triad Mendelian-consistency view: pass/warn/fail per segment, de novo candidate flag list.

This is the "case studies" view the manuscript supplementary cites.

## What NOT to do

- Do not start building before Stage 3 schema is locked (or at least placeholder-locked with project owner approval).
- Do not aggregate by family hub. Dyad / triad only.
- Do not try to re-derive Stage 3 functionality inside HPP. HPP consumes; Stage 3 produces.
- Do not extend HPP to be a cohort-level burden test. That's KBC's job (sibling spec, lives in `catfish-variant-analysis`).
- Do not promote HPP results to the manuscript headline. The headline is KBC.

## First message I would send

> Read `HPP_SPEC.md` in full. Then read this handoff's three open questions. Coordinate with the ngsPedigree Stage 3 design — confirm the placeholder schema in §3.2 matches Stage 3's planned output. If Stage 3 is not yet specced, start there before HPP. Then build MVP 1 (placeholder loader + parental hap builder + synthetic fixture) and stop for review.
