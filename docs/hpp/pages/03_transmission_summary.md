# Bloc 03 — `transmission_summary`

| | |
|---|---|
| **analysis_id** | `transmission_summary` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_transmission_summary` (Tables C + D unified) |
| **schemas** | `C_hpp_dyad_transmission_summary.schema.json` and `D_hpp_triad_transmission_summary.schema.json` |
| **MVP** | 4 — **built** |

## Goal

Per-relationship (dyad or triad) rollup: variant-transmission counts,
gene-status counts, Mendelian consistency. `relationship_type` is a
runtime parameter — one bloc, two table shapes.

## Output split (C vs D)

Atlas-core stores the layer as a single `hpp_transmission_summary`,
but the actual TSV is shape-specific:

| Triggered when | Schema | Notes |
|---|---|---|
| dyads | `C_hpp_dyad_transmission_summary` | partial Mendelian — only parent-hom sites are testable |
| triads | `D_hpp_triad_transmission_summary` | full Mendelian at every both-parents-called site + `n_de_novo_candidates` flag |

## Mendelian rules (SPEC §6 step 5)

**Triad — full test.** At every site where both parents have called
genotypes, compute the Mendelian-expected genotype set; if offspring
falls outside, increment `mendelian_inconsistent_sites`. If damaging
→ `mendelian_inconsistent_damaging_sites`. If offspring carries an
allele neither parent carries → `n_de_novo_candidates` (flag, NOT a
claim — could be Stage 3 error or genotyping noise at 9×).

**Dyad — partial.** Only sites where the parent is homozygous are
testable. `mendelian_consistency_status` reaches `pass` / `warn` /
`fail` on that subset; the full-site status is `untestable`.

## Family hubs are not used

Aggregation is strictly per dyad and per triad. The hatchery's
mixed-family hub structure is a confound for cohort-level work, not a
unit of analysis. There is no `family_transmission_summary` table.

## What's built

- `src/hpp/mendelian.py` — primitives:
  - `mendelian_expected_set(gp1, gp2)` Punnett-square enumeration
  - `triad_consistent(gp1, gp2, go)` full Mendelian check
  - `offspring_carries_novel_allele(...)` de-novo flag (NOT a claim)
  - `dyad_partial_consistent(parent_gt, offspring_gt)` parent-hom subset
  - `status_from_counts(...)` → `pass`/`warn`/`fail`/`untestable`
    (default thresholds: 0 → pass, 1–2 → warn, ≥3 → fail; HANDOFF
    open question #5 — calibrate against empirical data later)

- `src/hpp/summary.py` — emits:
  - `summarise_dyad()` → `TableCRow` (19 cols, matches schema C)
  - `summarise_triad()` → `TableDRow` (23 cols, matches schema D)

- Real-fixture coverage:
  - `synthetic_dyad` / `synthetic_triad` exercise the "pass" path.
  - `synthetic_triad_violation` is a deliberate-violation fixture
    covering the `fail` status, `inconsistent_damaging_sites`
    counting under T1 vs T3, and the `n_de_novo_candidates` flag.

## Open questions

- Threshold for downgrading `pass` → `warn` → `fail` on the dyad
  partial-Mendelian subset. Current draft: 0 inconsistencies → pass;
  1–2 → warn; ≥3 → fail. Wants empirical calibration from a small
  test run.
- Whether `n_de_novo_candidates` should be reported at all in the
  manuscript supplementary or only as an internal QC flag.
  Current draft: report with explicit caveat language.
