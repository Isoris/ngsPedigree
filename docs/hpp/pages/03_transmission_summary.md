# Bloc 03 — `transmission_summary`

| | |
|---|---|
| **analysis_id** | `transmission_summary` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_transmission_summary` (Tables C + D unified) |
| **schemas** | `C_hpp_dyad_transmission_summary.schema.json` and `D_hpp_triad_transmission_summary.schema.json` |
| **MVP** | 4 |

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

## What's NOT built yet (MVP 4)

Algorithm in `SPEC_HPP.md` §10.4 pseudocode; output writer hooks in
`src/hpp/io.py`; the summariser will land at `src/hpp/summary.py`
with the Mendelian check at `src/hpp/mendelian.py`.

## Open questions

- Threshold for downgrading `pass` → `warn` → `fail` on the dyad
  partial-Mendelian subset. Current draft: 0 inconsistencies → pass;
  1–2 → warn; ≥3 → fail. Wants empirical calibration from a small
  test run.
- Whether `n_de_novo_candidates` should be reported at all in the
  manuscript supplementary or only as an internal QC flag.
  Current draft: report with explicit caveat language.
