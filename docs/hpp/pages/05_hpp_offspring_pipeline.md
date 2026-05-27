# Bloc 05 — `hpp_offspring_pipeline` (CHAIN)

| | |
|---|---|
| **analysis_id** | `hpp_offspring_pipeline` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_offspring_gene_status` (headline) |
| **type** | CHAIN — composes atomic blocs 01–04 |
| **MVP** | full pipeline ships when MVP 1–4 are complete; MVP 5 cross-check is optional |

## Goal

User-runnable end-to-end HPP run. Atlas-core surfaces this as the
"big green button" for the relatedness_atlas / HPP brick.

## Chain

```
ngsPedigree Stage 2  ──▶ pedigree_relationship_set (dyad/triad table)
ngsPedigree Stage 3  ──▶ inheritance_map_set + parent_phase_map
MODULE_CONSERVATION  ──▶ joint_vcf + variant_master_scored

         │
         ▼
  Bloc 01: haplotype_projection
         │  produces hpp_offspring_haplotype_variants  (Table A)
         ▼
  Bloc 02: offspring_gene_status
         │  produces hpp_offspring_gene_status  (Table B — HEADLINE)
         │
         ├──▶ Bloc 03: transmission_summary (optional terminal)
         │         produces hpp_transmission_summary  (Tables C/D)
         │
         └──▶ Bloc 04: hpp_kbc_arrangement_crosscheck (optional, needs KBC)
                   produces hpp_kbc_crosscheck  (Table E)
```

## Runtime parameters

| Param | Type | Notes |
|---|---|---|
| `damaging_tier` | enum `T1` / `T2` / `T3` | KBC §1.8 tier; default `T1` |
| `relationship_type` | enum `dyad` / `triad` / `both` | default `both` |
| `bronze_policy` | enum `include` / `exclude` | default `include` |
| `candidate_intervals` | BED or `null` | default `null` (genome-wide); accepts inversion-atlas intervals |
| `with_transmission_summary` | bool | default `true` |
| `with_kbc_crosscheck` | bool | default `false` (KBC must be available) |

## Hand-off contract with atlas-core

The chain's `produces` is single-valued: `hpp_offspring_gene_status`.
Optional terminal blocs emit their own layers (`hpp_transmission_summary`,
`hpp_kbc_crosscheck`) which atlas-core surfaces as **derivatives** of
this chain, not as separate chain outputs.

## Cohort scope

226-sample pure *C. gariepinus* hatchery. No cross-species rows.
Restricted to offspring with a Stage 2 confirmed dyad / triad —
typically a subset of the 226.

## What's NOT built yet

The chain orchestrator lands when MVPs 2 + 3 + 4 are complete. Until
then, blocs are runnable individually but the wired-up
`hpp_offspring_pipeline` script does not exist.
