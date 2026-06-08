# HPP — Haplotype Projection from Pedigree

**Atlas:** `relatedness_atlas`
**Module:** `ngspedigree_hpp` v0.5.0-mvp1
**Status:** SPEC + MVP 1 only. Awaiting (a) Stage 3 inheritance-map
schema lock and (b) audit of `SPEC_HPP.md`.

Stage 4 of ngsPedigree. Consumes Stage 3 per-dyad / per-triad
inheritance maps and projects parental consequence-annotated variants
onto offspring haplotype copies. Per-individual structural-mutational
state at the limits of the inheritance map's confidence.

## Pages

| # | Page | Bloc | Output layer |
|---|---|---|---|
| 01 | [`01_haplotype_projection.md`](01_haplotype_projection.md) | `haplotype_projection` | `hpp_offspring_haplotype_variants` |
| 02 | [`02_offspring_gene_status.md`](02_offspring_gene_status.md) | `offspring_gene_status` | `hpp_offspring_gene_status` |
| 03 | [`03_transmission_summary.md`](03_transmission_summary.md) | `transmission_summary` | `hpp_transmission_summary` |
| 04 | [`04_hpp_kbc_arrangement_crosscheck.md`](04_hpp_kbc_arrangement_crosscheck.md) | `hpp_kbc_arrangement_crosscheck` | `hpp_kbc_crosscheck` |
| 05 | [`05_hpp_offspring_pipeline.md`](05_hpp_offspring_pipeline.md) | `hpp_offspring_pipeline` (CHAIN) | `hpp_offspring_gene_status` |
| 06 | [`06_inversion_polarization.md`](06_inversion_polarization.md) | `inversion_polarization` (ngsTracts hand-off) | `polarized_transmissions.out.json` |
| 10 | [`10_inputs_and_adapters.md`](10_inputs_and_adapters.md) | — | adapter contracts |
| 20 | [`20_output_schemas.md`](20_output_schemas.md) | — | JSON Schemas A–E |
| 99 | [`99_status.md`](99_status.md) | — | MVP roadmap + open questions |

## Atlas-core registry

The four JSONL blocks that surface these bricks in atlas-core's
Catalogue (page 4) are in `docs/hpp/atlas_core_registry/` and packaged
as `docs/hpp/atlas_core_registry.tar.gz` with the
`toolkit_registries/relatedness/01_registry/` target path baked in.

## Hard rules (from `HANDOFF.md`)

1. **Three-cohort rule absolute.** 226-sample pure *C. gariepinus* hatchery only.
2. **No family-hub aggregation.** Dyad / triad only.
3. **No cohort-level statistical claim.** That's KBC.
4. **No claim of fitness, balancing selection, or overdominance.**
5. **Confidence labels mandatory** on every assignment.
6. **Stage 3 placeholder schema is a contract.** Swap, don't change.
7. **No HWE-as-evidence.**

## Sibling specs

- `KBC_SPEC.md` — Karyotype Burden Contrast (cohort-level POD test, headline manuscript result)
- `HAPS_SPEC.md` — Haplotype-Aware Protein Scoring (post-manuscript)
- `UMBRELLA_README.md` — how the three siblings fit together
