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
| 07 | [`07_mtdna_maternal_validation.md`](07_mtdna_maternal_validation.md) | `mtdna_maternal_validation` (pedigree pre-flight) | `mtdna_validation` block in OUT JSON |
| 08 | [`08_backbone_verification.md`](08_backbone_verification.md) | full backbone (ngsRelate `.res` → ngsTracts), verified | — |
| 09 | [`09_synthetic_panels.md`](09_synthetic_panels.md) | synthetic-panel recovery tests (6 topology mixes) | recovery report |
| 10 | [`10_karyotype_catalogue.md`](10_karyotype_catalogue.md) | karyotype-catalogue JSON adapter (registry IN) | `KaryotypeCall` list |
| 11 | [`11_hemizygous_markers.md`](11_hemizygous_markers.md) | hemizygous DEL markers ("fake trio" direction) | `TriadVerdict` |
| 12 | [`12_sv_only_pedigree.md`](12_sv_only_pedigree.md) | broke-grad-student pipeline (Delly + Manta → pedigree + chromosome inheritance map) | `inheritance_segments.tsv` |
| 13 | [`13_lrr_enrichment.md`](13_lrr_enrichment.md) | family-based OR enrichment for candidate LRRs | `lrr_enrichment.tsv` |
| 14 | [`14_mendelian_segregation.md`](14_mendelian_segregation.md) | per-family × per-LRR Mendelian segregation (steps 7–8) | `family_lrr.tsv`, `lrr_summary.tsv` |
| 15 | [`15_framework.md`](15_framework.md) | framework summary: 12-step validation chain → blocs | — |
| 16/17 | [`16_full_pipeline_and_discovery.md`](16_full_pipeline_and_discovery.md) | one-command pipeline + de novo LRR discovery + arrangement-linkage classifier (situation 1) | end-to-end TSV/JSON bundle |
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
