# Atlas-core registry hand-off — HPP / `relatedness_atlas`

Four JSONL files to be dropped into atlas-core's
`toolkit_registries/relatedness/01_registry/`:

- `module_registry.jsonl` — 1 row (`ngspedigree_hpp`)
- `analysis_registry.jsonl` — 5 rows (4 atomic + 1 CHAIN)
- `analysis_modes.jsonl` — 5 rows (one per bloc, no fan-out by scope)
- `layer_registry.jsonl` — 4 output layers

Tarball with the atlas-core target path baked in:
`docs/hpp/atlas_core_registry.tar.gz`.

```
$ tar -tzf docs/hpp/atlas_core_registry.tar.gz
toolkit_registries/relatedness/01_registry/module_registry.jsonl
toolkit_registries/relatedness/01_registry/analysis_registry.jsonl
toolkit_registries/relatedness/01_registry/analysis_modes.jsonl
toolkit_registries/relatedness/01_registry/layer_registry.jsonl
```

## Smoke-test invariants (atlas-core enforces)

- every `analysis_modes.analysis_type` ∈ `analysis_registry.analysis_id` ✓
- every `analysis_modes.produces` is single-valued AND ∈ that registry row's declared `produces` ✓
- every `analysis_modes.module_name` ∈ `module_registry.module_name` ✓

## Atomic blocs

| analysis_id | produces (layer) |
|---|---|
| `haplotype_projection` | `hpp_offspring_haplotype_variants` |
| `offspring_gene_status` | `hpp_offspring_gene_status` |
| `transmission_summary` | `hpp_transmission_summary` |
| `hpp_kbc_arrangement_crosscheck` | `hpp_kbc_crosscheck` |

`relationship_type` (dyad | triad) is a runtime parameter on each — NOT
fanned out as separate registry rows, per the unified-ancestry pattern.

## Chain

| analysis_id | produces (headline) |
|---|---|
| `hpp_offspring_pipeline` | `hpp_offspring_gene_status` |

Chains: Stage 2 pedigree set → Stage 3 inheritance maps → projection →
gene status. Optional terminal blocs: `transmission_summary`,
`hpp_kbc_arrangement_crosscheck`.

## Cohort

226-sample pure *Clarias gariepinus* hatchery, ref `fClaHyb_Gar_LG`.
No cross-species rows.

## Backing module status

`ngspedigree_hpp` v0.5.0-mvp1: SPEC + MVP 1 only.
`installed: true`, `ready: false`. `stale_reason:
awaiting_stage3_real_schema_and_spec_audit`. MVPs 2–6 are forward-specced
but not implemented; registry rows are `status: experimental` so
atlas-core surfaces them as planned bricks.
