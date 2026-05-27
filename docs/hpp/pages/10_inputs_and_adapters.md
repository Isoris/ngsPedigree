# Inputs and adapters

HPP is built around swap-target adapters so the projection / gene-status /
Mendelian / cross-check logic does not change when an upstream
producer evolves. Each adapter has two implementations:
**placeholder** (used today, satisfies the contract on synthetic
fixtures) and **real** (swapped in when the upstream producer ships).

## Adapter surface

| Adapter | Placeholder | Real | Status |
|---|---|---|---|
| Stage 3 inheritance maps | `src/hpp/stage3_placeholder.py` | `src/hpp/stage3_real.py` | placeholder built; real stub raises `Stage3RealNotReadyError` |
| Variant annotation | `src/hpp/variant_master.PlaceholderVariantMaster` | `src/hpp/variant_master.load_variant_master` | placeholder built; real loader raises `VariantMasterNotReadyError` (MVP 2) |
| Joint VCF | `src/hpp/vcf_lite.read_vcf` | cyvcf2 / pysam adapter | fixture loader built; real loader = MVP 2 |
| KBC arrangement | `src/hpp/kbc_adapter.PlaceholderKbcAdapter` | `src/hpp/kbc_adapter.load_kbc_table_b` | placeholder built; real loader raises `KbcNotReadyError` (MVP 5) |

`src/hpp/io.py` carries the `Protocol` definitions every adapter
must satisfy and the `default_*_adapter()` wiring HPP uses at runtime.

## Input contracts (placeholder schemas)

| TSV | Schema | Notes |
|---|---|---|
| `inheritance_map_dyad.tsv` | `schemas/inheritance_map_dyad.placeholder.schema.json` | SPEC §3.2 |
| `inheritance_map_triad.tsv` | `schemas/inheritance_map_triad.placeholder.schema.json` | SPEC §3.2 |
| `parent_phase.tsv` | `schemas/parent_phase.placeholder.schema.json` | sidecar for HANDOFF open question #1 |

When ngsPedigree Stage 3 ships, the parent_phase sidecar is expected
to fold into Stage 3's primary output. The placeholder schema is
the contract until that swap.

## Variant_master fields HPP reads

Only these columns of `variant_master_scored.tsv` are touched:
`variant_id`, `gene_id`, `transcript_id`, `consequence`, `impact`,
`sift_class`, `vesm_llr`, `splice_subclass`.

No dependency on GERP / phastCons / phyloP / orthology-tier / CAFE.
EGO pipeline is decoupled.

## KBC fields HPP reads

Only these columns of `kbc_variant_arrangement_assignments.tsv` are
touched: `variant_id`, `inversion_id`, `pod_segment`,
`arrangement_background`, `assignment_confidence`.

## Sample-level inputs

| | Source |
|---|---|
| Sample metadata | project sample sheet (cohort flag, NAToRA pruning flag for KBC's primary mode) |
| ngsRelate kinship | not used by HPP — KBC's secondary mode only |
| ROH BEDs | optional; sanity-checks unambiguous projection in ROH-resident segments |

## Reference

`fClaHyb_Gar_LG.fa` — used only for variant normalisation sanity
checks in MVP 2+.

## Three-cohort rule

HPP operates exclusively on the 226-sample pure *C. gariepinus*
hatchery cohort. Adapters that touch joint-VCF data filter the
sample set to this cohort at the boundary.
