# Bloc 01 — `haplotype_projection`

| | |
|---|---|
| **analysis_id** | `haplotype_projection` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_offspring_haplotype_variants` (Table A) |
| **schema** | `ngsPedigree_v0.5.0/schemas/A_hpp_offspring_haplotype_variants.schema.json` |
| **MVP** | 2 (MVP 1 has the `parental_haplotype` substep) |

## Goal

For each ngsPedigree-confirmed dyad (P → O) or triad (P1 + P2 → O),
project the parents' consequence-annotated variants onto the
offspring's chromosome copies using the Stage 3 inheritance map.

## Inputs

| Dimension | Adapter | Source |
|---|---|---|
| `pedigree_relationship` | (Stage 2 reader) | ngsPedigree Stage 2 dyad/triad table |
| `inheritance_map` | `stage3_placeholder` → `stage3_real` (MVP 6 swap) | Stage 3 per-relationship map |
| `joint_vcf_genotypes` | `vcf_lite` (synthetic), real adapter (MVP 2) | MODULE_CONSERVATION STEP 03 |
| `parent_phase_map` | `stage3_placeholder.load_parent_phase` | sidecar placeholder; will move into Stage 3 output when locked |

`relationship_type` (`dyad` | `triad`) is a **runtime parameter**, not
a separate bloc. Family hubs are never the input — dyad / triad only.

## Algorithm

1. Per parent: build `(hap_1, hap_2, unphased)` variant lists.
   - Hom-alt (`1/1`) → on **both** haps.
   - Het (`0/1`) with `parent_phase` = 1 or 2 → assigned to that hap.
   - Het with no `parent_phase` entry → `unphased` bucket.
   - Hom-ref (`0/0`) / missing (`./.`) → skipped.
2. For each segment in the inheritance map, look up the
   `parental_hap_inherited` value (`1` / `2` / `ambiguous`).
3. Emit one Table A row per (offspring, variant, hap_copy):
   - if parent is hom-alt at that site → unambiguous, `projection_source = parent_homozygous`.
   - if parent is het and segment phase is known → `parent_heterozygous_phased`.
   - if parent is het and segment is `ambiguous` or variant is in `unphased` bucket → `parent_heterozygous_unphased`, `hap_copy = unassigned`.
4. Triads run the projection from both parents and tag rows by hap origin.

## Confidence (SPEC §5)

```
if segment_confidence == 'Bronze':
    confidence = 'low'
else:
    if projection_source == 'parent_homozygous':
        confidence = 'high'
    elif projection_source == 'parent_heterozygous_phased':
        confidence = 'high' if segment_confidence == 'Gold' else 'medium'
    else:                                                       # unphased
        confidence = 'unresolved'
```

## What's built

- `src/hpp/parental_haps.build_parental_hap_variants()` — step 1 (MVP 1).
- `src/hpp/project.project_dyad_to_offspring()` and `project_triad_to_offspring()` — steps 2–4 (MVP 2).
- `composite_confidence()` — SPEC §5 rule, fully tested.
- `TableARow` dataclass — column list verified to match `A_hpp_offspring_haplotype_variants.schema.json`.
- TSV emission via `io.write_tsv` on `TABLE_A_COLUMNS`.
- Synthetic dyad (4 expected rows at T1) + synthetic triad (5 expected rows at T1) fixtures with EXPECTED_TABLE_A.md known-truth.

## What's NOT built yet (MVP 2b)

- Real joint-VCF adapter (`vcf_lite` remains synthetic-fixture only).
- Dyad's unknown-parent side (offspring's hap-from-other-parent) — currently left unrepresented, not deduced from offspring GT.

## Open questions

- **#1** Where does parent-het phase live in Stage 3's output? MVP 1
  parks this in a sidecar `parent_phase.tsv`; real schema TBD.
- **#2** Bronze-segment default policy (include but tag, or exclude
  from headline counts?). Current config default: `include`.
