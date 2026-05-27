# Bloc 02 — `offspring_gene_status`

| | |
|---|---|
| **analysis_id** | `offspring_gene_status` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_offspring_gene_status` (Table B) |
| **schema** | `ngsPedigree_v0.5.0/schemas/B_hpp_offspring_gene_status.schema.json` |
| **MVP** | 3 |

## Goal

For each (offspring, gene, segment), classify into the seven-class
status enum based on the per-haplotype damaging variants from Table A.
This is the **manuscript-supplementary headline** product.

## Inputs

| Dimension | Adapter |
|---|---|
| `hpp_offspring_haplotype_variants` | output of Bloc 01 (Table A) |
| `variant_master_scored` | `variant_master.PlaceholderVariantMaster` → real loader at MVP 2 |

Damaging-variant tier (`T1` / `T2` / `T3` per KBC §1.8) is a **runtime
parameter**. Same data table, three classifications.

## Classification rules (SPEC §7)

```python
d1 = damaging variants on hap_from_P1 in (gene ∩ segment)
d2 = damaging variants on hap_from_P2 in (gene ∩ segment)
unresolved = damaging variants in (gene ∩ segment) with hap_copy = 'unassigned'

if unresolved:
    return 'partially_resolved' if d1 or d2 else 'unresolved'
if not d1 and not d2:
    return 'reference_like'
if d1 and not d2:
    return 'het_masked' if len(d1) == 1 else 'compound_het_cis'
if d2 and not d1:
    return 'het_masked' if len(d2) == 1 else 'compound_het_cis'
# both copies hit
if d1 == d2 and len(d1) == 1:
    return 'hom_exposed_same_variant'
return 'compound_het_trans'
```

The seven categories — same as KBC, but here resolved per-individual
with inheritance-map exactness rather than karyotype-class average.

## Confidence

Per gene-status row: the **minimum** confidence over all per-variant
rows feeding into either copy. Conservative by design — a single
Bronze segment downgrades the entire gene classification.

## What's NOT built yet (MVP 3)

Everything. Algorithm is in `SPEC_HPP.md` §10.3 pseudocode; adapter
hooks exist in `src/hpp/io.py`; the classifier function will land at
`src/hpp/gene_status.py`.

## Open questions

- Whether `compound_het_trans` rows should also get a Tier-4-style
  frame-aware coding consequence applied (apply the variants to CDS
  and translate). Current draft says no — that's HAPS's job
  (sibling spec, post-manuscript). HPP stays at the categorical level.
