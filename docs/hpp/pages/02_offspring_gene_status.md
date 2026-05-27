# Bloc 02 — `offspring_gene_status`

| | |
|---|---|
| **analysis_id** | `offspring_gene_status` |
| **module** | `ngspedigree_hpp` |
| **produces** | `hpp_offspring_gene_status` (Table B) |
| **schema** | `ngsPedigree_v0.5.0/schemas/B_hpp_offspring_gene_status.schema.json` |
| **MVP** | 3 — **built** |

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

## What's built

- `src/hpp/gene_status.classify_gene_status()` — groups Table A rows by
  (offspring, gene, segment), partitions by `hap_copy` + `allele_state`,
  applies the seven-class rule.
- `_classify(d1, d2, unresolved)` — pure rule; all 7 enum cases unit-tested.
- `TableBRow` dataclass — column list verified equal to schema B.
- Confidence aggregation = worst-of over contributing rows (SPEC §5).
  Synonymous-only groups still propagate worst confidence.
- Segment lookup re-derived from `chrom + pos` against the inheritance
  map (segments aren't carried on Table A per spec).
- Real-fixture coverage on dyad and triad fixtures at T1 and T3.

## Open questions

- Whether `compound_het_trans` rows should also get a Tier-4-style
  frame-aware coding consequence applied (apply the variants to CDS
  and translate). Current draft says no — that's HAPS's job
  (sibling spec, post-manuscript). HPP stays at the categorical level.
