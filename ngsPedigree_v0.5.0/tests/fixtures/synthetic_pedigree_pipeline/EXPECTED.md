# synthetic_pedigree_pipeline — known truth

End-to-end Stage 1/2 → polarization-IN-JSON conversion fixture.

## Setup

- **Hub H001** (`two_parents_with_sibship`, sex-known):
  S001 (father, band 0) + S002 (mother, band 2) → S003 (band 1).
- **Hub H002** (`two_parents_with_sibship`, blind mode):
  S004 (parent_a, band 1) + S005 (parent_b, band 2) →
  S006 (band 1), S007 (band 2), S008 (band 0).

Stage 2 flags `S004↔S007` for review.

## Expected dyads (from PO edges with directional roles)

Without Stage 2 filter — 10 dyads:

| parent | offspring | parent_sex |
|---|---|---|
| S001 | S003 | male |
| S002 | S003 | female |
| S004 | S006 | None (parent_a, sex unknown) |
| S005 | S006 | None |
| S004 | S007 | None |
| S005 | S007 | None |
| S004 | S008 | None |
| S005 | S008 | None |

(8 dyads — S006/S007/S008 each have 2 parents.)

Wait: that's 2 + 6 = 8 dyads. Sibling FS edges (S006↔S007, S006↔S008,
S007↔S008) are not PO so they don't enter dyads. Unrelated pairs are
not PO either. Total: 8 dyads.

With Stage 2 filter (`--stage2-edges`), S004→S007 is dropped → 7 dyads.

## Expected triads

H001: 1 triad (S001, S002, S003) with explicit sex.
H002: 3 triads (S004, S005, S006), (S004, S005, S007), (S004, S005, S008)
with parent_a→paternal, parent_b→maternal convention; a warning is
emitted that sex is unknown.

Without Stage 2 filter: 4 triads.
With Stage 2 filter: the S004↔S007 edge is dropped from the PO
neighbour set, so the triad (S004, S005, S007) loses paternal support
and is **dropped** — 3 triads remain.

## Karyotype calls

The TSV has an `inversion_id` column. With `--inversion-id inv_LG01_pod`,
S001..S008 are included (8 calls); S099 (inv_OTHER) is filtered out.

## Resulting polarization IN JSON

- `schema = "ngspedigree_karyotype_calls_in_v1"`
- `inversion_id = "inv_LG01_pod"`
- `polarity_hint = "band_0_is_REF"` (as supplied)
- `karyotype_calls`: 8 entries
- `dyads`: 8 (without --stage2) or 7 (with --stage2)
- `triads`: 4 (without --stage2) or 3 (with --stage2)
- Warnings include "sex unknown for the two parents" for H002.

## End-to-end (without --stage2, with --mtdna)

Run `06_build_polarization_input.py` → IN JSON, then run
`05_polarize_inversion.py --in <that> --mtdna ...`. This is the
verified full pipeline:

```
ngsRelate .res  →  Stage 1 (built)
                →  Stage 2 (built, optional QC filter)
                →  06_build_polarization_input.py  (built — closes the gap)
                →  IN JSON
                →  05_polarize_inversion.py [--mtdna]
                →  OUT JSON  →  ngsTracts
```
