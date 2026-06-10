# Bloc 10 — Karyotype-catalogue JSON adapter (registry IN)

| | |
|---|---|
| **module** | `src/hpp/karyotype_catalogue.py` |
| **schema** | `schemas/karyotype_catalogue.in.schema.json` (`ngspedigree_karyotype_catalogue_v1`) |
| **purpose** | forward-compatible IN-side contract for whole-genome karyotype calls supplied by the (not-yet-implemented) data registry |

## Why this exists

The data registry will eventually emit one whole-genome catalogue:
every karyotype call for every (sample, LRR-or-inversion) pair, in a
single long-table JSON. We define that contract now so that:

1. The polarization pipeline (bloc 06) has a clean upstream IN shape.
2. When the registry ships, we change one loader call — the rest of
   the pipeline stays the same.

## Wire format

```json
{
  "schema": "ngspedigree_karyotype_catalogue_v1",
  "rows": [
    {"chrom": "Chr1", "lrr_id": "LRR_001", "sample_id": "S001", "karyotype": "HOM1"},
    {"chrom": "Chr1", "lrr_id": "LRR_001", "sample_id": "S002", "karyotype": "HET"},
    {"chrom": "Chr2", "lrr_id": "LRR_003", "sample_id": "S001", "karyotype": "HOM2"}
  ]
}
```

Required per row: `chrom`, `lrr_id`, `sample_id`, `karyotype`.

Optional per row:
- `confidence` — `high` (default), `medium`, `low`.
- `inversion_id` — overrides `lrr_id` when the polarization side uses
  a different identifier.

## Label → band mapping

| Catalogue label | PC1 band | Polarization arrangement (under `band_0_is_REF`) |
|---|---|---|
| `HOM1` | 0 | HOM_REF |
| `HET`  | 1 | HET |
| `HOM2` | 2 | HOM_INV |

## API

```python
from hpp.karyotype_catalogue import load_catalogue

cat = load_catalogue("registry_output.json")

# inventory
cat.lrrs()                # ['LRR_001', 'LRR_002', ...]
cat.samples()             # ['S001', 'S002', ...]
cat.coverage()            # {'LRR_001': n, ...}

# per-LRR slice → KaryotypeCall list, ready for polarize()
calls = cat.filter_to_inversion("LRR_001",
                                  sample_whitelist=["S001", "S002"])

# polarization-IN-JSON-shaped array for the karyotype_calls field
from hpp.karyotype_catalogue import catalogue_calls_to_in_json_array
arr = catalogue_calls_to_in_json_array(cat, "LRR_001")
```

## Validation

The loader is strict:
- bad / missing `schema` → `KaryotypeCatalogueError`
- missing required fields → `KaryotypeCatalogueError`
- `karyotype` not in `{HOM1, HET, HOM2}` → `KaryotypeCatalogueError`
- duplicate `sample_id` at one LRR → `KaryotypeCatalogueError`

## Round-trip

`write_catalogue(path, catalogue)` round-trips. Useful for synthesizing
test fixtures or replaying registry data.

## Status

Spec + loader + tests are built. The hand-off in the pipeline is:

```
data registry (not yet built)
        ↓  ngspedigree_karyotype_catalogue_v1 JSON
hpp.karyotype_catalogue.load_catalogue()
        ↓  KaryotypeCatalogue
.filter_to_inversion(LRR_or_inversion_id)
        ↓  List[KaryotypeCall]
hpp.inversion_polarization.polarize()  (bloc 06)
```

When the librarian (registry implementer) ships, only the producer
changes. The pipeline downstream is untouched.
