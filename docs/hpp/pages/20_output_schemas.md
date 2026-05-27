# Output schemas — Tables A through E

All HPP outputs are TSV with a sidecar JSON Schema. The TSV
extension follows the project convention; the JSON Schema is the
machine-readable contract.

| Table | Layer ID | Schema file | Bloc that produces it |
|---|---|---|---|
| A | `hpp_offspring_haplotype_variants` | `schemas/A_hpp_offspring_haplotype_variants.schema.json` | `haplotype_projection` |
| B | `hpp_offspring_gene_status` | `schemas/B_hpp_offspring_gene_status.schema.json` | `offspring_gene_status` |
| C | `hpp_dyad_transmission_summary` (part of `hpp_transmission_summary`) | `schemas/C_hpp_dyad_transmission_summary.schema.json` | `transmission_summary` (dyad) |
| D | `hpp_triad_transmission_summary` (part of `hpp_transmission_summary`) | `schemas/D_hpp_triad_transmission_summary.schema.json` | `transmission_summary` (triad) |
| E | `hpp_kbc_crosscheck` | `schemas/E_hpp_kbc_crosscheck.schema.json` | `hpp_kbc_arrangement_crosscheck` |

Counts:

| Table | Columns | Rows-per |
|---|---:|---|
| A | 20 | (offspring × variant × hap_copy) |
| B | 13 | (offspring × gene × segment) |
| C | 19 | (dyad) |
| D | 23 | (triad) |
| E | 10 | (offspring × variant × inversion) |

## Common cross-table fields

- `damaging_tier` — `T1` / `T2` / `T3` (KBC §1.8). Carried on every row.
- `confidence` — `high` / `medium` / `low` / `unresolved`. Composite per SPEC §5.
- `segment_confidence` — `Gold` / `Silver` / `Bronze`. From Stage 3.

## Loading schemas at runtime

```python
import json, pathlib
SCHEMA_DIR = pathlib.Path("ngsPedigree_v0.5.0/schemas")
table_A_schema = json.loads((SCHEMA_DIR / "A_hpp_offspring_haplotype_variants.schema.json").read_text())
declared_columns = [c["enum"][0] if "enum" in c else c["const"]
                    for c in table_A_schema["properties"]["columns"]["items"]["enum"]]
```

(Or — simpler — use the `enum` block directly:
`table_A_schema["properties"]["columns"]["items"]["enum"]`.)

## Hand-off into atlas-core's `layer_registry`

Each of the five output layers also has a row in
`docs/hpp/atlas_core_registry/layer_registry.jsonl`. Atlas-core's
smoke test enforces that every `analysis_modes.produces` value
resolves to a row in `layer_registry`.
