# ngsPedigree output schema

Three artifacts per run, written to `--outdir`.

---

## `pairwise_relationship_classification.tsv`

Edge-level annotations. One row per pair from the input `.res`.

| Column | Type | Description |
|--------|------|-------------|
| `sample_a` | str | First sample ID (resolved from `a` index via `samples.txt`) |
| `sample_b` | str | Second sample ID |
| `a` | int | Original ngsRelate `a` index (0-based) |
| `b` | int | Original ngsRelate `b` index |
| `nSites` | int | Number of sites compared (from ngsRelate, if present) |
| `theta` | float | Kinship coefficient (from ngsRelate) |
| `IBS0` | float | Proportion of sites IBS=0 (from ngsRelate) |
| `KING` | float | KING-robust kinship estimator (if present) |
| `R0` | float | IBS-derived ratio (if present) |
| `R1` | float | IBS-derived ratio (if present) |
| `J9` | float | Jacquard 9 component (if present) |
| `edge_class` | str | One of: `duplicate_or_clone`, `parent_offspring`, `full_sibling`, `ambiguous_first_degree`, `second_degree`, `third_degree`, `unrelated`, `undetermined` |
| `confidence` | str | `high`, `medium`, or `low` |
| `reasons` | str | Semicolon-separated diagnostic notes (e.g. `low_n_sites_500`, `king_low_for_first_degree_0.140`) |

### Notes on use

- The table is **complete** — every pair from `.res` has a row, including
  unrelated pairs. Filter on `edge_class` to get a subgraph.
- For Stage 2 (per-chromosome): additional columns will be appended in-place,
  named `edge_class_LG01`, `edge_class_LG02`, etc.
- For reverse queries ("candidate parents of CGA187"), filter:
  `(sample_a == 'CGA187' OR sample_b == 'CGA187') AND edge_class == 'parent_offspring'`.

---

## `family_hub_roster.tsv`

Node-level annotations. One row per sample.

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | str | Sample ID |
| `hub_id` | str | Hub identifier (`H001`, `H002`, ... or `_isolated`) |
| `hub_type` | str | One of: `two_parents_with_sibship`, `parent_with_sibship`, `sibship_only`, `po_dyad_only`, `duplicate_pair`, `mixed_or_complex`, `singleton`, `isolated` |
| `hub_size` | int | Number of samples in the hub |
| `possible_role` | str | See "Roles" below |
| `role_confidence` | str | `high`, `medium`, or `low` |
| `reason` | str | Free-text reason for the role assignment |

### Roles

Blind mode (no `--sex`):

| Role | Meaning |
|------|---------|
| `forced_parent` | Has PO edges to ≥3 mutually-FS nodes; must be the parent |
| `likely_parent` | Has PO edges to exactly 2 mutually-FS nodes |
| `parent_a`, `parent_b` | One of two parents in a `two_parents_with_sibship` hub |
| `possible_offspring` | PO-related to one or two inferred parents |
| `possible_full_sib` | In a `sibship_only` hub |
| `ambiguous_first_degree_PO` | Single PO edge, can't tell which is parent |
| `ambiguous_in_parent_hub` | In a parent_with_sibship hub but not PO to the inferred parent |
| `ambiguous_in_two_parent_hub` | In a two-parent hub but not PO to both parents |
| `ambiguous_first_degree` | In a mixed/complex hub, no clean role |
| `ambiguous` | Catch-all for unexpected configurations |
| `duplicate` | In a `duplicate_pair` hub |
| `isolated` | No first-degree edges |

With `--sex` provided:

| Promoted role | Promotion rule |
|---------------|----------------|
| `mother` | Was `parent_a`/`parent_b`/`forced_parent`/`likely_parent` AND `sex=female` |
| `father` | Same role conditions AND `sex=male` |

If both parents in a two-parent hub are same-sex (impossible in *C. gariepinus*),
the roles are NOT promoted; instead, a warning is appended to the `reason`
column (`same_sex_male_parent_pair_unusual` or `same_sex_female_parent_pair_unusual`).

---

## `ngspedigree_run_envelope.json`

Self-describing wrapper. Schema name: `ngspedigree_run_envelope_v1`.

```json
{
  "schema": "ngspedigree_run_envelope_v1",
  "produced_by": {
    "tool": "STEP_PED_01_annotate_relationships",
    "version": "v0.1.0",
    "params": {
      "thresholds": { "theta_first": 0.177, "theta_second": 0.0884,
                       "theta_third": 0.0442, "theta_dup_min": 0.45,
                       "ibs0_po_max": 0.005 },
      "sex_provided": false,
      "n_sex_samples_provided": 0
    }
  },
  "inputs": {
    "res_path": "...",
    "samples_path": "...",
    "sex_path": null
  },
  "run_id": "cohort_226_full_v1",
  "mode": "blind",
  "computed_at": "2026-05-10T...",
  "n_samples": 226,
  "n_pairs": 25425,
  "n_pairs_by_class": {
    "unrelated": 25000,
    "parent_offspring": 80,
    "full_sibling": 250,
    "duplicate_or_clone": 2,
    "second_degree": 50,
    "third_degree": 30,
    "ambiguous_first_degree": 13
  },
  "n_hubs_by_type": {
    "two_parents_with_sibship": 5,
    "parent_with_sibship": 12,
    "sibship_only": 20,
    "po_dyad_only": 8,
    "duplicate_pair": 1,
    "mixed_or_complex": 1
  },
  "n_multi_node_components": 47,
  "n_singletons": 21,
  "n_label_promotions_via_sex": 0,
  "artifacts": {
    "pairwise_relationship_classification": "pairwise_relationship_classification.tsv",
    "family_hub_roster": "family_hub_roster.tsv"
  },
  "schema_status": "v1_draft",
  "downstream": {
    "next_step": "STEP_PED_02_per_chromosome_annotation.py",
    "consumers": [
      "family_segregation_gate (Slice 3)",
      "STEP_PED_03_inheritance_map.py",
      "ngsTracts (NCO/DCO calling)"
    ]
  }
}
```

## Versioning

Schema version is in the `schema` field (currently `ngspedigree_run_envelope_v1`).
When the schema changes incompatibly, the version bumps and old consumers
should refuse to read the new shape.
