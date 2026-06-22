# Bloc 11 — Hemizygous DEL markers ("fake trio" direction inference)

| | |
|---|---|
| **module** | `src/hpp/hemizygous_check.py` |
| **schema** | `schemas/del_markers.in.schema.json` (`ngspedigree_del_markers_v1`) |
| **purpose** | Mendelian DEL-transmission scoring to discriminate which of three candidate samples is the offspring in a putative trio |

## Why a separate bloc from inversion polarization (bloc 06)

The math is the same biallelic Punnett enumeration. The biological
interpretation and the short-read caveats are different enough to keep
the modules apart:

| | inversion polarization (06) | hemizygous DEL (11) |
|---|---|---|
| marker | PC1-band karyotype call (HOM1/HET/HOM2) | structural deletion genotype (0/0, 0/1, 1/1) |
| polarity hint required | yes (external) | no — REF/DEL is intrinsic |
| `1/1` reliability | high (homozygous arrangements call cleanly) | **fragile** (short-read homozygous DEL can mimic mapping dropout) |
| SNPs inside the marker | (n/a — inversion is large) | **must be masked from SNP Mendelian** to avoid allele-dropout artefacts |

## The "fake trio" use case

θ relatedness identifies first-degree dyads but cannot orient triads:
`P1+P2 → C` is kinship-symmetric with `P1+C → P2` and `P2+C → P1`.
DEL Mendelian errors discriminate the three permutations — the true
child has the lowest DEL error rate against the inferred parents.

```python
from hpp.hemizygous_check import load_del_markers, best_direction

calls = load_del_markers("del_markers.json")
v = best_direction(
    a="S001", b="S002", c="S003",
    del_calls_by_marker=calls,
    informative_marker_floor=5,   # don't pick a winner from noise
    min_margin=0.05,              # require at least 5% lower error than 2nd-best
)
print(v.best_direction)   # e.g. "S001+S002->S003"
print(v.best_error_rate)  # 0.0
print(v.margin)           # gap to second-best direction
```

`best_direction` returns `None` for the winner when no direction
meets the informative-marker floor or the error-rate margin is too
small — refusing to call a winner from noise.

## Transmission table (Punnett)

| P1 × P2 | Allowed child |
|---|---|
| 0/0 × 0/0 | {0/0} |
| 0/0 × 0/1 | {0/0, 0/1} |
| 0/0 × 1/1 | {0/1} |
| 0/1 × 0/1 | {0/0, 0/1, 1/1} |
| 0/1 × 1/1 | {0/1, 1/1} |
| 1/1 × 1/1 | {1/1} |

Symmetric in (P1, P2).

## Two reliability modes

- **`strict_hom_del=True` (default)** — trust 1/1 calls; mark
  (0/0 × 1/1 → 0/0) and (parent=0/0, offspring=1/1) etc. as
  hard incompatibilities.
- **`strict_hom_del=False`** — treat 1/1 as ambiguous. Any marker
  with a 1/1 GT in any of (P1, P2, child) becomes uninformative
  rather than incompatible. Use this when short-read mapping is
  unreliable in DEL regions.

## What this bloc does NOT do

- It does **not** call DELs from reads — input is the per-(marker,
  sample) genotype call.
- It does **not** mask SNPs inside DEL intervals from SNP-Mendelian
  checks — that masking belongs to whatever feeds the SNP Mendelian
  step. The IN-JSON contract carries optional `chrom`/`start`/`end`
  per marker so consumers can do the masking themselves.
- It does **not** infer paternity/maternity directly. It scores
  triad-direction compatibility; sex assignment still comes from
  upstream (sex map, or mtDNA pre-flight from bloc 07).

## Input contract

`schemas/del_markers.in.schema.json` (`ngspedigree_del_markers_v1`):

```json
{
  "schema": "ngspedigree_del_markers_v1",
  "rows": [
    {"marker_id": "DEL_001", "sample_id": "S001", "genotype": "0/1",
     "chrom": "Chr1", "start": 100000, "end": 105000,
     "depth_ratio": 0.51, "confidence": "high"},
    ...
  ]
}
```

Required: `marker_id`, `sample_id`, `genotype`.
Optional: `chrom`, `start`, `end`, `depth_ratio`, `confidence`, `notes`.

`load_del_markers(path)` returns the nested map
`{marker_id: {sample_id: genotype_str}}` that the triad and dyad
scorers consume.

## How it fits into the backbone

```
ngsRelate .res        →  Stage 1     →  candidate dyads (θ + IBS0)
                           ↓
                       ambiguous direction
                           ↓
del_markers.in.json   →  bloc 11  ──→  best_direction() picks the trio
                                       offspring with lowest error rate
                           ↓
                       confirmed triad   →  pedigree_extract  →  polarization
```

Bloc 11 is optional — when DEL markers exist, it lets the backbone
break direction ties; when they don't, the pipeline falls back to
roster-based extraction (blind-mode parent_a/parent_b convention).
