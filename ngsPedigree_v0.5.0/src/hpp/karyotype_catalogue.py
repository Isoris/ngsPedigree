"""
Karyotype-catalogue JSON adapter.

Read shape (the data-registry contract, as imagined for the
not-yet-hired librarian): one whole-genome catalogue carrying every
karyotype call for every (sample, inversion-or-LRR) pair. The
catalogue is supplied as JSON wrapping a long row table.

Wire format (``ngspedigree_karyotype_catalogue_v1``):

```json
{
  "schema": "ngspedigree_karyotype_catalogue_v1",
  "rows": [
    {"chrom": "Chr1", "lrr_id": "LRR_001", "sample_id": "S001", "karyotype": "HOM1"},
    {"chrom": "Chr1", "lrr_id": "LRR_001", "sample_id": "S002", "karyotype": "HET"},
    {"chrom": "Chr1", "lrr_id": "LRR_002", "sample_id": "S001", "karyotype": "HOM2"},
    {"chrom": "Chr2", "lrr_id": "LRR_003", "sample_id": "S001", "karyotype": "HOM1"},
    ...
  ]
}
```

The ``karyotype`` label uses the registry's nomenclature:

  - ``HOM1`` ↔ band 0 (lower PC1 cluster)
  - ``HET``  ↔ band 1 (middle / heterozygous)
  - ``HOM2`` ↔ band 2 (upper PC1 cluster)

Optional per-row fields: ``confidence`` (high/medium/low — default
"high"), ``inversion_id`` (an explicit polarization-side identifier
when it differs from ``lrr_id``).

Mapping to the polarization IN JSON:

The catalogue is a whole-genome table; the polarization bloc consumes
per-inversion karyotype calls. ``filter_to_inversion`` extracts the
subset for one LRR / inversion ID and converts each row to a
``KaryotypeCall`` ready to be fed into ``build_polarization_input`` or
``polarize`` directly.

This contract is forward-compatible: it defines the registry's output
shape before the registry exists. When the librarian (registry
implementer) ships, the loader stays the same; only the upstream
producer changes.
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

from .inversion_polarization import KaryotypeCall

CATALOGUE_SCHEMA_VERSION = "ngspedigree_karyotype_catalogue_v1"

# Registry label → PC1 band integer.
KARYOTYPE_LABEL_TO_BAND = {
    "HOM1": 0,
    "HET":  1,
    "HOM2": 2,
}
BAND_TO_LABEL = {v: k for k, v in KARYOTYPE_LABEL_TO_BAND.items()}

ALLOWED_LABELS = set(KARYOTYPE_LABEL_TO_BAND)


class KaryotypeCatalogueError(ValueError):
    pass


@dataclass(frozen=True)
class CatalogueRow:
    chrom: str
    lrr_id: str
    sample_id: str
    karyotype: str            # HOM1 / HET / HOM2
    band: int                 # 0 / 1 / 2
    confidence: str = "high"
    inversion_id: Optional[str] = None  # falls back to lrr_id when None


@dataclass
class KaryotypeCatalogue:
    rows: List[CatalogueRow]

    # ------------------------------------------------------------------
    # Indexing.
    # ------------------------------------------------------------------

    def lrrs(self) -> List[str]:
        seen: List[str] = []
        seen_set: Set[str] = set()
        for r in self.rows:
            if r.lrr_id not in seen_set:
                seen.append(r.lrr_id)
                seen_set.add(r.lrr_id)
        return seen

    def chroms(self) -> List[str]:
        return sorted({r.chrom for r in self.rows})

    def samples(self) -> List[str]:
        seen: List[str] = []
        seen_set: Set[str] = set()
        for r in self.rows:
            if r.sample_id not in seen_set:
                seen.append(r.sample_id)
                seen_set.add(r.sample_id)
        return seen

    def rows_for_lrr(self, lrr_or_inversion: str) -> List[CatalogueRow]:
        out: List[CatalogueRow] = []
        for r in self.rows:
            if r.lrr_id == lrr_or_inversion:
                out.append(r)
                continue
            if r.inversion_id == lrr_or_inversion:
                out.append(r)
        return out

    def filter_to_inversion(
        self,
        lrr_or_inversion_id: str,
        *,
        sample_whitelist: Optional[Iterable[str]] = None,
    ) -> List[KaryotypeCall]:
        """Convert the rows for one LRR (or inversion_id) into a list of
        ``KaryotypeCall`` records suitable for ``polarize`` and the IN
        JSON. ``sample_whitelist`` (if given) restricts to those samples.
        """
        wl = set(sample_whitelist) if sample_whitelist is not None else None
        out: List[KaryotypeCall] = []
        seen_sample: Set[str] = set()
        for r in self.rows_for_lrr(lrr_or_inversion_id):
            if wl is not None and r.sample_id not in wl:
                continue
            if r.sample_id in seen_sample:
                raise KaryotypeCatalogueError(
                    f"duplicate karyotype call for sample {r.sample_id!r} "
                    f"at LRR {lrr_or_inversion_id!r}"
                )
            seen_sample.add(r.sample_id)
            out.append(KaryotypeCall(
                sample_id=r.sample_id,
                inversion_id=lrr_or_inversion_id,
                band=r.band,
                confidence=r.confidence,
            ))
        return out

    # ------------------------------------------------------------------
    # Stats / debug.
    # ------------------------------------------------------------------

    def coverage(self) -> Dict[str, int]:
        """{lrr_id: n_samples_with_calls}."""
        counter: Dict[str, int] = defaultdict(int)
        for r in self.rows:
            counter[r.lrr_id] += 1
        return dict(counter)


# ----------------------------------------------------------------------
# Loader.
# ----------------------------------------------------------------------


def _row_loc(idx: int, base: str = "") -> str:
    return f"{base}/rows/{idx}" if base else f"rows[{idx}]"


def load_catalogue(path: str | Path) -> KaryotypeCatalogue:
    """Load a `ngspedigree_karyotype_catalogue_v1` JSON file."""
    path = Path(path)
    with open(path) as fh:
        doc = json.load(fh)
    return parse_catalogue(doc, source=str(path))


def parse_catalogue(
    doc: dict, *, source: str = "<dict>"
) -> KaryotypeCatalogue:
    schema = doc.get("schema")
    if schema != CATALOGUE_SCHEMA_VERSION:
        raise KaryotypeCatalogueError(
            f"{source}: schema must be {CATALOGUE_SCHEMA_VERSION!r}; got {schema!r}"
        )
    rows_raw = doc.get("rows")
    if rows_raw is None:
        raise KaryotypeCatalogueError(f"{source}: missing 'rows' field")
    rows: List[CatalogueRow] = []
    for i, r in enumerate(rows_raw):
        loc = _row_loc(i, source)
        for req in ("chrom", "lrr_id", "sample_id", "karyotype"):
            if req not in r:
                raise KaryotypeCatalogueError(
                    f"{loc}: missing required field {req!r}"
                )
        label = r["karyotype"]
        if label not in ALLOWED_LABELS:
            raise KaryotypeCatalogueError(
                f"{loc}: karyotype must be in {sorted(ALLOWED_LABELS)}; "
                f"got {label!r}"
            )
        rows.append(CatalogueRow(
            chrom=str(r["chrom"]),
            lrr_id=str(r["lrr_id"]),
            sample_id=str(r["sample_id"]),
            karyotype=label,
            band=KARYOTYPE_LABEL_TO_BAND[label],
            confidence=str(r.get("confidence") or "high"),
            inversion_id=(str(r["inversion_id"])
                          if r.get("inversion_id") else None),
        ))
    return KaryotypeCatalogue(rows=rows)


def write_catalogue(
    path: str | Path, catalogue: KaryotypeCatalogue, *, indent: int = 2,
) -> Path:
    """Round-trip writer — useful for synthesizing fixtures."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "schema": CATALOGUE_SCHEMA_VERSION,
        "rows": [
            {
                "chrom": r.chrom,
                "lrr_id": r.lrr_id,
                "sample_id": r.sample_id,
                "karyotype": r.karyotype,
                **({"confidence": r.confidence} if r.confidence != "high" else {}),
                **({"inversion_id": r.inversion_id} if r.inversion_id else {}),
            }
            for r in catalogue.rows
        ],
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=indent)
    return path


# ----------------------------------------------------------------------
# Conversion to the polarization IN JSON shape.
# ----------------------------------------------------------------------


def catalogue_calls_to_in_json_array(
    catalogue: KaryotypeCatalogue, lrr_or_inversion_id: str,
) -> List[dict]:
    """Flatten the catalogue's per-LRR slice into the
    ``karyotype_calls`` array shape used by the polarization IN JSON
    (``ngspedigree_karyotype_calls_in_v1``).
    """
    return [
        {
            "sample_id": c.sample_id,
            "band": c.band,
            "confidence": c.confidence,
        }
        for c in catalogue.filter_to_inversion(lrr_or_inversion_id)
    ]
