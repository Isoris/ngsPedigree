"""
ngsTracts JSON adapters — IN and OUT.

IN  (consumed by ngsPedigree):
    karyotype_calls.in.json  — per-inversion karyotype band table
                                (band ∈ {0, 1, 2}) plus a polarity hint
                                supplied by the caller (typically derived
                                from the reference assembly's arrangement
                                at the inversion).

OUT (consumed by ngsTracts):
    polarized_transmissions.out.json — chosen polarity, per-(parent,
                                offspring) transmitted-arrangement calls,
                                and the Mendelian drive test. ngsTracts
                                uses transmitted_arrangement as the phase
                                polarity for marker-level CO/NCO/DCO
                                scanning over the inversion interval.

Both files validate against the JSON Schemas in ../schemas/:
    karyotype_calls.in.schema.json
    polarized_transmissions.out.schema.json

The schemas themselves are conventional ("this is the field set") rather
than enforced at load time — the IN loader validates required fields and
enum values and raises a typed exception on contract violation.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .inversion_polarization import (
    ALLOWED_BANDS,
    ALLOWED_POLARITY_HINTS,
    DriveStats,
    DyadPair,
    KaryotypeCall,
    PolarizationResult,
    TransmissionCall,
    TriadTrio,
)

NGSTRACTS_IN_SCHEMA_VERSION = "ngspedigree_karyotype_calls_in_v1"
NGSTRACTS_OUT_SCHEMA_VERSION = "ngspedigree_polarized_transmissions_v1"


class NgsTractsAdapterError(ValueError):
    pass


# ----------------------------------------------------------------------
# IN adapter.
# ----------------------------------------------------------------------


def _require(d: dict, key: str, where: str):
    if key not in d:
        raise NgsTractsAdapterError(f"{where}: missing required field {key!r}")
    return d[key]


def load_karyotype_calls(path: str | Path) -> Tuple[
    str, str, List[KaryotypeCall], List[DyadPair], List[TriadTrio]
]:
    """Load an IN-side JSON file.

    Returns:
        (inversion_id, polarity_hint, karyotype_calls, dyads, triads)
    """
    path = Path(path)
    with open(path) as fh:
        doc = json.load(fh)

    schema = _require(doc, "schema", f"{path}")
    if schema != NGSTRACTS_IN_SCHEMA_VERSION:
        raise NgsTractsAdapterError(
            f"{path}: schema must be {NGSTRACTS_IN_SCHEMA_VERSION!r}; "
            f"got {schema!r}"
        )

    inversion_id = _require(doc, "inversion_id", f"{path}")
    polarity_hint = _require(doc, "polarity_hint", f"{path}")
    if polarity_hint not in ALLOWED_POLARITY_HINTS:
        raise NgsTractsAdapterError(
            f"{path}: polarity_hint must be in {ALLOWED_POLARITY_HINTS}; "
            f"got {polarity_hint!r}"
        )

    calls: List[KaryotypeCall] = []
    for i, c in enumerate(_require(doc, "karyotype_calls", f"{path}")):
        loc = f"{path}#/karyotype_calls/{i}"
        band = int(_require(c, "band", loc))
        if band not in ALLOWED_BANDS:
            raise NgsTractsAdapterError(
                f"{loc}: band must be in {ALLOWED_BANDS}; got {band!r}"
            )
        calls.append(KaryotypeCall(
            sample_id=str(_require(c, "sample_id", loc)),
            inversion_id=inversion_id,
            band=band,
            confidence=str(c.get("confidence", "high")),
        ))

    dyads: List[DyadPair] = []
    for i, d in enumerate(doc.get("dyads", [])):
        loc = f"{path}#/dyads/{i}"
        dyads.append(DyadPair(
            parent_sample_id=str(_require(d, "parent_sample_id", loc)),
            offspring_sample_id=str(_require(d, "offspring_sample_id", loc)),
            parent_sex=(str(d["parent_sex"]) if d.get("parent_sex") else None),
        ))

    triads: List[TriadTrio] = []
    for i, t in enumerate(doc.get("triads", [])):
        loc = f"{path}#/triads/{i}"
        triads.append(TriadTrio(
            paternal_sample_id=str(_require(t, "paternal_sample_id", loc)),
            maternal_sample_id=str(_require(t, "maternal_sample_id", loc)),
            offspring_sample_id=str(_require(t, "offspring_sample_id", loc)),
        ))

    return inversion_id, polarity_hint, calls, dyads, triads


# ----------------------------------------------------------------------
# OUT adapter.
# ----------------------------------------------------------------------


def write_polarized_transmissions(
    path: str | Path,
    *,
    result: PolarizationResult,
    transmissions: Sequence[TransmissionCall],
    drive_stats: DriveStats,
    intended_consumer: str = "ngsTracts",
    mtdna_validation: Optional[Dict] = None,
    extra_metadata: Optional[Dict] = None,
) -> Path:
    """Emit the OUT JSON. ngsTracts consumes this and starts its
    marker-level scan using transmitted_arrangement as the phase polarity
    per (parent, offspring, inversion).

    ``mtdna_validation``, when supplied, records the pre-flight
    maternal-lineage check that filtered the dyad/triad input set
    before polarization. Consumers can use it to distinguish nuclear
    Mendelian failures (visible in ``polarization.incompatible_*``)
    from maternal-pedigree failures (visible in
    ``mtdna_validation.checks``).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    doc = {
        "schema": NGSTRACTS_OUT_SCHEMA_VERSION,
        "intended_consumer": intended_consumer,
        "inversion_id": result.inversion_id,
        "polarization": result.to_dict(),
        "transmissions": [t.to_dict() for t in transmissions],
        "drive_stats": drive_stats.to_dict(),
    }
    if mtdna_validation is not None:
        doc["mtdna_validation"] = mtdna_validation
    if extra_metadata:
        doc["metadata"] = extra_metadata

    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)
    return path
