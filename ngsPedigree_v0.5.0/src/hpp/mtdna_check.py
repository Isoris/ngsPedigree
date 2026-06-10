"""
mtDNA maternal-lineage validation.

Mitochondrial haplotypes are inherited maternally with no paternal
contribution and no recombination. The implication for pedigree work
is asymmetric:

  - Different mtDNA haplotype between a candidate mother and her
    proposed offspring is **strong rejection** of maternity. Common
    causes: wrong-mother label, sample swap, or tank-mixing in
    hatchery data.
  - The same mtDNA haplotype is **compatible** with the proposed
    mother, but is not unique proof — multiple females in the cohort
    can share a haplotype.
  - mtDNA does **not** inform paternity (males do not transmit mtDNA).
  - mtDNA does **not** inform nuclear-inversion polarization
    (REF vs INV is a nuclear question; mtDNA carries no signal about it).

This module therefore implements layer 2 of the inversion-inheritance
stack — pedigree validation — and runs as a pre-flight check before
the polarization bloc. Maternal-(dyad/triad) relationships flagged
``incompatible`` here are excluded from the polarization input set so
they do not contaminate the Mendelian compatibility counts or the
transmission calls. Paternal dyads are never filtered by mtDNA — it is
not informative for them.

Input contract: a TSV with required columns ``sample_id``,
``mtdna_haplotype``, and optional ``mtdna_sequence``, ``mtdna_n_sites``.
When sequences of equal length are present on both members of a pair,
a Hamming distance is computed and a small distance (≤ threshold,
default 2) reverses the haplotype-label mismatch to ``compatible`` —
this is for the case where two clusters carry near-identical sequence
and the label disagreement is essentially a clustering boundary.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .inversion_polarization import DyadPair, TriadTrio


HAMMING_THRESHOLD_DEFAULT = 2

REQUIRED_COLUMNS = ("sample_id", "mtdna_haplotype")
OPTIONAL_COLUMNS = ("mtdna_sequence", "mtdna_n_sites")


class MtdnaContractError(ValueError):
    pass


@dataclass(frozen=True)
class MtdnaRecord:
    sample_id: str
    haplotype: str
    sequence: Optional[str] = None
    n_sites: Optional[int] = None


@dataclass(frozen=True)
class MtdnaCheck:
    mother_sample_id: str
    offspring_sample_id: str
    relationship_type: str        # "triad" | "maternal_dyad"
    mother_haplotype: Optional[str]
    offspring_haplotype: Optional[str]
    distance: Optional[int]
    status: str                   # "compatible" | "incompatible" | "ambiguous"
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# Loader.
# ----------------------------------------------------------------------


def load_mtdna_haplotypes(path: str | Path) -> Dict[str, MtdnaRecord]:
    """Load the mtDNA TSV and return ``{sample_id: MtdnaRecord}``."""
    path = Path(path)
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        if not header_line:
            raise MtdnaContractError(f"{path}: empty file")
        header = header_line.split("\t")
        for c in REQUIRED_COLUMNS:
            if c not in header:
                raise MtdnaContractError(f"{path}: missing required column {c!r}")

        out: Dict[str, MtdnaRecord] = {}
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                raise MtdnaContractError(
                    f"{path}:{n}: expected {len(header)} fields, got {len(parts)}"
                )
            row = dict(zip(header, parts))
            sid = row["sample_id"]
            hap = row["mtdna_haplotype"]
            if not sid or not hap:
                raise MtdnaContractError(
                    f"{path}:{n}: sample_id and mtdna_haplotype must be non-empty"
                )
            seq = row.get("mtdna_sequence") or None
            n_sites_raw = row.get("mtdna_n_sites")
            n_sites = int(n_sites_raw) if n_sites_raw else None
            if sid in out and (out[sid].haplotype != hap or out[sid].sequence != seq):
                raise MtdnaContractError(
                    f"{path}:{n}: contradictory mtDNA record for {sid!r}"
                )
            out[sid] = MtdnaRecord(sample_id=sid, haplotype=hap,
                                   sequence=seq, n_sites=n_sites)
    return out


# ----------------------------------------------------------------------
# Pair check.
# ----------------------------------------------------------------------


def _hamming(a: str, b: str) -> Optional[int]:
    if len(a) != len(b):
        return None
    return sum(1 for x, y in zip(a, b) if x != y)


def check_pair(
    mother_sample_id: str,
    offspring_sample_id: str,
    mt_records: Dict[str, MtdnaRecord],
    *,
    relationship_type: str = "triad",
    hamming_threshold: int = HAMMING_THRESHOLD_DEFAULT,
) -> MtdnaCheck:
    """Compare mother and offspring mtDNA records.

    Decision order:
      1. either record missing                   -> ambiguous
      2. haplotype labels equal                  -> compatible (distance 0)
      3. labels differ, sequences both present,
         equal length, Hamming <= threshold      -> compatible (close)
      4. labels differ, sequences both present,
         Hamming > threshold                     -> incompatible
      5. labels differ, no usable sequence pair  -> incompatible
    """
    m = mt_records.get(mother_sample_id)
    o = mt_records.get(offspring_sample_id)
    if m is None or o is None:
        missing = [
            sid for sid, rec in
            [(mother_sample_id, m), (offspring_sample_id, o)] if rec is None
        ]
        return MtdnaCheck(
            mother_sample_id=mother_sample_id,
            offspring_sample_id=offspring_sample_id,
            relationship_type=relationship_type,
            mother_haplotype=(m.haplotype if m else None),
            offspring_haplotype=(o.haplotype if o else None),
            distance=None,
            status="ambiguous",
            note=f"no mtDNA record for {','.join(missing)}",
        )
    if m.haplotype == o.haplotype:
        return MtdnaCheck(
            mother_sample_id=mother_sample_id,
            offspring_sample_id=offspring_sample_id,
            relationship_type=relationship_type,
            mother_haplotype=m.haplotype,
            offspring_haplotype=o.haplotype,
            distance=0,
            status="compatible",
            note="haplotype label match",
        )
    # Labels differ.
    if m.sequence and o.sequence:
        d = _hamming(m.sequence, o.sequence)
        if d is not None and d <= hamming_threshold:
            return MtdnaCheck(
                mother_sample_id=mother_sample_id,
                offspring_sample_id=offspring_sample_id,
                relationship_type=relationship_type,
                mother_haplotype=m.haplotype,
                offspring_haplotype=o.haplotype,
                distance=d,
                status="compatible",
                note=(f"label differs but Hamming distance {d} ≤ "
                      f"{hamming_threshold}"),
            )
        return MtdnaCheck(
            mother_sample_id=mother_sample_id,
            offspring_sample_id=offspring_sample_id,
            relationship_type=relationship_type,
            mother_haplotype=m.haplotype,
            offspring_haplotype=o.haplotype,
            distance=d,
            status="incompatible",
            note=(f"haplotype labels differ; Hamming={d}, "
                  f"threshold={hamming_threshold}"),
        )
    return MtdnaCheck(
        mother_sample_id=mother_sample_id,
        offspring_sample_id=offspring_sample_id,
        relationship_type=relationship_type,
        mother_haplotype=m.haplotype,
        offspring_haplotype=o.haplotype,
        distance=None,
        status="incompatible",
        note="haplotype labels differ; no sequence available for distance",
    )


# ----------------------------------------------------------------------
# Pedigree-wide check.
# ----------------------------------------------------------------------


def check_pedigree(
    mt_records: Dict[str, MtdnaRecord],
    dyads: Sequence[DyadPair],
    triads: Sequence[TriadTrio],
    *,
    hamming_threshold: int = HAMMING_THRESHOLD_DEFAULT,
) -> List[MtdnaCheck]:
    """Run the mtDNA check on every (mother, offspring) pair appearing in
    the triad list or as a maternal dyad. (Triad maternal-side wins
    over an overlapping maternal-dyad row — they refer to the same
    pair.) Paternal dyads are never checked because mtDNA carries no
    paternity signal."""
    out: List[MtdnaCheck] = []
    seen: set = set()
    for tr in triads:
        key = (tr.maternal_sample_id, tr.offspring_sample_id)
        if key in seen:
            continue
        out.append(check_pair(
            tr.maternal_sample_id, tr.offspring_sample_id, mt_records,
            relationship_type="triad",
            hamming_threshold=hamming_threshold,
        ))
        seen.add(key)
    for dy in dyads:
        if dy.parent_sex != "female":
            continue
        key = (dy.parent_sample_id, dy.offspring_sample_id)
        if key in seen:
            continue
        out.append(check_pair(
            dy.parent_sample_id, dy.offspring_sample_id, mt_records,
            relationship_type="maternal_dyad",
            hamming_threshold=hamming_threshold,
        ))
        seen.add(key)
    return out


# ----------------------------------------------------------------------
# Filter step.
# ----------------------------------------------------------------------


def filter_by_mtdna(
    dyads: Sequence[DyadPair],
    triads: Sequence[TriadTrio],
    checks: Sequence[MtdnaCheck],
) -> Tuple[List[DyadPair], List[TriadTrio]]:
    """Drop triads + maternal dyads where the mtDNA check returned
    ``incompatible``. Compatible / ambiguous pairs are kept (ambiguous
    cases are reported but never filter, since absence of mtDNA is not
    evidence of pedigree error)."""
    incompatible = {
        (c.mother_sample_id, c.offspring_sample_id)
        for c in checks
        if c.status == "incompatible"
    }
    filt_triads = [
        t for t in triads
        if (t.maternal_sample_id, t.offspring_sample_id) not in incompatible
    ]
    filt_dyads = [
        d for d in dyads
        if not (
            d.parent_sex == "female"
            and (d.parent_sample_id, d.offspring_sample_id) in incompatible
        )
    ]
    return filt_dyads, filt_triads


# ----------------------------------------------------------------------
# Summary block for the OUT JSON.
# ----------------------------------------------------------------------


def build_validation_block(
    checks: Sequence[MtdnaCheck],
    *,
    n_triads_excluded: int,
    n_dyads_excluded: int,
    hamming_threshold: int = HAMMING_THRESHOLD_DEFAULT,
) -> dict:
    return {
        "supplied": True,
        "n_relationships_checked": len(checks),
        "n_compatible": sum(1 for c in checks if c.status == "compatible"),
        "n_incompatible": sum(1 for c in checks if c.status == "incompatible"),
        "n_ambiguous": sum(1 for c in checks if c.status == "ambiguous"),
        "n_triads_excluded": n_triads_excluded,
        "n_dyads_excluded": n_dyads_excluded,
        "hamming_threshold": hamming_threshold,
        "checks": [c.to_dict() for c in checks],
    }


def build_not_supplied_block() -> dict:
    return {"supplied": False}
