"""
Cross-caller SV catalogue merging (Delly + Manta).

A simple reciprocal-overlap matcher that unifies DEL calls from two
(or more) callers into one marker set per cohort. The merger never
*calls* a deletion — it only decides which Delly and Manta records
refer to the same underlying variant.

Default rule: two records merge if
  - same chromosome
  - left-breakpoint distance ≤ bp_tolerance (default 500 bp)
  - right-breakpoint distance ≤ bp_tolerance
  - reciprocal overlap ≥ reciprocal_overlap (default 0.5)

Genotype reconciliation per merged marker:
  - if both callers agree → keep the consensus
  - if they disagree (e.g. 0/0 vs 0/1) → emit "./." (uncertain) by
    default, or take Delly's call (or Manta's) under preferred_caller
  - missing on one side → use the other side
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .vcf_sv import SvRecord


@dataclass(frozen=True)
class MergedMarker:
    marker_id: str                # cohort-wide id, e.g. "DEL_Chr1_100000_105000"
    chrom: str
    start: int                    # consensus left breakpoint
    end: int                      # consensus right breakpoint
    sv_length: int
    callers: Tuple[str, ...]      # ("delly",) or ("manta",) or ("delly","manta")
    n_callers: int
    qual_max: Optional[float]
    contributing_record_ids: Tuple[str, ...]
    genotypes: Dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["callers"] = list(self.callers)
        d["contributing_record_ids"] = list(self.contributing_record_ids)
        return d


# ----------------------------------------------------------------------
# Matching predicate.
# ----------------------------------------------------------------------


def _reciprocal_overlap(a: SvRecord, b: SvRecord) -> float:
    if a.chrom != b.chrom:
        return 0.0
    inter = max(0, min(a.end, b.end) - max(a.pos, b.pos))
    if inter == 0:
        return 0.0
    len_a = max(1, a.end - a.pos)
    len_b = max(1, b.end - b.pos)
    return inter / max(len_a, len_b)


def _bp_distance_ok(a: SvRecord, b: SvRecord, tol: int) -> bool:
    return (abs(a.pos - b.pos) <= tol) and (abs(a.end - b.end) <= tol)


def records_match(a: SvRecord, b: SvRecord,
                  *, bp_tolerance: int = 500,
                  reciprocal_overlap: float = 0.5) -> bool:
    if a.chrom != b.chrom or a.sv_type != b.sv_type:
        return False
    if not _bp_distance_ok(a, b, bp_tolerance):
        return False
    return _reciprocal_overlap(a, b) >= reciprocal_overlap


# ----------------------------------------------------------------------
# Per-marker genotype reconciliation.
# ----------------------------------------------------------------------


def _reconcile_gt(
    delly_gt: str, manta_gt: str, *,
    preferred_caller: Optional[str] = None,
) -> str:
    if delly_gt == "./." and manta_gt == "./.":
        return "./."
    if delly_gt == "./.":
        return manta_gt
    if manta_gt == "./.":
        return delly_gt
    if delly_gt == manta_gt:
        return delly_gt
    if preferred_caller == "delly":
        return delly_gt
    if preferred_caller == "manta":
        return manta_gt
    return "./."


# ----------------------------------------------------------------------
# Merge across two record lists.
# ----------------------------------------------------------------------


def merge_two_callers(
    delly_records: Sequence[SvRecord],
    manta_records: Sequence[SvRecord],
    *,
    bp_tolerance: int = 500,
    reciprocal_overlap: float = 0.5,
    preferred_caller: Optional[str] = None,
) -> List[MergedMarker]:
    """Merge a Delly DEL record list with a Manta DEL record list.

    The Delly side dominates ordering: each Delly record either finds a
    Manta partner or stays solo. Manta records not claimed by a Delly
    match are emitted as solo entries with `callers=("manta",)`.

    The merger is deterministic — chrom-sorted iteration, first-match.
    """
    delly_sorted = sorted(delly_records, key=lambda r: (r.chrom, r.pos, r.end))
    manta_sorted = sorted(manta_records, key=lambda r: (r.chrom, r.pos, r.end))
    manta_claimed: Set[int] = set()
    merged: List[MergedMarker] = []

    for d in delly_sorted:
        partner_idx: Optional[int] = None
        for i, m in enumerate(manta_sorted):
            if i in manta_claimed:
                continue
            if m.chrom != d.chrom:
                continue
            if m.pos > d.end + bp_tolerance:
                break    # sorted; no further Manta record on this chrom can match
            if records_match(d, m,
                             bp_tolerance=bp_tolerance,
                             reciprocal_overlap=reciprocal_overlap):
                partner_idx = i
                break
        if partner_idx is None:
            merged.append(_record_to_marker(d))
            continue
        m = manta_sorted[partner_idx]
        manta_claimed.add(partner_idx)
        merged.append(_combine_records(d, m, preferred_caller=preferred_caller))

    for i, m in enumerate(manta_sorted):
        if i in manta_claimed:
            continue
        merged.append(_record_to_marker(m))

    merged.sort(key=lambda x: (x.chrom, x.start, x.end))
    return merged


def _record_to_marker(r: SvRecord) -> MergedMarker:
    return MergedMarker(
        marker_id=f"DEL_{r.chrom}_{r.pos}_{r.end}",
        chrom=r.chrom,
        start=r.pos,
        end=r.end,
        sv_length=r.svlen or (r.end - r.pos),
        callers=(r.caller,),
        n_callers=1,
        qual_max=r.qual,
        contributing_record_ids=(r.raw_id,),
        genotypes=dict(r.genotypes),
    )


def _combine_records(
    d: SvRecord, m: SvRecord, *, preferred_caller: Optional[str],
) -> MergedMarker:
    start = (d.pos + m.pos) // 2
    end = (d.end + m.end) // 2
    quals = [q for q in (d.qual, m.qual) if q is not None]
    qual = max(quals) if quals else None
    samples = set(d.genotypes) | set(m.genotypes)
    gts = {
        s: _reconcile_gt(
            d.genotypes.get(s, "./."),
            m.genotypes.get(s, "./."),
            preferred_caller=preferred_caller,
        )
        for s in samples
    }
    return MergedMarker(
        marker_id=f"DEL_{d.chrom}_{start}_{end}",
        chrom=d.chrom,
        start=start,
        end=end,
        sv_length=end - start,
        callers=tuple(sorted({d.caller, m.caller})),
        n_callers=2,
        qual_max=qual,
        contributing_record_ids=(d.raw_id, m.raw_id),
        genotypes=gts,
    )


# ----------------------------------------------------------------------
# Cohort utility: marker × sample genotype matrix.
# ----------------------------------------------------------------------


def to_genotype_matrix(
    markers: Sequence[MergedMarker],
) -> Dict[str, Dict[str, str]]:
    """Return ``{marker_id: {sample_id: genotype}}`` — the shape that
    feeds del_relatedness and the hemizygous-marker scorer."""
    return {m.marker_id: dict(m.genotypes) for m in markers}


def cohort_samples(markers: Sequence[MergedMarker]) -> List[str]:
    seen: List[str] = []
    seen_set: Set[str] = set()
    for m in markers:
        for s in m.genotypes:
            if s not in seen_set:
                seen.append(s)
                seen_set.add(s)
    return seen
