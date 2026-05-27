"""
HPP MVP 3 — per-(offspring, gene, segment) status classifier.

Implements SPEC_HPP.md §7 seven-class enum on top of Table A rows.

Inputs:
  - Table A rows (output of project.py)
  - Inheritance-map segments (for seg_start / seg_end lookup, since
    Table A does not carry segment coordinates by spec)
  - Damaging-tier (T1 | T2 | T3) — must match the tier under which
    Table A was produced

Output:
  - List[TableBRow]; column list matches schemas/B_hpp_offspring_gene_status.schema.json
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .project import TableARow
from .stage3_placeholder import DyadSegment, TriadSegment
from .variant_master import PlaceholderVariantMaster


# ----------------------------------------------------------------------
# Table B row.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class TableBRow:
    offspring_sample_id: str
    gene_id: str
    transcript_id: Optional[str]
    chrom: str
    seg_start: int
    seg_end: int
    hap_from_P1_damaging_variants: str  # ";"-joined variant_ids
    hap_from_P2_damaging_variants: str
    predicted_gene_status: str
    inside_inversion: Optional[str]            # MVP 3: always None
    inversion_karyotype_class: Optional[str]   # MVP 3: always None
    confidence: str
    damaging_tier: str

    def to_dict(self) -> dict:
        return asdict(self)


TABLE_B_COLUMNS = [f.name for f in TableBRow.__dataclass_fields__.values()]


# ----------------------------------------------------------------------
# Confidence ordering — SPEC §5.
# ----------------------------------------------------------------------

_CONFIDENCE_RANK = {"high": 0, "medium": 1, "low": 2, "unresolved": 3}
_CONFIDENCE_INV = {v: k for k, v in _CONFIDENCE_RANK.items()}


def _worst_confidence(values: Iterable[str]) -> str:
    """Return the *worst* (lowest-quality) confidence in the iterable."""
    seen = list(values)
    if not seen:
        return "high"
    return _CONFIDENCE_INV[max(_CONFIDENCE_RANK[v] for v in seen)]


# ----------------------------------------------------------------------
# Segment lookup.
# ----------------------------------------------------------------------


SegmentRange = Tuple[str, int, int]   # chrom, seg_start, seg_end


def _segment_for_pos(
    segment_ranges: Sequence[SegmentRange], chrom: str, pos: int
) -> Optional[SegmentRange]:
    """Find the unique segment containing 1-based VCF pos (=> pos0=pos-1)."""
    pos0 = pos - 1
    for r in segment_ranges:
        if r[0] == chrom and r[1] <= pos0 < r[2]:
            return r
    return None


def _ranges_from_segments(
    segments: Sequence[DyadSegment] | Sequence[TriadSegment],
) -> List[SegmentRange]:
    out: List[SegmentRange] = []
    for s in segments:
        out.append((s.chrom, s.seg_start, s.seg_end))
    return out


# ----------------------------------------------------------------------
# Seven-class rule — SPEC §7.
# ----------------------------------------------------------------------


def _classify(d1: List[str], d2: List[str], unresolved: List[str]) -> str:
    """Pure classification from variant-id sets on each copy."""
    if unresolved:
        return "partially_resolved" if (d1 or d2) else "unresolved"
    if not d1 and not d2:
        return "reference_like"
    if d1 and not d2:
        return "het_masked" if len(d1) == 1 else "compound_het_cis"
    if d2 and not d1:
        return "het_masked" if len(d2) == 1 else "compound_het_cis"
    # both copies hit
    if len(d1) == 1 and len(d2) == 1 and d1[0] == d2[0]:
        return "hom_exposed_same_variant"
    return "compound_het_trans"


# ----------------------------------------------------------------------
# Public entry point.
# ----------------------------------------------------------------------


def classify_gene_status(
    *,
    table_a_rows: Sequence[TableARow],
    segments: Sequence[DyadSegment] | Sequence[TriadSegment],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str = "T1",
) -> List[TableBRow]:
    """Classify every (offspring, gene, segment) group from Table A.

    Synonymous-only rows contribute to *confidence* aggregation but not
    to d1 / d2 / unresolved sets; a group consisting entirely of
    synonymous rows classifies as ``reference_like``.

    Rows with ``gene_id is None`` are silently skipped (no gene to
    classify against).
    """
    ranges = _ranges_from_segments(segments)

    # group rows by (offspring, gene, segment_range)
    groups: Dict[Tuple[str, str, SegmentRange], List[TableARow]] = {}
    for row in table_a_rows:
        if row.gene_id is None:
            continue
        seg = _segment_for_pos(ranges, row.chrom, row.pos)
        if seg is None:
            continue
        key = (row.offspring_sample_id, row.gene_id, seg)
        groups.setdefault(key, []).append(row)

    out: List[TableBRow] = []
    for (offspring, gene, seg), rows in groups.items():
        d1: List[str] = []
        d2: List[str] = []
        unresolved: List[str] = []

        for r in rows:
            if not variant_master.is_damaging(r.variant_id, tier=damaging_tier):
                continue
            if r.hap_copy == "unassigned":
                unresolved.append(r.variant_id)
                continue
            if r.allele_state != "alt":
                continue   # ref-on-this-hap row — not carrying the damaging allele
            if r.hap_copy == "hap_from_P1":
                d1.append(r.variant_id)
            elif r.hap_copy == "hap_from_P2":
                d2.append(r.variant_id)

        status = _classify(d1, d2, unresolved)

        # Confidence from contributing rows. If the group is reference_like
        # (no damaging variants), aggregate over all rows so synonymous
        # confidence is reflected.
        contributing = [
            r for r in rows
            if variant_master.is_damaging(r.variant_id, tier=damaging_tier)
        ] or rows
        conf = _worst_confidence(r.confidence for r in contributing)

        # Use the first row's gene/transcript metadata.
        rep = rows[0]
        out.append(TableBRow(
            offspring_sample_id=offspring,
            gene_id=gene,
            transcript_id=rep.transcript_id,
            chrom=seg[0],
            seg_start=seg[1],
            seg_end=seg[2],
            hap_from_P1_damaging_variants=";".join(d1),
            hap_from_P2_damaging_variants=";".join(d2),
            predicted_gene_status=status,
            inside_inversion=None,
            inversion_karyotype_class=None,
            confidence=conf,
            damaging_tier=damaging_tier,
        ))
    return out
