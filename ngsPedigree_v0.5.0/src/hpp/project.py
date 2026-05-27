"""
HPP MVP 2 — dyad and triad projection into Table A.

Implements ``project_dyad_to_offspring`` and ``project_triad_to_offspring``
from SPEC_HPP.md §10.2 and §6 Step 3. Emits one
``TableARow`` per (offspring, variant, hap_copy).

Damaging-or-synonymous filter follows the spec: variants the
``VariantMasterAdapter`` knows about as damaging (under the selected
KBC §1.8 tier) OR as synonymous controls are kept; everything else is
dropped.

Confidence rules per SPEC §5.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable, List, Optional, Sequence

from .parental_haps import ParentHapVariants
from .stage3_placeholder import DyadSegment, TriadSegment
from .variant_master import PlaceholderVariantMaster, VariantAnnotation
from .vcf_lite import VariantRecord


# ----------------------------------------------------------------------
# Table A row.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class TableARow:
    offspring_sample_id: str
    relationship_id: str
    relationship_type: str   # "dyad" | "triad"
    variant_id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene_id: Optional[str]
    transcript_id: Optional[str]
    consequence: Optional[str]
    impact: Optional[str]
    sift_class: Optional[str]
    vesm_llr: Optional[float]
    splice_subclass: Optional[str]
    hap_copy: str            # hap_from_P1 | hap_from_P2 | unassigned
    allele_state: str        # ref | alt | unknown
    projection_source: str   # parent_homozygous | parent_heterozygous_phased | parent_heterozygous_unphased
    segment_confidence: str  # Gold | Silver | Bronze
    confidence: str          # high | medium | low | unresolved

    def to_dict(self) -> dict:
        return asdict(self)


TABLE_A_COLUMNS = [f.name for f in TableARow.__dataclass_fields__.values()]


# ----------------------------------------------------------------------
# Confidence — SPEC §5.
# ----------------------------------------------------------------------


def composite_confidence(projection_source: str, segment_confidence: str) -> str:
    if segment_confidence == "Bronze":
        return "low"
    if projection_source == "parent_homozygous":
        return "high"
    if projection_source == "parent_heterozygous_phased":
        return "high" if segment_confidence == "Gold" else "medium"
    # parent_heterozygous_unphased
    return "unresolved"


# ----------------------------------------------------------------------
# Damaging-or-synonymous filter (SPEC §10.2).
# ----------------------------------------------------------------------


def _ann_is_synonymous(ann: Optional[VariantAnnotation]) -> bool:
    return ann is not None and ann.consequence == "synonymous_variant"


def _keep_variant(
    variant_id: str,
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str,
) -> bool:
    ann = variant_master.lookup(variant_id)
    if ann is None:
        return False
    if _ann_is_synonymous(ann):
        return True
    return variant_master.is_damaging(variant_id, tier=damaging_tier)


def _ann_fields(
    variant_id: str, variant_master: PlaceholderVariantMaster
) -> dict:
    ann = variant_master.lookup(variant_id)
    if ann is None:
        return dict(gene_id=None, transcript_id=None, consequence=None,
                    impact=None, sift_class=None, vesm_llr=None,
                    splice_subclass=None)
    return dict(
        gene_id=ann.gene_id,
        transcript_id=ann.transcript_id,
        consequence=ann.consequence,
        impact=ann.impact,
        sift_class=ann.sift_class,
        vesm_llr=ann.vesm_llr,
        splice_subclass=ann.splice_subclass,
    )


# ----------------------------------------------------------------------
# Variant-in-segment helper.
# ----------------------------------------------------------------------


def _variants_in_segment(
    variants: Iterable[VariantRecord], chrom: str, seg_start: int, seg_end: int
) -> List[VariantRecord]:
    """0-based half-open [seg_start, seg_end) on `chrom`. VCF pos is 1-based."""
    out = []
    for v in variants:
        if v.chrom != chrom:
            continue
        pos0 = v.pos - 1
        if seg_start <= pos0 < seg_end:
            out.append(v)
    return out


# ----------------------------------------------------------------------
# Single-parent projection — used by both dyad and triad code paths.
# ----------------------------------------------------------------------


def _project_one_parent(
    *,
    offspring_sample_id: str,
    relationship_id: str,
    relationship_type: str,
    parent_sample_id: str,
    parent_haps: ParentHapVariants,
    inherited_hap: str,           # "1" | "2" | "ambiguous"
    segment_confidence: str,      # Gold | Silver | Bronze
    hap_copy: str,                # hap_from_P1 | hap_from_P2
    segment_variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str,
) -> List[TableARow]:
    rows: List[TableARow] = []
    for v in segment_variants:
        if not _keep_variant(v.variant_id, variant_master, damaging_tier):
            continue
        gt = v.genotypes.get(parent_sample_id)
        if gt in (None, "./.", "0/0"):
            continue

        ann = _ann_fields(v.variant_id, variant_master)

        if gt == "1/1":
            projection_source = "parent_homozygous"
            allele_state = "alt"
        elif gt in ("0/1", "1/0"):
            if inherited_hap == "ambiguous":
                projection_source = "parent_heterozygous_unphased"
                allele_state = "unknown"
                hap_copy_eff = "unassigned"
            elif v in parent_haps.unphased:
                # parent is het, no phase info for this variant
                projection_source = "parent_heterozygous_unphased"
                allele_state = "unknown"
                hap_copy_eff = "unassigned"
            else:
                projection_source = "parent_heterozygous_phased"
                on_inherited = parent_haps.carries_on_hap(
                    v.variant_id, int(inherited_hap)
                )
                allele_state = "alt" if on_inherited else "ref"
                hap_copy_eff = hap_copy
        else:
            raise ValueError(
                f"_project_one_parent: unrecognised GT {gt!r} for "
                f"{parent_sample_id} at {v.variant_id}"
            )

        # Resolved-hap rows keep the requested hap_copy; ambiguous rows
        # overwrite to "unassigned" above.
        if gt == "1/1":
            hap_copy_eff = hap_copy

        rows.append(TableARow(
            offspring_sample_id=offspring_sample_id,
            relationship_id=relationship_id,
            relationship_type=relationship_type,
            variant_id=v.variant_id,
            chrom=v.chrom,
            pos=v.pos,
            ref=v.ref,
            alt=v.alt,
            hap_copy=hap_copy_eff,
            allele_state=allele_state,
            projection_source=projection_source,
            segment_confidence=segment_confidence,
            confidence=composite_confidence(projection_source, segment_confidence),
            **ann,
        ))
    return rows


# ----------------------------------------------------------------------
# Public entry points.
# ----------------------------------------------------------------------


def project_dyad_to_offspring(
    *,
    segments: Sequence[DyadSegment],
    parent_haps: ParentHapVariants,
    all_variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str = "T1",
) -> List[TableARow]:
    """Project a single parent onto its offspring across one dyad's segments.

    Only emits ``hap_from_P1`` rows (the known parent side). The
    unknown-parent side stays unrepresented at MVP 2; downstream
    consumers should not assume ``hap_from_P2`` rows exist for dyad
    projections.
    """
    if not segments:
        return []
    dyad_id = segments[0].dyad_id
    offspring_sample_id = segments[0].offspring_sample_id
    parent_sample_id = segments[0].parent_sample_id

    rows: List[TableARow] = []
    for seg in segments:
        if seg.dyad_id != dyad_id:
            raise ValueError(
                f"project_dyad_to_offspring: mixed dyad_ids: "
                f"{dyad_id} vs {seg.dyad_id}"
            )
        seg_vars = _variants_in_segment(
            all_variants, seg.chrom, seg.seg_start, seg.seg_end
        )
        rows.extend(_project_one_parent(
            offspring_sample_id=offspring_sample_id,
            relationship_id=dyad_id,
            relationship_type="dyad",
            parent_sample_id=parent_sample_id,
            parent_haps=parent_haps,
            inherited_hap=seg.parental_hap_inherited,
            segment_confidence=seg.segment_confidence,
            hap_copy="hap_from_P1",
            segment_variants=seg_vars,
            variant_master=variant_master,
            damaging_tier=damaging_tier,
        ))
    return rows


def project_triad_to_offspring(
    *,
    segments: Sequence[TriadSegment],
    paternal_haps: ParentHapVariants,
    maternal_haps: ParentHapVariants,
    all_variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str = "T1",
) -> List[TableARow]:
    """Project both parents onto the offspring across one triad's segments.

    Paternal projections emit ``hap_copy = hap_from_P1``, maternal
    emit ``hap_copy = hap_from_P2``.
    """
    if not segments:
        return []
    triad_id = segments[0].triad_id
    offspring_sample_id = segments[0].offspring_sample_id
    paternal_sample_id = segments[0].paternal_sample_id
    maternal_sample_id = segments[0].maternal_sample_id

    rows: List[TableARow] = []
    for seg in segments:
        if seg.triad_id != triad_id:
            raise ValueError(
                f"project_triad_to_offspring: mixed triad_ids: "
                f"{triad_id} vs {seg.triad_id}"
            )
        seg_vars = _variants_in_segment(
            all_variants, seg.chrom, seg.seg_start, seg.seg_end
        )
        rows.extend(_project_one_parent(
            offspring_sample_id=offspring_sample_id,
            relationship_id=triad_id,
            relationship_type="triad",
            parent_sample_id=paternal_sample_id,
            parent_haps=paternal_haps,
            inherited_hap=seg.paternal_hap_inherited,
            segment_confidence=seg.segment_confidence,
            hap_copy="hap_from_P1",
            segment_variants=seg_vars,
            variant_master=variant_master,
            damaging_tier=damaging_tier,
        ))
        rows.extend(_project_one_parent(
            offspring_sample_id=offspring_sample_id,
            relationship_id=triad_id,
            relationship_type="triad",
            parent_sample_id=maternal_sample_id,
            parent_haps=maternal_haps,
            inherited_hap=seg.maternal_hap_inherited,
            segment_confidence=seg.segment_confidence,
            hap_copy="hap_from_P2",
            segment_variants=seg_vars,
            variant_master=variant_master,
            damaging_tier=damaging_tier,
        ))
    return rows
