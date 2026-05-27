"""
HPP MVP 4 — transmission summaries (Tables C, D).

Per-dyad and per-triad rollups of Table A + Table B counts plus the
Mendelian consistency status.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import List, Optional, Sequence

from .gene_status import TableBRow
from .mendelian import (
    dyad_partial_consistent,
    is_called,
    offspring_carries_novel_allele,
    status_from_counts,
    triad_consistent,
)
from .project import TableARow
from .stage3_placeholder import DyadSegment, TriadSegment
from .variant_master import PlaceholderVariantMaster
from .vcf_lite import VariantRecord


# ----------------------------------------------------------------------
# Table C row (dyad).
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class TableCRow:
    dyad_id: str
    parent_sample_id: str
    offspring_sample_id: str
    n_segments_total: int
    n_segments_Gold: int
    n_segments_Silver: int
    n_segments_Bronze: int
    n_damaging_variants_in_parent: int
    n_damaging_variants_transmitted: int
    n_damaging_variants_resolved: int
    n_damaging_variants_unresolved: int
    n_genes_het_masked: int
    n_genes_hom_exposed_same_variant: int
    n_genes_compound_het_trans: int
    n_genes_compound_het_cis: int
    n_genes_partially_resolved: int
    mendelian_consistency_status: str
    mendelian_inconsistent_sites: int
    damaging_tier: str

    def to_dict(self) -> dict:
        return asdict(self)


TABLE_C_COLUMNS = [f.name for f in TableCRow.__dataclass_fields__.values()]


# ----------------------------------------------------------------------
# Table D row (triad).
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class TableDRow:
    triad_id: str
    paternal_sample_id: str
    maternal_sample_id: str
    offspring_sample_id: str
    n_segments_total: int
    n_segments_Gold: int
    n_segments_Silver: int
    n_segments_Bronze: int
    n_damaging_variants_in_paternal: int
    n_damaging_variants_in_maternal: int
    n_damaging_variants_transmitted: int
    n_damaging_variants_resolved: int
    n_damaging_variants_unresolved: int
    n_genes_het_masked: int
    n_genes_hom_exposed_same_variant: int
    n_genes_compound_het_trans: int
    n_genes_compound_het_cis: int
    n_genes_partially_resolved: int
    mendelian_consistency_status: str
    mendelian_inconsistent_sites: int
    mendelian_inconsistent_damaging_sites: int
    n_de_novo_candidates: int
    damaging_tier: str

    def to_dict(self) -> dict:
        return asdict(self)


TABLE_D_COLUMNS = [f.name for f in TableDRow.__dataclass_fields__.values()]


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------


def _segment_counts(segments) -> dict:
    out = {"total": len(segments), "Gold": 0, "Silver": 0, "Bronze": 0}
    for s in segments:
        out[s.segment_confidence] += 1
    return out


def _count_damaging_in_parent(
    parent_id: str,
    variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str,
) -> int:
    n = 0
    for v in variants:
        gt = v.genotypes.get(parent_id)
        if gt in (None, "./.", "0/0"):
            continue
        if variant_master.is_damaging(v.variant_id, tier=damaging_tier):
            n += 1
    return n


def _gene_status_counts(table_b_rows: Sequence[TableBRow]) -> dict:
    out = {
        "het_masked": 0,
        "hom_exposed_same_variant": 0,
        "compound_het_trans": 0,
        "compound_het_cis": 0,
        "partially_resolved": 0,
    }
    for r in table_b_rows:
        if r.predicted_gene_status in out:
            out[r.predicted_gene_status] += 1
    return out


def _projection_counts(
    table_a_rows: Sequence[TableARow],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str,
) -> dict:
    """transmitted / resolved / unresolved counts on damaging variants only."""
    transmitted = resolved = unresolved = 0
    for r in table_a_rows:
        if not variant_master.is_damaging(r.variant_id, tier=damaging_tier):
            continue
        if r.allele_state == "alt":
            transmitted += 1
        if r.projection_source in ("parent_homozygous",
                                   "parent_heterozygous_phased"):
            resolved += 1
        elif r.projection_source == "parent_heterozygous_unphased":
            unresolved += 1
    return {"transmitted": transmitted, "resolved": resolved,
            "unresolved": unresolved}


# ----------------------------------------------------------------------
# Dyad summary.
# ----------------------------------------------------------------------


def summarise_dyad(
    *,
    segments: Sequence[DyadSegment],
    table_a_rows: Sequence[TableARow],
    table_b_rows: Sequence[TableBRow],
    all_variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str = "T1",
) -> TableCRow:
    if not segments:
        raise ValueError("summarise_dyad: no segments")
    dyad_id = segments[0].dyad_id
    parent_id = segments[0].parent_sample_id
    offspring_id = segments[0].offspring_sample_id

    segc = _segment_counts(segments)
    projc = _projection_counts(table_a_rows, variant_master, damaging_tier)
    gsc = _gene_status_counts(table_b_rows)

    # Partial Mendelian — parent-homozygous subset.
    testable = 0
    inconsistent = 0
    for v in all_variants:
        gp = v.genotypes.get(parent_id)
        go = v.genotypes.get(offspring_id)
        verdict = dyad_partial_consistent(gp, go)
        if verdict is None:
            continue
        testable += 1
        if not verdict:
            inconsistent += 1
    status = status_from_counts(
        inconsistent=inconsistent, testable=testable, mode="dyad"
    )

    return TableCRow(
        dyad_id=dyad_id,
        parent_sample_id=parent_id,
        offspring_sample_id=offspring_id,
        n_segments_total=segc["total"],
        n_segments_Gold=segc["Gold"],
        n_segments_Silver=segc["Silver"],
        n_segments_Bronze=segc["Bronze"],
        n_damaging_variants_in_parent=_count_damaging_in_parent(
            parent_id, all_variants, variant_master, damaging_tier
        ),
        n_damaging_variants_transmitted=projc["transmitted"],
        n_damaging_variants_resolved=projc["resolved"],
        n_damaging_variants_unresolved=projc["unresolved"],
        n_genes_het_masked=gsc["het_masked"],
        n_genes_hom_exposed_same_variant=gsc["hom_exposed_same_variant"],
        n_genes_compound_het_trans=gsc["compound_het_trans"],
        n_genes_compound_het_cis=gsc["compound_het_cis"],
        n_genes_partially_resolved=gsc["partially_resolved"],
        mendelian_consistency_status=status,
        mendelian_inconsistent_sites=inconsistent,
        damaging_tier=damaging_tier,
    )


# ----------------------------------------------------------------------
# Triad summary.
# ----------------------------------------------------------------------


def summarise_triad(
    *,
    segments: Sequence[TriadSegment],
    table_a_rows: Sequence[TableARow],
    table_b_rows: Sequence[TableBRow],
    all_variants: Sequence[VariantRecord],
    variant_master: PlaceholderVariantMaster,
    damaging_tier: str = "T1",
) -> TableDRow:
    if not segments:
        raise ValueError("summarise_triad: no segments")
    triad_id = segments[0].triad_id
    paternal_id = segments[0].paternal_sample_id
    maternal_id = segments[0].maternal_sample_id
    offspring_id = segments[0].offspring_sample_id

    segc = _segment_counts(segments)
    projc = _projection_counts(table_a_rows, variant_master, damaging_tier)
    gsc = _gene_status_counts(table_b_rows)

    # Full Mendelian on both-parents-called sites.
    testable = 0
    inconsistent = 0
    inconsistent_damaging = 0
    de_novo = 0
    for v in all_variants:
        gp1 = v.genotypes.get(paternal_id)
        gp2 = v.genotypes.get(maternal_id)
        go = v.genotypes.get(offspring_id)
        if not (is_called(gp1) and is_called(gp2) and is_called(go)):
            continue
        testable += 1
        if not triad_consistent(gp1, gp2, go):
            inconsistent += 1
            if variant_master.is_damaging(v.variant_id, tier=damaging_tier):
                inconsistent_damaging += 1
            if offspring_carries_novel_allele(gp1, gp2, go):
                de_novo += 1
    status = status_from_counts(
        inconsistent=inconsistent, testable=testable, mode="triad"
    )

    return TableDRow(
        triad_id=triad_id,
        paternal_sample_id=paternal_id,
        maternal_sample_id=maternal_id,
        offspring_sample_id=offspring_id,
        n_segments_total=segc["total"],
        n_segments_Gold=segc["Gold"],
        n_segments_Silver=segc["Silver"],
        n_segments_Bronze=segc["Bronze"],
        n_damaging_variants_in_paternal=_count_damaging_in_parent(
            paternal_id, all_variants, variant_master, damaging_tier
        ),
        n_damaging_variants_in_maternal=_count_damaging_in_parent(
            maternal_id, all_variants, variant_master, damaging_tier
        ),
        n_damaging_variants_transmitted=projc["transmitted"],
        n_damaging_variants_resolved=projc["resolved"],
        n_damaging_variants_unresolved=projc["unresolved"],
        n_genes_het_masked=gsc["het_masked"],
        n_genes_hom_exposed_same_variant=gsc["hom_exposed_same_variant"],
        n_genes_compound_het_trans=gsc["compound_het_trans"],
        n_genes_compound_het_cis=gsc["compound_het_cis"],
        n_genes_partially_resolved=gsc["partially_resolved"],
        mendelian_consistency_status=status,
        mendelian_inconsistent_sites=inconsistent,
        mendelian_inconsistent_damaging_sites=inconsistent_damaging,
        n_de_novo_candidates=de_novo,
        damaging_tier=damaging_tier,
    )
