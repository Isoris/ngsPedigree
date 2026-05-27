"""
Parental haplotype variant builder — HPP MVP 1.

Implements ``build_parental_hap_variants`` from SPEC_HPP.md §10.1.

For each parent in a dyad / triad, walks the joint VCF and groups the
parent's variant calls into two per-hap lists:

  - 1/1 sites land on BOTH hap 1 and hap 2 (homozygous).
  - 0/1 sites land on hap 1 OR hap 2 according to ``parent_phase``;
    sites with no phase entry are returned as ``unphased`` so the
    projection step can emit them with
    ``projection_source = parent_heterozygous_unphased``.
  - 0/0 and ./. sites are skipped.

The result is a ``ParentHapVariants`` record with:
  - ``hap1``: list of VariantRecord assigned to hap 1
  - ``hap2``: list of VariantRecord assigned to hap 2
  - ``unphased``: list of VariantRecord that are heterozygous and
    have no parent_phase entry
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List

from .stage3_placeholder import PhaseKey
from .vcf_lite import VariantRecord


@dataclass
class ParentHapVariants:
    parent_sample_id: str
    hap1: List[VariantRecord] = field(default_factory=list)
    hap2: List[VariantRecord] = field(default_factory=list)
    unphased: List[VariantRecord] = field(default_factory=list)

    def n_total(self) -> int:
        return len(self.hap1) + len(self.hap2) + len(self.unphased)

    def carries_on_hap(self, variant_id: str, hap_no: int) -> bool:
        lst = self.hap1 if hap_no == 1 else self.hap2
        return any(v.variant_id == variant_id for v in lst)


def build_parental_hap_variants(
    parent_sample_id: str,
    variants: Iterable[VariantRecord],
    parent_phase: Dict[PhaseKey, str],
) -> ParentHapVariants:
    out = ParentHapVariants(parent_sample_id=parent_sample_id)
    for v in variants:
        gt = v.genotypes.get(parent_sample_id)
        if gt in (None, "./.", "0/0"):
            continue
        if gt == "1/1":
            out.hap1.append(v)
            out.hap2.append(v)
            continue
        if gt in ("0/1", "1/0"):
            hap = parent_phase.get((parent_sample_id, v.variant_id))
            if hap == "1":
                out.hap1.append(v)
            elif hap == "2":
                out.hap2.append(v)
            else:
                out.unphased.append(v)
            continue
        raise ValueError(
            f"build_parental_hap_variants: unrecognised GT {gt!r} for "
            f"{parent_sample_id} at {v.variant_id}"
        )
    return out


def build_for_parents(
    parent_sample_ids: Iterable[str],
    variants: Iterable[VariantRecord],
    parent_phase: Dict[PhaseKey, str],
) -> Dict[str, ParentHapVariants]:
    variants = list(variants)
    return {
        pid: build_parental_hap_variants(pid, variants, parent_phase)
        for pid in parent_sample_ids
    }
