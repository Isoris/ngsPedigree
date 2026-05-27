"""
HPP MVP 4 — Mendelian consistency.

Triads:
  Full test at every site where both parents have called genotypes.
  Offspring must inherit one allele from each parent.

Dyads:
  Partial test — only sites where the parent is homozygous can be
  checked (the inherited allele is fixed).

Genotype-string convention follows the VCF reader: "0/0", "0/1",
"1/1", "./." (after vcf_lite normalisation of phased "|" → "/").
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Set


def _alleles(gt: str) -> Set[str]:
    """{'0'} for hom-ref, {'0', '1'} for het, {'1'} for hom-alt."""
    return set(gt.split("/"))


def mendelian_expected_set(gp1: str, gp2: str) -> Set[str]:
    """Return canonical /-joined sorted offspring genotypes consistent with
    a single allele drawn from each parent."""
    expected: Set[str] = set()
    for a1 in _alleles(gp1):
        for a2 in _alleles(gp2):
            expected.add("/".join(sorted([a1, a2])))
    return expected


def is_called(gt: Optional[str]) -> bool:
    return gt is not None and gt != "./."


def triad_consistent(gp1: str, gp2: str, go: str) -> bool:
    if not (is_called(gp1) and is_called(gp2) and is_called(go)):
        return True   # uncallable sites are not inconsistencies
    return go in mendelian_expected_set(gp1, gp2)


def offspring_carries_novel_allele(gp1: str, gp2: str, go: str) -> bool:
    if not (is_called(gp1) and is_called(gp2) and is_called(go)):
        return False
    parental = _alleles(gp1) | _alleles(gp2)
    return bool(_alleles(go) - parental)


def dyad_partial_consistent(parent_gt: str, offspring_gt: str) -> Optional[bool]:
    """Partial Mendelian check on parent-homozygous subset only.

    Returns:
      True  — consistent
      False — inconsistent (offspring missing the parent-fixed allele)
      None  — untestable (parent het or either GT missing)
    """
    if not (is_called(parent_gt) and is_called(offspring_gt)):
        return None
    if parent_gt == "1/1":
        return "1" in _alleles(offspring_gt)
    if parent_gt == "0/0":
        return "0" in _alleles(offspring_gt)
    return None


# ----------------------------------------------------------------------
# Status from counts.
# ----------------------------------------------------------------------


def status_from_counts(
    *, inconsistent: int, testable: int, mode: str = "triad"
) -> str:
    """Map inconsistency counts to the enum used in Tables C and D.

    Triad: testable across all both-parents-called sites.
    Dyad : testable across parent-homozygous sites only.

    Default thresholds (HANDOFF open question #5; calibrate later):
      0 inconsistencies            → pass
      1-2 inconsistencies          → warn
      ≥ 3 inconsistencies          → fail
      0 testable sites             → untestable
    """
    if testable == 0:
        return "untestable"
    if inconsistent == 0:
        return "pass"
    if inconsistent <= 2:
        return "warn"
    return "fail"
