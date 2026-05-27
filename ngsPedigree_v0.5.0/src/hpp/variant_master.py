"""
variant_master_scored.tsv adapter — MVP 2 stub.

HPP consumes a small subset of MODULE_CONSERVATION STEP 16 columns:
  - variant_id        (chrom:pos:ref:alt)
  - gene_id, transcript_id
  - consequence, impact   (SnpEff)
  - sift_class            (SIFT4G)
  - vesm_llr              (VESM)
  - splice_subclass       (splice module)

Damaging-set definition follows KBC §1.8 three tiers. HPP receives the
*tier selection* as a runtime parameter; the same data table can be
classified under T1, T2, or T3.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

# KBC §1.8 — Tier 1 high-confidence LoF consequence tokens (SnpEff).
TIER_1_CONSEQUENCES = frozenset({
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "splice_donor_variant",
    "splice_acceptor_variant",
})

# Tier 2 adds validated splice Class A subclasses (read from
# splice_subclass column). The Class A enum is project-local and
# checked against the §1.9 six-check splice-validation gate before use.
TIER_2_SPLICE_CLASS_A = frozenset({
    "splice_donor_5th_base_classA",
    "splice_branch_classA",
    "polypyrimidine_tract_classA",
})

# Tier 3 — model-dependent missense. SIFT4G "deleterious" plus VESM LLR
# below a hard threshold.
TIER_3_VESM_THRESHOLD = -7.0
TIER_3_SIFT_DAMAGING = frozenset({"deleterious", "deleterious_low_confidence"})


@dataclass(frozen=True)
class VariantAnnotation:
    variant_id: str
    gene_id: Optional[str]
    transcript_id: Optional[str]
    consequence: Optional[str]
    impact: Optional[str]
    sift_class: Optional[str]
    vesm_llr: Optional[float]
    splice_subclass: Optional[str]


class VariantMasterNotReadyError(NotImplementedError):
    pass


class PlaceholderVariantMaster:
    """In-memory placeholder used by HPP unit tests and synthetic fixtures.

    The real adapter reads variant_master_scored.tsv (TSV with the columns
    listed in MODULE_CONSERVATION STEP 16). MVP 2 will implement it
    against the real table.
    """

    def __init__(self, records: Optional[dict[str, VariantAnnotation]] = None):
        self._records: dict[str, VariantAnnotation] = records or {}

    def add(self, ann: VariantAnnotation) -> None:
        self._records[ann.variant_id] = ann

    def lookup(self, variant_id: str) -> Optional[VariantAnnotation]:
        return self._records.get(variant_id)

    def is_damaging(self, variant_id: str, tier: str = "T1") -> bool:
        ann = self._records.get(variant_id)
        if ann is None:
            return False
        if tier not in {"T1", "T2", "T3"}:
            raise ValueError(f"unknown damaging tier: {tier!r}")

        if ann.consequence in TIER_1_CONSEQUENCES:
            return True
        if tier == "T1":
            return False

        if ann.splice_subclass in TIER_2_SPLICE_CLASS_A:
            return True
        if tier == "T2":
            return False

        # T3 — model-dependent missense
        if ann.sift_class in TIER_3_SIFT_DAMAGING:
            return True
        if ann.vesm_llr is not None and ann.vesm_llr <= TIER_3_VESM_THRESHOLD:
            return True
        return False


def load_variant_master(path) -> PlaceholderVariantMaster:
    raise VariantMasterNotReadyError(
        "variant_master.load_variant_master: real-table loader is MVP 2. "
        "Construct PlaceholderVariantMaster directly for unit tests."
    )
