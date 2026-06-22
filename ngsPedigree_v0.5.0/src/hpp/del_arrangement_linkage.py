"""
Per-DEL × per-LRR arrangement-linkage classifier (the "situation 1"
table).

For each DEL marker that sits inside (or near) a candidate LRR, compute
the DEL allele frequency stratified by the host's PCA-band arrangement
class (HOM_REF / HET / HOM_INV). Classify each DEL as:

  - ``arrangement_0_marker``  DEL allele tracks HOM_REF arrangement
  - ``arrangement_1_marker``  DEL allele tracks HOM_INV arrangement
  - ``arrangement_1_marker_hom_depleted``
                              DEL tracks arrangement 1, but no/few
                              HOM_INV samples are observed (the
                              "situation 1" pattern — hemizygous-only
                              because the homokaryotype is absent or
                              depleted)
  - ``arrangement_0_marker_hom_depleted``  symmetric for arrangement 0
  - ``unlinked``               DEL frequency similar across all classes
  - ``ambiguous``              insufficient samples in one or more classes

The output is a long-table per (DEL, LRR) pair. Wording follows the
user's "do not overclaim lethality" rule: an absent homokaryotype is
reported as `hom_depleted`, never as `lethal`.

What linkage means here
-----------------------
A DEL allele lives on a specific haplotype of one of the two
arrangement classes. In a heterozygous (HET) sample, the DEL is
hemizygous (one copy carries it, one does not). In the homozygous
arrangement that carries the DEL, the DEL is homozygous (when present
on both copies). The host arrangement's PC1 band thus predicts the
DEL state at a tight expected ratio:

  Arrangement 1-linked DEL:
      HOM_REF samples (arr 0/0) → DEL absent      (rate ≈ 0)
      HET    samples (arr 0/1) → DEL hemizygous  (rate ≈ 0.5 per allele)
      HOM_INV samples (arr 1/1) → DEL homozygous  (rate ≈ 1.0 per allele)

So the per-class DEL allele frequency (counts of DEL alleles / 2 × n_samples)
should rise from 0 → 0.5 → 1.0 across the three classes for an
arrangement-1-linked DEL.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .del_inheritance import DelMarkerLocus


# DEL allele frequency per arrangement class
@dataclass(frozen=True)
class DelLinkageRecord:
    del_id: str
    lrr_id: str
    chrom: str
    del_start: int
    del_end: int
    # per-class DEL frequencies (DEL allele count / 2 * n_samples)
    n_hom_ref: int
    n_het: int
    n_hom_inv: int
    del_freq_hom_ref: Optional[float]
    del_freq_het: Optional[float]
    del_freq_hom_inv: Optional[float]
    interpretation: str
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def _classify_arrangement_linkage(
    freq_hom_ref: Optional[float],
    freq_het: Optional[float],
    freq_hom_inv: Optional[float],
    n_hom_ref: int,
    n_het: int,
    n_hom_inv: int,
    *,
    low_threshold: float = 0.10,
    high_threshold: float = 0.85,
    min_per_class: int = 2,
) -> Tuple[str, str]:
    """Return (interpretation, notes)."""
    # Need at least HET observations and at least one homokaryotype.
    n_classes_observed = sum(
        1 for n in (n_hom_ref, n_het, n_hom_inv) if n >= min_per_class
    )
    if n_classes_observed < 2:
        return ("ambiguous",
                "insufficient samples in one or more arrangement classes")

    fhr = freq_hom_ref if freq_hom_ref is not None else None
    fh = freq_het if freq_het is not None else None
    fhi = freq_hom_inv if freq_hom_inv is not None else None

    # Helper: is the value clearly "low" / "high" / "intermediate"?
    def _level(v):
        if v is None:
            return None
        if v <= low_threshold:
            return "low"
        if v >= high_threshold:
            return "high"
        return "intermediate"

    L_ref = _level(fhr)
    L_het = _level(fh)
    L_inv = _level(fhi)

    # arrangement-1 linked: low in hom_ref, high in hom_inv
    if L_ref == "low" and L_inv == "high":
        return ("arrangement_1_marker", "")
    if L_ref == "high" and L_inv == "low":
        return ("arrangement_0_marker", "")

    # arrangement-1 linked with hom_inv depleted: low in hom_ref,
    # heterozygous-only (intermediate in HET), hom_inv class absent
    if (L_ref == "low" and L_het == "intermediate"
            and (n_hom_inv < min_per_class or fhi is None)):
        return ("arrangement_1_marker_hom_depleted",
                "no/few HOM_INV samples observed; arrangement-linked DEL "
                "only seen hemizygous; cannot distinguish biological "
                "depletion from low frequency or technical hom_DEL "
                "miscall — interpret as depletion candidate, not lethality")
    if (L_inv == "low" and L_het == "intermediate"
            and (n_hom_ref < min_per_class or fhr is None)):
        return ("arrangement_0_marker_hom_depleted",
                "no/few HOM_REF samples observed; arrangement-linked DEL "
                "only seen hemizygous; cannot distinguish biological "
                "depletion from low frequency or technical hom_DEL "
                "miscall — interpret as depletion candidate, not lethality")

    # all three intermediate → unlinked / noisy
    if all(L in ("intermediate", None) for L in (L_ref, L_het, L_inv)):
        return ("unlinked",
                "DEL frequency similar across arrangement classes; "
                "DEL is not arrangement-linked")

    return ("ambiguous", "unexpected frequency pattern across classes")


# ----------------------------------------------------------------------
# Per-(DEL, LRR) classifier.
# ----------------------------------------------------------------------


def _gt_to_count(gt: str) -> Optional[int]:
    """Convert a genotype string to a count of alt alleles (DEL count)."""
    if gt == "0/0":
        return 0
    if gt == "0/1":
        return 1
    if gt == "1/1":
        return 2
    return None   # missing / unknown


def classify_del_linkage(
    *,
    del_id: str,
    lrr_id: str,
    chrom: str,
    del_start: int,
    del_end: int,
    sample_arrangement: Dict[str, str],     # sample_id → HOM_REF/HET/HOM_INV
    del_genotypes: Dict[str, str],          # sample_id → GT for this DEL
    low_threshold: float = 0.10,
    high_threshold: float = 0.85,
    min_per_class: int = 2,
) -> DelLinkageRecord:
    counts: Dict[str, List[int]] = {"HOM_REF": [], "HET": [], "HOM_INV": []}
    for sid, arr in sample_arrangement.items():
        if arr not in counts:
            continue
        gt = del_genotypes.get(sid)
        if gt is None:
            continue
        c = _gt_to_count(gt)
        if c is None:
            continue
        counts[arr].append(c)

    def _freq(xs: List[int]) -> Optional[float]:
        if not xs:
            return None
        return sum(xs) / (2 * len(xs))

    fhr = _freq(counts["HOM_REF"])
    fh = _freq(counts["HET"])
    fhi = _freq(counts["HOM_INV"])
    n_hr = len(counts["HOM_REF"])
    n_h = len(counts["HET"])
    n_hi = len(counts["HOM_INV"])

    interp, notes = _classify_arrangement_linkage(
        fhr, fh, fhi, n_hr, n_h, n_hi,
        low_threshold=low_threshold,
        high_threshold=high_threshold,
        min_per_class=min_per_class,
    )
    return DelLinkageRecord(
        del_id=del_id, lrr_id=lrr_id, chrom=chrom,
        del_start=del_start, del_end=del_end,
        n_hom_ref=n_hr, n_het=n_h, n_hom_inv=n_hi,
        del_freq_hom_ref=fhr, del_freq_het=fh, del_freq_hom_inv=fhi,
        interpretation=interp, notes=notes,
    )


# ----------------------------------------------------------------------
# Batch driver: classify every DEL inside each LRR.
# ----------------------------------------------------------------------


def classify_all(
    *,
    lrrs: Sequence,             # objects with .lrr_id, .chrom, .start, .end
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    sample_arrangement_by_lrr: Dict[str, Dict[str, str]],
    low_threshold: float = 0.10,
    high_threshold: float = 0.85,
    min_per_class: int = 2,
) -> List[DelLinkageRecord]:
    """Classify every DEL marker that falls inside each LRR using that
    LRR's per-sample arrangement assignment.

    ``sample_arrangement_by_lrr[lrr_id][sample_id]`` → "HOM_REF" | "HET" |
    "HOM_INV". (Different LRRs can have different arrangement-class
    assignments per sample.)
    """
    out: List[DelLinkageRecord] = []
    for lrr in lrrs:
        arr_map = sample_arrangement_by_lrr.get(lrr.lrr_id, {})
        for loc in loci:
            if loc.chrom != lrr.chrom:
                continue
            if not (lrr.start <= loc.midpoint < lrr.end):
                continue
            gts = genotype_matrix.get(loc.marker_id, {})
            rec = classify_del_linkage(
                del_id=loc.marker_id, lrr_id=lrr.lrr_id, chrom=loc.chrom,
                del_start=loc.start, del_end=loc.end,
                sample_arrangement=arr_map,
                del_genotypes=gts,
                low_threshold=low_threshold,
                high_threshold=high_threshold,
                min_per_class=min_per_class,
            )
            out.append(rec)
    return out


# ----------------------------------------------------------------------
# TSV writer.
# ----------------------------------------------------------------------


def write_linkage_tsv(path, records: Sequence[DelLinkageRecord]) -> None:
    cols = ["del_id", "lrr_id", "chrom", "del_start", "del_end",
            "n_hom_ref", "n_het", "n_hom_inv",
            "del_freq_hom_ref", "del_freq_het", "del_freq_hom_inv",
            "interpretation", "notes"]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    def _fmt(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.4f}"
        return str(v)

    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            row = [getattr(r, c) for c in cols]
            fh.write("\t".join(_fmt(v) for v in row) + "\n")
