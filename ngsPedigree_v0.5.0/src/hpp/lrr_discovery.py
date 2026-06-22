"""
De novo LRR / candidate-inversion discovery from DEL-marker correlations.

When no curated LRR list is available, search the merged DEL catalogue
for genomic regions where DEL genotypes are strongly haplotype-linked
across samples. The signal:

  - A real inversion creates two arrangement haplotypes; DELs sitting
    on one arrangement segregate together.
  - In a sliding window across each chromosome, that linkage shows up
    as elevated pairwise Pearson correlation between DEL genotype
    vectors (each vector is the per-sample DEL dosage 0/1/2).

The detector is intentionally conservative — its output is a candidate
LRR list, the same shape that `--list_of_LRR` accepts. False positives
are caught downstream by the family-based enrichment (bloc 13) and
Mendelian segregation analysis (bloc 14).
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .del_inheritance import DelMarkerLocus


@dataclass(frozen=True)
class CandidateLRR:
    lrr_id: str
    chrom: str
    start: int
    end: int
    n_markers: int
    mean_pairwise_correlation: float
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# Pearson correlation (stdlib).
# ----------------------------------------------------------------------


def _gt_to_count(gt: str) -> Optional[int]:
    if gt == "0/0":
        return 0
    if gt == "0/1":
        return 1
    if gt == "1/1":
        return 2
    return None


def pearson_corr(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    n = len(xs)
    if n < 3:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    den = math.sqrt(var_x * var_y)
    return num / den if den > 0 else None


def marker_pair_correlation(
    marker_a: str,
    marker_b: str,
    genotype_matrix: Dict[str, Dict[str, str]],
    samples: Sequence[str],
) -> Optional[float]:
    a_calls = genotype_matrix.get(marker_a, {})
    b_calls = genotype_matrix.get(marker_b, {})
    xs: List[float] = []
    ys: List[float] = []
    for s in samples:
        ca = _gt_to_count(a_calls.get(s, "./."))
        cb = _gt_to_count(b_calls.get(s, "./."))
        if ca is None or cb is None:
            continue
        xs.append(ca)
        ys.append(cb)
    return pearson_corr(xs, ys)


def window_mean_correlation(
    marker_ids: Sequence[str],
    genotype_matrix: Dict[str, Dict[str, str]],
    samples: Sequence[str],
) -> Tuple[Optional[float], int]:
    """Mean pairwise |correlation| over all (i, j) marker pairs in the
    window. Returns (mean, n_pairs_evaluated)."""
    if len(marker_ids) < 2:
        return (None, 0)
    total = 0.0
    n = 0
    for i, m_a in enumerate(marker_ids):
        for m_b in marker_ids[i + 1:]:
            c = marker_pair_correlation(m_a, m_b, genotype_matrix, samples)
            if c is None:
                continue
            total += abs(c)
            n += 1
    if n == 0:
        return (None, 0)
    return (total / n, n)


# ----------------------------------------------------------------------
# Sliding-window scanner.
# ----------------------------------------------------------------------


def discover_candidate_lrrs(
    *,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    samples: Sequence[str],
    window_size: int = 1_000_000,
    step: Optional[int] = None,
    min_markers_per_window: int = 4,
    correlation_threshold: float = 0.50,
    merge_adjacent: bool = True,
    max_merge_gap: int = 500_000,
) -> List[CandidateLRR]:
    """Walk each chromosome in sliding windows of ``window_size`` and
    emit candidate LRRs where the mean pairwise |correlation| between
    DEL markers exceeds ``correlation_threshold`` with at least
    ``min_markers_per_window`` markers.

    Adjacent (or nearly adjacent, within ``max_merge_gap``) windows are
    merged into one candidate LRR interval when ``merge_adjacent`` is True.
    """
    step = step or window_size // 2

    # Group loci by chromosome.
    by_chrom: Dict[str, List[DelMarkerLocus]] = {}
    for l in loci:
        by_chrom.setdefault(l.chrom, []).append(l)
    for v in by_chrom.values():
        v.sort(key=lambda x: x.midpoint)

    candidates: List[CandidateLRR] = []
    counter = 0
    for chrom, chr_loci in sorted(by_chrom.items()):
        if not chr_loci:
            continue
        chrom_start = chr_loci[0].midpoint
        chrom_end = chr_loci[-1].midpoint + 1
        win_lo = max(0, chrom_start - 1)
        while win_lo < chrom_end:
            win_hi = win_lo + window_size
            in_window = [l.marker_id for l in chr_loci
                         if win_lo <= l.midpoint < win_hi]
            if len(in_window) >= min_markers_per_window:
                mean_corr, n_pairs = window_mean_correlation(
                    in_window, genotype_matrix, samples,
                )
                if mean_corr is not None and mean_corr >= correlation_threshold:
                    counter += 1
                    candidates.append(CandidateLRR(
                        lrr_id=f"cLRR_{counter:04d}",
                        chrom=chrom, start=win_lo, end=win_hi,
                        n_markers=len(in_window),
                        mean_pairwise_correlation=mean_corr,
                    ))
            win_lo += step

    if not merge_adjacent:
        return candidates

    # Merge overlapping or near-adjacent candidates per chromosome.
    merged: List[CandidateLRR] = []
    by_chrom_cands: Dict[str, List[CandidateLRR]] = {}
    for c in candidates:
        by_chrom_cands.setdefault(c.chrom, []).append(c)
    for chrom in sorted(by_chrom_cands):
        cs = sorted(by_chrom_cands[chrom], key=lambda x: x.start)
        cur = cs[0]
        for nxt in cs[1:]:
            if nxt.start <= cur.end + max_merge_gap:
                # merge
                new_start = min(cur.start, nxt.start)
                new_end = max(cur.end, nxt.end)
                new_corr = (cur.mean_pairwise_correlation
                            + nxt.mean_pairwise_correlation) / 2
                new_n = cur.n_markers + nxt.n_markers
                cur = CandidateLRR(
                    lrr_id=cur.lrr_id, chrom=chrom,
                    start=new_start, end=new_end,
                    n_markers=new_n,
                    mean_pairwise_correlation=new_corr,
                    notes="merged",
                )
            else:
                merged.append(cur)
                cur = nxt
        merged.append(cur)

    # Re-id sequentially after merging.
    out: List[CandidateLRR] = []
    for i, c in enumerate(merged, start=1):
        out.append(CandidateLRR(
            lrr_id=f"cLRR_{i:04d}",
            chrom=c.chrom, start=c.start, end=c.end,
            n_markers=c.n_markers,
            mean_pairwise_correlation=c.mean_pairwise_correlation,
            notes=c.notes,
        ))
    return out


def write_candidate_lrr_tsv(path, cands: Sequence[CandidateLRR]) -> None:
    """Output in the same shape that `--list_of_LRR` consumes (plus
    extra QC columns)."""
    cols = ["lrr_id", "chrom", "start", "end",
            "n_markers", "mean_pairwise_correlation", "notes"]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in cands:
            row = [c.lrr_id, c.chrom, c.start, c.end,
                   c.n_markers, f"{c.mean_pairwise_correlation:.4f}",
                   c.notes]
            fh.write("\t".join(str(v) for v in row) + "\n")
