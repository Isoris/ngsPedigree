"""
Cross-LRR comparative analysis (situation 1 → genome-wide test).

For every candidate LRR, summarise the cohort genotype distribution
(arrangement class counts, FIS, het excess, hom/het ratios), compute a
between-arrangement divergence proxy (dXY across the DEL markers that
fall inside the LRR), and fit a simple linear regression across LRRs:

    log(min(hom_het_ratio)) ~ dXY_between_arrangements
    FIS                     ~ dXY_between_arrangements

The hypothesis is the user-stated one: more divergent LRRs should show
stronger homokaryotype depletion and heterokaryotype enrichment under
progressive pseudo-overdominance. dXY is treated as a **relative
divergence proxy**, not an absolute age estimate (dXY is influenced
by ancestral diversity, mutation rate, Ne, and selection — wording
follows the user's thesis-safe framing).

Tests follow the bloc-13 / 14 honesty rules:
  - any per-LRR metric undefined → null (no fabrication)
  - regression p-value is asymptotic (normal-approximation t-test);
    flagged as such for n_LRRs < 30
  - never claims age caused the depletion; the result supports an
    evolutionary model, it does not prove cause-and-effect.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ----------------------------------------------------------------------
# Records.
# ----------------------------------------------------------------------


@dataclass
class LRRComparativeMetrics:
    lrr_id: str
    chrom: str
    start: int
    end: int
    length: int
    n_samples_called: int
    n_hom0: int
    n_het: int
    n_hom1: int
    freq_hom0: Optional[float]
    freq_het: Optional[float]
    freq_hom1: Optional[float]
    p_allele_1: Optional[float]              # frequency of arrangement-1 allele
    H_obs: Optional[float]
    H_exp: Optional[float]
    FIS: Optional[float]
    heterokaryotype_enrichment: Optional[float]   # H_obs / H_exp
    hom0_het_ratio: Optional[float]
    hom1_het_ratio: Optional[float]
    min_hom_het_ratio: Optional[float]
    missing_hom_class: Optional[str]              # "HOM_REF" | "HOM_INV" | None
    n_markers_for_dxy: int
    dxy_between_arrangements: Optional[float]
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ComparativeRegression:
    y_name: str
    x_name: str
    n_lrrs: int
    slope: Optional[float]
    intercept: Optional[float]
    r: Optional[float]
    r_squared: Optional[float]
    t_stat: Optional[float]
    p_value_asymptotic: Optional[float]
    asymptotic_caveat: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# FIS, het excess.
# ----------------------------------------------------------------------


def genotype_class_counts(
    arrangement_by_sample: Dict[str, str],
) -> Tuple[int, int, int]:
    n_hom0 = sum(1 for a in arrangement_by_sample.values() if a == "HOM_REF")
    n_het = sum(1 for a in arrangement_by_sample.values() if a == "HET")
    n_hom1 = sum(1 for a in arrangement_by_sample.values() if a == "HOM_INV")
    return n_hom0, n_het, n_hom1


def compute_fis(n_hom0: int, n_het: int, n_hom1: int) -> Optional[float]:
    """FIS = 1 - H_obs / H_exp, where H_exp = 2 p (1 - p) and p is the
    allele frequency of arrangement-1. Returns None when H_exp == 0."""
    n_total = n_hom0 + n_het + n_hom1
    if n_total == 0:
        return None
    # Allele frequencies: p = (2*n_hom1 + n_het) / (2 * n_total)
    p = (2 * n_hom1 + n_het) / (2 * n_total)
    H_obs = n_het / n_total
    H_exp = 2 * p * (1 - p)
    if H_exp == 0:
        return None
    return 1.0 - H_obs / H_exp


def compute_het_excess(n_hom0: int, n_het: int, n_hom1: int) -> Optional[float]:
    """H_obs / H_exp. Values > 1 → het excess; < 1 → het deficit."""
    n_total = n_hom0 + n_het + n_hom1
    if n_total == 0:
        return None
    p = (2 * n_hom1 + n_het) / (2 * n_total)
    H_obs = n_het / n_total
    H_exp = 2 * p * (1 - p)
    if H_exp == 0:
        return None
    return H_obs / H_exp


def compute_hom_het_ratios(
    n_hom0: int, n_het: int, n_hom1: int,
) -> Tuple[float, float, float, Optional[str]]:
    """Haldane-corrected (+0.5) hom/het ratios. Returns:
      hom0_het_ratio, hom1_het_ratio, min_ratio, missing_class
    """
    r0 = (n_hom0 + 0.5) / (n_het + 0.5)
    r1 = (n_hom1 + 0.5) / (n_het + 0.5)
    r_min = min(r0, r1)
    missing = None
    if n_hom0 == 0 and n_hom1 > 0:
        missing = "HOM_REF"
    elif n_hom1 == 0 and n_hom0 > 0:
        missing = "HOM_INV"
    elif n_hom0 == 0 and n_hom1 == 0:
        missing = "BOTH"
    return r0, r1, r_min, missing


# ----------------------------------------------------------------------
# Between-arrangement dXY from biallelic marker genotypes.
# ----------------------------------------------------------------------


def _gt_to_dosage(gt: str) -> Optional[int]:
    return {"0/0": 0, "0/1": 1, "1/1": 2}.get(gt)


def compute_dxy_between_arrangements(
    *,
    arr0_samples: Sequence[str],
    arr1_samples: Sequence[str],
    marker_ids: Sequence[str],
    genotype_matrix: Dict[str, Dict[str, str]],
) -> Tuple[Optional[float], int]:
    """Between-group dXY across biallelic markers.

    For each marker, compute allele frequencies p_0 (in arr0 samples) and
    p_1 (in arr1 samples). Per-marker dXY contribution:

        d_marker = p_0 * (1 - p_1) + p_1 * (1 - p_0)

    Average across informative markers (those where at least 1 sample
    per group has a called genotype).

    Returns (dxy, n_markers_used).
    """
    if not arr0_samples or not arr1_samples:
        return (None, 0)
    contributions: List[float] = []
    for marker_id in marker_ids:
        gts = genotype_matrix.get(marker_id, {})
        # arr0 allele frequency
        d0 = [_gt_to_dosage(gts.get(s, "./.")) for s in arr0_samples]
        d0 = [x for x in d0 if x is not None]
        d1 = [_gt_to_dosage(gts.get(s, "./.")) for s in arr1_samples]
        d1 = [x for x in d1 if x is not None]
        if not d0 or not d1:
            continue
        p0 = sum(d0) / (2 * len(d0))
        p1 = sum(d1) / (2 * len(d1))
        contributions.append(p0 * (1 - p1) + p1 * (1 - p0))
    if not contributions:
        return (None, 0)
    return (sum(contributions) / len(contributions), len(contributions))


# ----------------------------------------------------------------------
# Per-LRR metric driver.
# ----------------------------------------------------------------------


def compute_lrr_metrics(
    *,
    lrr,                                                  # LRRInterval-like
    arrangement_by_sample: Dict[str, str],                # "HOM_REF"/"HET"/"HOM_INV"
    marker_ids_in_lrr: Sequence[str],
    genotype_matrix: Dict[str, Dict[str, str]],
) -> LRRComparativeMetrics:
    n_hom0, n_het, n_hom1 = genotype_class_counts(arrangement_by_sample)
    n_total = n_hom0 + n_het + n_hom1
    freq_hom0 = (n_hom0 / n_total) if n_total else None
    freq_het = (n_het / n_total) if n_total else None
    freq_hom1 = (n_hom1 / n_total) if n_total else None
    p = ((2 * n_hom1 + n_het) / (2 * n_total)) if n_total else None
    H_obs = (n_het / n_total) if n_total else None
    H_exp = (2 * p * (1 - p)) if (p is not None) else None
    fis = compute_fis(n_hom0, n_het, n_hom1)
    het_exc = compute_het_excess(n_hom0, n_het, n_hom1)
    r0, r1, r_min, missing = compute_hom_het_ratios(n_hom0, n_het, n_hom1)

    arr0 = [s for s, a in arrangement_by_sample.items() if a == "HOM_REF"]
    arr1 = [s for s, a in arrangement_by_sample.items() if a == "HOM_INV"]
    dxy, n_used = compute_dxy_between_arrangements(
        arr0_samples=arr0, arr1_samples=arr1,
        marker_ids=marker_ids_in_lrr,
        genotype_matrix=genotype_matrix,
    )

    return LRRComparativeMetrics(
        lrr_id=lrr.lrr_id, chrom=lrr.chrom, start=lrr.start, end=lrr.end,
        length=lrr.end - lrr.start,
        n_samples_called=n_total,
        n_hom0=n_hom0, n_het=n_het, n_hom1=n_hom1,
        freq_hom0=freq_hom0, freq_het=freq_het, freq_hom1=freq_hom1,
        p_allele_1=p, H_obs=H_obs, H_exp=H_exp,
        FIS=fis, heterokaryotype_enrichment=het_exc,
        hom0_het_ratio=r0, hom1_het_ratio=r1,
        min_hom_het_ratio=r_min, missing_hom_class=missing,
        n_markers_for_dxy=n_used, dxy_between_arrangements=dxy,
    )


# ----------------------------------------------------------------------
# Cross-LRR regression.
# ----------------------------------------------------------------------


def _normal_two_sided_pvalue(z: float) -> float:
    """Two-sided p-value for a standard normal statistic z."""
    return math.erfc(abs(z) / math.sqrt(2))


def linear_regression(
    xs: Sequence[float], ys: Sequence[float],
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float],
           Optional[float], Optional[float]]:
    """Simple OLS: returns (slope, intercept, r, r_squared, t_stat,
    asymptotic_p_value). All None when n < 3 or zero variance."""
    n = len(xs)
    if n < 3:
        return (None, None, None, None, None, None)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return (None, None, None, None, None, None)
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    r = sxy / math.sqrt(sxx * syy)
    r2 = r ** 2
    # t = r * sqrt(n - 2) / sqrt(1 - r^2)
    denom = math.sqrt(max(1e-15, 1.0 - r2))
    t_stat = r * math.sqrt(n - 2) / denom
    p = _normal_two_sided_pvalue(t_stat)
    return (slope, intercept, r, r2, t_stat, p)


def fit_comparative_regression(
    metrics: Sequence[LRRComparativeMetrics],
    *,
    y_kind: str = "log_min_hom_het_ratio",
) -> ComparativeRegression:
    """Fit y ~ dXY across the LRR set.

    ``y_kind`` selects the response variable:
      - ``log_min_hom_het_ratio``  → log of min_hom_het_ratio (Haldane corrected)
      - ``fis``                     → FIS
      - ``heterokaryotype_enrichment`` → H_obs / H_exp
    """
    pairs: List[Tuple[float, float]] = []
    for m in metrics:
        if m.dxy_between_arrangements is None:
            continue
        if y_kind == "log_min_hom_het_ratio":
            if m.min_hom_het_ratio is None or m.min_hom_het_ratio <= 0:
                continue
            y = math.log(m.min_hom_het_ratio)
        elif y_kind == "fis":
            if m.FIS is None:
                continue
            y = m.FIS
        elif y_kind == "heterokaryotype_enrichment":
            if m.heterokaryotype_enrichment is None:
                continue
            y = m.heterokaryotype_enrichment
        else:
            raise ValueError(f"unknown y_kind {y_kind!r}")
        pairs.append((m.dxy_between_arrangements, y))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]
    slope, intercept, r, r2, t_stat, p_val = linear_regression(xs, ys)
    caveat = ""
    if len(pairs) < 30 and p_val is not None:
        caveat = ("asymptotic normal-approximation p-value with n < 30 — "
                  "treat as indicative; compute exact t-test against "
                  "Student t_{n-2} externally if needed")
    return ComparativeRegression(
        y_name=y_kind, x_name="dxy_between_arrangements",
        n_lrrs=len(pairs),
        slope=slope, intercept=intercept,
        r=r, r_squared=r2, t_stat=t_stat,
        p_value_asymptotic=p_val,
        asymptotic_caveat=caveat,
    )


# ----------------------------------------------------------------------
# Output writers.
# ----------------------------------------------------------------------


def write_lrr_summary_tsv(
    path, metrics: Sequence[LRRComparativeMetrics],
) -> None:
    cols = list(LRRComparativeMetrics.__dataclass_fields__.keys())
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    def _fmt(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.6f}"
        if isinstance(v, bool):
            return "true" if v else "false"
        return str(v)

    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in metrics:
            fh.write("\t".join(_fmt(getattr(m, c)) for c in cols) + "\n")


def write_regression_summary_tsv(
    path, regressions: Sequence[ComparativeRegression],
) -> None:
    cols = ["y_name", "x_name", "n_lrrs", "slope", "intercept",
            "r", "r_squared", "t_stat", "p_value_asymptotic",
            "asymptotic_caveat"]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for reg in regressions:
            row = [getattr(reg, c) for c in cols]
            fh.write("\t".join(
                "" if v is None else (f"{v:.6f}" if isinstance(v, float) else str(v))
                for v in row
            ) + "\n")
