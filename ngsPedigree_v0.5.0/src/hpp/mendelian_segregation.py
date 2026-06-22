"""
Per-family × per-LRR Mendelian-segregation analysis (framework steps 7–8).

For each confirmed family (a triad with one or more offspring sharing
the same parental pair) and each candidate LRR / inversion-like region,
classify the parental cross type, count observed offspring arrangement
genotypes, and test goodness-of-fit against the Mendelian expectation.

Cross types and informativeness
-------------------------------

Parental arrangements (HOM_REF / HET / HOM_INV) collapse to seven
unordered cross types. Each carries an informativeness tag:

  Cross              Expected offspring          Informative for
  HOM_REF × HOM_REF  100% HOM_REF                 — (uninformative)
  HOM_REF × HET      50% HOM_REF / 50% HET        segregation
  HOM_REF × HOM_INV  100% HET                     inheritance validation only
  HET × HET          25% HOM_REF / 50% HET / 25% HOM_INV   distortion / overdominance
  HET × HOM_INV      50% HET / 50% HOM_INV        segregation
  HOM_INV × HOM_INV  100% HOM_INV                 — (uninformative)

The most informative cross for testing pseudo-overdominance is HET × HET,
because both homozygotes are expected and homozygote depletion is the
biological signature.

Tests
-----

  - 2-class expected (50/50)        → exact two-sided binomial against p=0.5
  - 3-class expected (25/50/25)     → chi-square goodness-of-fit with df=2,
                                        analytic p-value: P(chi2 >= x) = exp(-x/2)
  - fixed (100% of one class)        → 0 expected variance; we just check
                                        whether all observed match the expected

Interpretation categories per (family, LRR):

  ``fixed_inheritance_validation``  HOM_REF × HOM_INV → all HET observed
  ``fixed_consistent``               other fixed crosses, all observed match
  ``fixed_violation``                fixed cross with any off-class offspring
                                       (pedigree error candidate)
  ``segregation_consistent``         informative cross, p ≥ 0.05
  ``segregation_distorted``          informative cross, p < 0.05
                                       and minor allele depleted
  ``homozygote_depletion``           HET × HET specifically, both homs
                                       under-represented (pseudo-overdominance candidate)
  ``small_n``                        n_offspring < min_offspring_for_test
  ``ambiguous``                      one or both parental genotypes
                                       are ambiguous or uncalled
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# Arrangement genotypes — match the polarization bloc's labels.
HOM_REF = "HOM_REF"
HET = "HET"
HOM_INV = "HOM_INV"

# Mapping band → arrangement under band_0_is_REF polarity.
BAND_TO_ARR_REF = {0: HOM_REF, 1: HET, 2: HOM_INV}


@dataclass(frozen=True)
class FamilyLrrRecord:
    family_id: str
    lrr_id: str
    paternal_sample_id: str
    maternal_sample_id: str
    paternal_arrangement: str
    maternal_arrangement: str
    n_offspring: int
    obs_HOM_REF: int
    obs_HET: int
    obs_HOM_INV: int
    cross_type: str
    expected_dist: Dict[str, float]
    test_used: str                    # "binomial" | "chi2" | "fixed_check" | "none"
    test_stat: Optional[float]
    p_value: Optional[float]
    interpretation: str
    notes: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


# ----------------------------------------------------------------------
# Cross-type classification.
# ----------------------------------------------------------------------


_CROSS_EXPECTED: Dict[Tuple[str, str], Dict[str, float]] = {
    (HOM_REF, HOM_REF): {HOM_REF: 1.0},
    (HOM_REF, HET):     {HOM_REF: 0.5, HET: 0.5},
    (HOM_REF, HOM_INV): {HET: 1.0},
    (HET, HET):         {HOM_REF: 0.25, HET: 0.5, HOM_INV: 0.25},
    (HET, HOM_INV):     {HET: 0.5, HOM_INV: 0.5},
    (HOM_INV, HOM_INV): {HOM_INV: 1.0},
}


def _canonical_cross(p1: str, p2: str) -> Tuple[str, str]:
    """Order-independent cross key (so (HOM_REF, HET) == (HET, HOM_REF))."""
    order = {HOM_REF: 0, HET: 1, HOM_INV: 2}
    if order[p1] <= order[p2]:
        return (p1, p2)
    return (p2, p1)


def cross_type(p1: str, p2: str) -> str:
    k = _canonical_cross(p1, p2)
    return f"{k[0]}_x_{k[1]}"


def expected_distribution(p1: str, p2: str) -> Dict[str, float]:
    return dict(_CROSS_EXPECTED[_canonical_cross(p1, p2)])


def is_informative_cross(p1: str, p2: str) -> bool:
    """True iff the cross produces more than one offspring class."""
    return len(expected_distribution(p1, p2)) > 1


# ----------------------------------------------------------------------
# Statistical tests (stdlib).
# ----------------------------------------------------------------------


def binomial_two_sided_pvalue(k: int, n: int, p: float = 0.5) -> float:
    if n == 0:
        return 1.0
    if not 0 <= k <= n:
        raise ValueError(f"k must be in [0, n]; got k={k}, n={n}")
    if p == 0.5:
        m = min(k, n - k)
        tail = sum(math.comb(n, j) for j in range(0, m + 1)) * (0.5 ** n)
        return min(1.0, 2.0 * tail)
    p_obs = math.comb(n, k) * (p ** k) * ((1 - p) ** (n - k))
    total = 0.0
    for j in range(n + 1):
        p_j = math.comb(n, j) * (p ** j) * ((1 - p) ** (n - j))
        if p_j <= p_obs + 1e-15:
            total += p_j
    return min(1.0, total)


def chi2_gof_df2(observed: Sequence[float], expected: Sequence[float]) -> Tuple[float, float]:
    """Chi-square goodness-of-fit for df=2 (i.e. 3 classes). Returns
    (chi2_stat, p_value). p computed analytically: P(X >= x) = exp(-x/2)."""
    chi2 = 0.0
    for o, e in zip(observed, expected):
        if e == 0:
            if o != 0:
                return (float("inf"), 0.0)
            continue
        chi2 += (o - e) ** 2 / e
    p_val = math.exp(-chi2 / 2.0) if chi2 >= 0 else 1.0
    return chi2, p_val


# ----------------------------------------------------------------------
# Family × LRR analysis.
# ----------------------------------------------------------------------


def analyze_family_lrr(
    *,
    family_id: str,
    lrr_id: str,
    paternal_sample_id: str,
    maternal_sample_id: str,
    paternal_arrangement: str,
    maternal_arrangement: str,
    offspring_arrangements: Sequence[str],
    min_offspring_for_test: int = 4,
    distortion_alpha: float = 0.05,
) -> FamilyLrrRecord:
    n_off = len(offspring_arrangements)
    obs = {HOM_REF: 0, HET: 0, HOM_INV: 0}
    for a in offspring_arrangements:
        if a in obs:
            obs[a] += 1

    # Handle ambiguous parents up-front.
    if (paternal_arrangement not in (HOM_REF, HET, HOM_INV)
            or maternal_arrangement not in (HOM_REF, HET, HOM_INV)):
        return FamilyLrrRecord(
            family_id=family_id, lrr_id=lrr_id,
            paternal_sample_id=paternal_sample_id,
            maternal_sample_id=maternal_sample_id,
            paternal_arrangement=paternal_arrangement,
            maternal_arrangement=maternal_arrangement,
            n_offspring=n_off,
            obs_HOM_REF=obs[HOM_REF], obs_HET=obs[HET], obs_HOM_INV=obs[HOM_INV],
            cross_type="unresolved",
            expected_dist={},
            test_used="none",
            test_stat=None, p_value=None,
            interpretation="ambiguous",
            notes="parental arrangement uncalled or ambiguous",
        )

    ctype = cross_type(paternal_arrangement, maternal_arrangement)
    expected = expected_distribution(paternal_arrangement, maternal_arrangement)
    expected_counts = [expected.get(k, 0.0) * n_off for k in (HOM_REF, HET, HOM_INV)]
    observed_counts = [obs[HOM_REF], obs[HET], obs[HOM_INV]]

    # Fixed (single-class) cross.
    if len(expected) == 1:
        target = next(iter(expected))
        all_match = (obs[target] == n_off)
        if n_off == 0:
            interpretation = "small_n"
            test_used = "none"
            stat = p_val = None
        elif all_match:
            interpretation = (
                "fixed_inheritance_validation"
                if ctype == f"{HOM_REF}_x_{HOM_INV}" else "fixed_consistent"
            )
            test_used = "fixed_check"
            stat = 0.0
            p_val = 1.0
        else:
            interpretation = "fixed_violation"
            test_used = "fixed_check"
            stat = float(n_off - obs[target])
            p_val = 0.0
        return FamilyLrrRecord(
            family_id=family_id, lrr_id=lrr_id,
            paternal_sample_id=paternal_sample_id,
            maternal_sample_id=maternal_sample_id,
            paternal_arrangement=paternal_arrangement,
            maternal_arrangement=maternal_arrangement,
            n_offspring=n_off,
            obs_HOM_REF=obs[HOM_REF], obs_HET=obs[HET], obs_HOM_INV=obs[HOM_INV],
            cross_type=ctype,
            expected_dist=expected,
            test_used=test_used,
            test_stat=stat, p_value=p_val,
            interpretation=interpretation,
        )

    if n_off < min_offspring_for_test:
        return FamilyLrrRecord(
            family_id=family_id, lrr_id=lrr_id,
            paternal_sample_id=paternal_sample_id,
            maternal_sample_id=maternal_sample_id,
            paternal_arrangement=paternal_arrangement,
            maternal_arrangement=maternal_arrangement,
            n_offspring=n_off,
            obs_HOM_REF=obs[HOM_REF], obs_HET=obs[HET], obs_HOM_INV=obs[HOM_INV],
            cross_type=ctype,
            expected_dist=expected,
            test_used="none",
            test_stat=None, p_value=None,
            interpretation="small_n",
            notes=f"n_offspring={n_off} < min_offspring_for_test={min_offspring_for_test}",
        )

    # 2-class informative cross → binomial.
    if len(expected) == 2:
        keys = sorted(expected, key=lambda k: -expected[k])
        # In a 2-class cross with equal expectations, p=0.5 for either class.
        # Use the "first" class count for the binomial.
        k_count = obs[keys[0]]
        p_val = binomial_two_sided_pvalue(k_count, n_off, p=0.5)
        stat = float(k_count)
        if p_val < distortion_alpha:
            interpretation = "segregation_distorted"
        else:
            interpretation = "segregation_consistent"
        return FamilyLrrRecord(
            family_id=family_id, lrr_id=lrr_id,
            paternal_sample_id=paternal_sample_id,
            maternal_sample_id=maternal_sample_id,
            paternal_arrangement=paternal_arrangement,
            maternal_arrangement=maternal_arrangement,
            n_offspring=n_off,
            obs_HOM_REF=obs[HOM_REF], obs_HET=obs[HET], obs_HOM_INV=obs[HOM_INV],
            cross_type=ctype,
            expected_dist=expected,
            test_used="binomial",
            test_stat=stat, p_value=p_val,
            interpretation=interpretation,
        )

    # 3-class HET x HET → chi-square df=2.
    chi2, p_val = chi2_gof_df2(observed_counts, expected_counts)
    # Pseudo-overdominance signature: both homozygotes under-represented.
    interpretation = "segregation_consistent"
    if p_val < distortion_alpha:
        hom_obs = obs[HOM_REF] + obs[HOM_INV]
        hom_exp = expected_counts[0] + expected_counts[2]
        if hom_obs < hom_exp:
            interpretation = "homozygote_depletion"
        else:
            interpretation = "segregation_distorted"
    return FamilyLrrRecord(
        family_id=family_id, lrr_id=lrr_id,
        paternal_sample_id=paternal_sample_id,
        maternal_sample_id=maternal_sample_id,
        paternal_arrangement=paternal_arrangement,
        maternal_arrangement=maternal_arrangement,
        n_offspring=n_off,
        obs_HOM_REF=obs[HOM_REF], obs_HET=obs[HET], obs_HOM_INV=obs[HOM_INV],
        cross_type=ctype,
        expected_dist=expected,
        test_used="chi2",
        test_stat=chi2, p_value=p_val,
        interpretation=interpretation,
    )


# ----------------------------------------------------------------------
# Cohort-level aggregation across families for one LRR.
# ----------------------------------------------------------------------


@dataclass
class LrrSegregationSummary:
    lrr_id: str
    n_families: int
    n_informative_families: int
    n_segregation_consistent: int
    n_segregation_distorted: int
    n_homozygote_depletion: int
    n_fixed_inheritance_validation: int
    n_fixed_consistent: int
    n_fixed_violation: int
    n_small_n: int
    n_ambiguous: int
    # Cohort-wide HET × HET aggregate counts.
    het_het_obs_HOM_REF: int = 0
    het_het_obs_HET: int = 0
    het_het_obs_HOM_INV: int = 0
    het_het_n_total: int = 0
    het_het_chi2: Optional[float] = None
    het_het_p_value: Optional[float] = None
    het_het_interpretation: str = "none"

    def to_dict(self) -> dict:
        return asdict(self)


def summarise_lrr(
    lrr_id: str,
    records: Sequence[FamilyLrrRecord],
    *,
    distortion_alpha: float = 0.05,
) -> LrrSegregationSummary:
    s = LrrSegregationSummary(
        lrr_id=lrr_id, n_families=len(records),
        n_informative_families=0,
        n_segregation_consistent=0, n_segregation_distorted=0,
        n_homozygote_depletion=0, n_fixed_inheritance_validation=0,
        n_fixed_consistent=0, n_fixed_violation=0,
        n_small_n=0, n_ambiguous=0,
    )
    for r in records:
        if r.interpretation == "segregation_consistent":
            s.n_segregation_consistent += 1; s.n_informative_families += 1
        elif r.interpretation == "segregation_distorted":
            s.n_segregation_distorted += 1; s.n_informative_families += 1
        elif r.interpretation == "homozygote_depletion":
            s.n_homozygote_depletion += 1; s.n_informative_families += 1
        elif r.interpretation == "fixed_inheritance_validation":
            s.n_fixed_inheritance_validation += 1
        elif r.interpretation == "fixed_consistent":
            s.n_fixed_consistent += 1
        elif r.interpretation == "fixed_violation":
            s.n_fixed_violation += 1
        elif r.interpretation == "small_n":
            s.n_small_n += 1
        elif r.interpretation == "ambiguous":
            s.n_ambiguous += 1
        # Aggregate HET × HET offspring counts.
        if r.cross_type == f"{HET}_x_{HET}":
            s.het_het_obs_HOM_REF += r.obs_HOM_REF
            s.het_het_obs_HET += r.obs_HET
            s.het_het_obs_HOM_INV += r.obs_HOM_INV
            s.het_het_n_total += r.n_offspring
    if s.het_het_n_total > 0:
        expected = [
            0.25 * s.het_het_n_total, 0.5 * s.het_het_n_total,
            0.25 * s.het_het_n_total,
        ]
        observed = [s.het_het_obs_HOM_REF, s.het_het_obs_HET, s.het_het_obs_HOM_INV]
        chi2, p_val = chi2_gof_df2(observed, expected)
        s.het_het_chi2 = chi2
        s.het_het_p_value = p_val
        if p_val < distortion_alpha:
            hom_obs = s.het_het_obs_HOM_REF + s.het_het_obs_HOM_INV
            hom_exp = 0.5 * s.het_het_n_total
            s.het_het_interpretation = (
                "cohort_homozygote_depletion"
                if hom_obs < hom_exp else "cohort_distorted"
            )
        else:
            s.het_het_interpretation = "cohort_consistent"
    return s


# ----------------------------------------------------------------------
# TSV writers.
# ----------------------------------------------------------------------


def write_family_lrr_tsv(path, records: Sequence[FamilyLrrRecord]) -> None:
    cols = ["family_id", "lrr_id",
            "paternal_sample_id", "maternal_sample_id",
            "paternal_arrangement", "maternal_arrangement",
            "n_offspring", "obs_HOM_REF", "obs_HET", "obs_HOM_INV",
            "cross_type", "expected_HOM_REF", "expected_HET", "expected_HOM_INV",
            "test_used", "test_stat", "p_value",
            "interpretation", "notes"]
    from pathlib import Path
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            exp = r.expected_dist
            row = [
                r.family_id, r.lrr_id,
                r.paternal_sample_id, r.maternal_sample_id,
                r.paternal_arrangement, r.maternal_arrangement,
                r.n_offspring, r.obs_HOM_REF, r.obs_HET, r.obs_HOM_INV,
                r.cross_type,
                exp.get(HOM_REF, 0.0), exp.get(HET, 0.0), exp.get(HOM_INV, 0.0),
                r.test_used,
                "" if r.test_stat is None else f"{r.test_stat:.4f}",
                "" if r.p_value is None else f"{r.p_value:.6f}",
                r.interpretation, r.notes,
            ]
            fh.write("\t".join(
                f"{v:.4f}" if isinstance(v, float) else str(v) for v in row
            ) + "\n")


def write_lrr_summary_tsv(path,
                          summaries: Sequence[LrrSegregationSummary]) -> None:
    cols = ["lrr_id", "n_families", "n_informative_families",
            "n_segregation_consistent", "n_segregation_distorted",
            "n_homozygote_depletion",
            "n_fixed_inheritance_validation", "n_fixed_consistent",
            "n_fixed_violation",
            "n_small_n", "n_ambiguous",
            "het_het_n_total", "het_het_obs_HOM_REF", "het_het_obs_HET",
            "het_het_obs_HOM_INV",
            "het_het_chi2", "het_het_p_value", "het_het_interpretation"]
    from pathlib import Path
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in summaries:
            row = [getattr(s, c) for c in cols]
            fh.write("\t".join(
                f"{v:.4f}" if isinstance(v, float) else
                ("" if v is None else str(v))
                for v in row
            ) + "\n")
