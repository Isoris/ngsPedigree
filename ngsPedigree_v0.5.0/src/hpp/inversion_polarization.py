"""
Inversion polarization + transmitted-arrangement calling (ngsPedigree → ngsTracts).

Scope: layers 2-3 of the four-layer inversion-inheritance stack. This
module does NOT call recombination tracts (that is ngsTracts' job); it
emits the polarized transmissions that ngsTracts consumes via the JSON
adapter in ngstracts_io.py.

Inputs:
  - Karyotype calls: (sample_id, inversion_id, band ∈ {0, 1, 2})
    where band 1 is the heterozygous middle PC1 cluster and bands 0/2
    are the two homozygote bands. The arrangement label (REF vs INV)
    of bands 0 and 2 is set by an external polarity hint (typically
    derived from the reference assembly's arrangement); see §1.4 of
    docs/hpp/pages/06_inversion_polarization.md for why pedigree data
    alone cannot choose polarity.
  - Dyads (parent → offspring) and triads (P1 + P2 → O) from
    ngsPedigree Stage 2.

Outputs:
  - PolarizationResult: chosen polarity, dyad/triad contradiction counts
    under BOTH orientations (reported for QC; equal up to per-individual
    confidence weighting), per-transmission calls, and the Mendelian
    drive test against H0 = 0.5.

Pedigree-data symmetry
----------------------
Under a global flip of band-0 ↔ band-2 labels, the Mendelian
compatibility predicate is invariant: a (parent-band, offspring-band)
pair that is impossible under one orientation is impossible under the
flipped orientation as well. Pure genotype data therefore cannot pick a
polarity. The polarity hint anchors this; the per-orientation
contradiction counts in PolarizationResult are useful as a pedigree-
data quality measurement, not as a polarity selector.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


# ----------------------------------------------------------------------
# Genotype-level enums.
# ----------------------------------------------------------------------

# Band labels are 0 / 1 / 2 (PC1 lower / middle / upper).
# Arrangement labels are REF / INV (the polarized assignment).
ALLELE_REF = "REF"
ALLELE_INV = "INV"
ARRANGEMENT_HOM_REF = "HOM_REF"
ARRANGEMENT_HET = "HET"
ARRANGEMENT_HOM_INV = "HOM_INV"

ALLOWED_BANDS = {0, 1, 2}
ALLOWED_POLARITY_HINTS = {"band_0_is_REF", "band_0_is_INV"}


def _arrangement_for_band(band: int, polarity: str) -> str:
    """Map a PC1 band (0/1/2) to an arrangement label under a polarity."""
    if band == 1:
        return ARRANGEMENT_HET
    if polarity == "band_0_is_REF":
        return ARRANGEMENT_HOM_REF if band == 0 else ARRANGEMENT_HOM_INV
    if polarity == "band_0_is_INV":
        return ARRANGEMENT_HOM_INV if band == 0 else ARRANGEMENT_HOM_REF
    raise ValueError(f"unknown polarity: {polarity!r}")


def _alleles(arr: str) -> Tuple[str, str]:
    """Return the two alleles of an arrangement genotype."""
    if arr == ARRANGEMENT_HOM_REF:
        return (ALLELE_REF, ALLELE_REF)
    if arr == ARRANGEMENT_HOM_INV:
        return (ALLELE_INV, ALLELE_INV)
    return (ALLELE_REF, ALLELE_INV)


# ----------------------------------------------------------------------
# Records.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class KaryotypeCall:
    sample_id: str
    inversion_id: str
    band: int
    confidence: str = "high"   # high | medium | low

    def __post_init__(self):
        if self.band not in ALLOWED_BANDS:
            raise ValueError(f"band must be in {ALLOWED_BANDS}; got {self.band!r}")


@dataclass(frozen=True)
class DyadPair:
    parent_sample_id: str
    offspring_sample_id: str
    parent_sex: Optional[str] = None   # "male" | "female" | None


@dataclass(frozen=True)
class TriadTrio:
    paternal_sample_id: str
    maternal_sample_id: str
    offspring_sample_id: str


@dataclass(frozen=True)
class CompatibilityCount:
    n_tested: int
    n_compatible: int
    n_incompatible: int
    n_ambiguous: int   # band missing / sample uncalled

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class TransmissionCall:
    inversion_id: str
    parent_sample_id: str
    offspring_sample_id: str
    parent_arrangement: str          # HOM_REF | HET | HOM_INV
    offspring_arrangement: str
    transmitted_arrangement: str     # REF | INV | ambiguous
    relationship_type: str           # "dyad" | "triad"
    co_parent_arrangement: Optional[str] = None   # triad only
    parent_sex: Optional[str] = None
    informative_for_drive: bool = False  # True only for HET parents
                                          # whose transmission is unambiguous

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class PolarizationResult:
    inversion_id: str
    polarity_hint: str
    chosen_polarity: str
    n_samples_called: int
    n_samples_uncalled: int
    band_counts: Dict[int, int]
    dyad_compat: Dict[str, CompatibilityCount]   # polarity -> count
    triad_compat: Dict[str, CompatibilityCount]
    polarities_symmetric: bool   # True iff the two orientations agree
    incompatible_dyads: List[Tuple[str, str]] = field(default_factory=list)
    incompatible_triads: List[Tuple[str, str, str]] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "inversion_id": self.inversion_id,
            "polarity_hint": self.polarity_hint,
            "chosen_polarity": self.chosen_polarity,
            "n_samples_called": self.n_samples_called,
            "n_samples_uncalled": self.n_samples_uncalled,
            "band_counts": {str(k): v for k, v in self.band_counts.items()},
            "dyad_compat": {k: v.to_dict() for k, v in self.dyad_compat.items()},
            "triad_compat": {k: v.to_dict() for k, v in self.triad_compat.items()},
            "polarities_symmetric": self.polarities_symmetric,
            "incompatible_dyads": [list(t) for t in self.incompatible_dyads],
            "incompatible_triads": [list(t) for t in self.incompatible_triads],
        }


@dataclass
class DriveStats:
    n_informative_transmissions: int
    n_REF_transmitted: int
    n_INV_transmitted: int
    INV_transmission_rate: Optional[float]
    binomial_pvalue: Optional[float]
    paternal_n: int
    paternal_INV_rate: Optional[float]
    paternal_pvalue: Optional[float]
    maternal_n: int
    maternal_INV_rate: Optional[float]
    maternal_pvalue: Optional[float]

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# Mendelian compatibility tables (in arrangement-label space).
# ----------------------------------------------------------------------

_TRIAD_ALLOWED_OFFSPRING: Dict[Tuple[str, str], Set[str]] = {
    (ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_REF): {ARRANGEMENT_HOM_REF},
    (ARRANGEMENT_HOM_REF, ARRANGEMENT_HET): {ARRANGEMENT_HOM_REF, ARRANGEMENT_HET},
    (ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV): {ARRANGEMENT_HET},
    (ARRANGEMENT_HET, ARRANGEMENT_HET):
        {ARRANGEMENT_HOM_REF, ARRANGEMENT_HET, ARRANGEMENT_HOM_INV},
    (ARRANGEMENT_HET, ARRANGEMENT_HOM_INV): {ARRANGEMENT_HET, ARRANGEMENT_HOM_INV},
    (ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_INV): {ARRANGEMENT_HOM_INV},
}


def _triad_allowed(parent1_arr: str, parent2_arr: str) -> Set[str]:
    key = tuple(sorted([parent1_arr, parent2_arr],
                       key=lambda a: (a != ARRANGEMENT_HOM_REF,
                                      a != ARRANGEMENT_HET,
                                      a != ARRANGEMENT_HOM_INV)))
    return _TRIAD_ALLOWED_OFFSPRING.get(tuple(key), set())


def _dyad_compatible(parent_arr: str, offspring_arr: str) -> bool:
    """A dyad is compatible if the parent's allele set has a non-empty
    intersection with the offspring's allele set."""
    p_alleles = set(_alleles(parent_arr))
    o_alleles = set(_alleles(offspring_arr))
    # Parent transmits one allele; offspring must carry it.
    return bool(p_alleles & o_alleles)


def _triad_compatible(p1_arr: str, p2_arr: str, off_arr: str) -> bool:
    return off_arr in _triad_allowed(p1_arr, p2_arr)


# ----------------------------------------------------------------------
# Polarization scoring.
# ----------------------------------------------------------------------


def _band_map(calls: Sequence[KaryotypeCall]) -> Dict[str, int]:
    """sample_id -> band, restricted to a single inversion."""
    return {c.sample_id: c.band for c in calls}


def _score_polarity(
    *,
    polarity: str,
    band_by_sample: Dict[str, int],
    dyads: Sequence[DyadPair],
    triads: Sequence[TriadTrio],
) -> Tuple[CompatibilityCount, CompatibilityCount,
           List[Tuple[str, str]], List[Tuple[str, str, str]]]:
    # Dyads.
    d_compat = d_incompat = d_amb = 0
    incompat_d: List[Tuple[str, str]] = []
    for dy in dyads:
        pb = band_by_sample.get(dy.parent_sample_id)
        ob = band_by_sample.get(dy.offspring_sample_id)
        if pb is None or ob is None:
            d_amb += 1
            continue
        p_arr = _arrangement_for_band(pb, polarity)
        o_arr = _arrangement_for_band(ob, polarity)
        if _dyad_compatible(p_arr, o_arr):
            d_compat += 1
        else:
            d_incompat += 1
            incompat_d.append((dy.parent_sample_id, dy.offspring_sample_id))

    # Triads.
    t_compat = t_incompat = t_amb = 0
    incompat_t: List[Tuple[str, str, str]] = []
    for tr in triads:
        pb1 = band_by_sample.get(tr.paternal_sample_id)
        pb2 = band_by_sample.get(tr.maternal_sample_id)
        ob = band_by_sample.get(tr.offspring_sample_id)
        if pb1 is None or pb2 is None or ob is None:
            t_amb += 1
            continue
        a1 = _arrangement_for_band(pb1, polarity)
        a2 = _arrangement_for_band(pb2, polarity)
        ao = _arrangement_for_band(ob, polarity)
        if _triad_compatible(a1, a2, ao):
            t_compat += 1
        else:
            t_incompat += 1
            incompat_t.append((tr.paternal_sample_id,
                               tr.maternal_sample_id,
                               tr.offspring_sample_id))

    n_d = d_compat + d_incompat + d_amb
    n_t = t_compat + t_incompat + t_amb
    return (
        CompatibilityCount(n_tested=n_d, n_compatible=d_compat,
                           n_incompatible=d_incompat, n_ambiguous=d_amb),
        CompatibilityCount(n_tested=n_t, n_compatible=t_compat,
                           n_incompatible=t_incompat, n_ambiguous=t_amb),
        incompat_d,
        incompat_t,
    )


def polarize(
    *,
    inversion_id: str,
    karyotype_calls: Sequence[KaryotypeCall],
    dyads: Sequence[DyadPair],
    triads: Sequence[TriadTrio],
    polarity_hint: str,
) -> PolarizationResult:
    """Compute compatibility counts under both orientations and choose
    the polarity from the hint.

    polarity_hint must be one of {"band_0_is_REF", "band_0_is_INV"}.
    The hint is the chosen polarity; the dual orientation is reported
    only for QC (it should be symmetric).
    """
    if polarity_hint not in ALLOWED_POLARITY_HINTS:
        raise ValueError(
            f"polarity_hint must be in {ALLOWED_POLARITY_HINTS}; got {polarity_hint!r}"
        )

    calls = [c for c in karyotype_calls if c.inversion_id == inversion_id]
    band_by = _band_map(calls)

    band_counts = {0: 0, 1: 0, 2: 0}
    for c in calls:
        band_counts[c.band] += 1

    d_a, t_a, inc_d_a, inc_t_a = _score_polarity(
        polarity="band_0_is_REF",
        band_by_sample=band_by, dyads=dyads, triads=triads,
    )
    d_b, t_b, inc_d_b, inc_t_b = _score_polarity(
        polarity="band_0_is_INV",
        band_by_sample=band_by, dyads=dyads, triads=triads,
    )

    chosen_polarity = polarity_hint
    inc_d = inc_d_a if chosen_polarity == "band_0_is_REF" else inc_d_b
    inc_t = inc_t_a if chosen_polarity == "band_0_is_REF" else inc_t_b

    # Symmetry check — the two orientations should agree exactly on
    # incompatibility counts (see module docstring).
    symmetric = (
        d_a.n_incompatible == d_b.n_incompatible
        and t_a.n_incompatible == t_b.n_incompatible
    )

    return PolarizationResult(
        inversion_id=inversion_id,
        polarity_hint=polarity_hint,
        chosen_polarity=chosen_polarity,
        n_samples_called=len(calls),
        n_samples_uncalled=0,   # callers can pass only called samples here
        band_counts=band_counts,
        dyad_compat={"band_0_is_REF": d_a, "band_0_is_INV": d_b},
        triad_compat={"band_0_is_REF": t_a, "band_0_is_INV": t_b},
        polarities_symmetric=symmetric,
        incompatible_dyads=inc_d,
        incompatible_triads=inc_t,
    )


# ----------------------------------------------------------------------
# Transmission calling.
# ----------------------------------------------------------------------


def _infer_dyad_transmission(parent_arr: str, offspring_arr: str) -> str:
    """Return REF / INV / ambiguous for the parent's transmitted allele."""
    if parent_arr == ARRANGEMENT_HOM_REF:
        return ALLELE_REF
    if parent_arr == ARRANGEMENT_HOM_INV:
        return ALLELE_INV
    # parent_arr == HET
    if offspring_arr == ARRANGEMENT_HOM_REF:
        return ALLELE_REF
    if offspring_arr == ARRANGEMENT_HOM_INV:
        return ALLELE_INV
    return "ambiguous"


def _infer_triad_transmission(
    p1_arr: str, p2_arr: str, off_arr: str
) -> Tuple[str, str]:
    """Return (p1_transmitted, p2_transmitted), each REF/INV/ambiguous."""
    if not _triad_compatible(p1_arr, p2_arr, off_arr):
        return ("contradiction", "contradiction")
    if off_arr == ARRANGEMENT_HOM_REF:
        return (ALLELE_REF, ALLELE_REF)
    if off_arr == ARRANGEMENT_HOM_INV:
        return (ALLELE_INV, ALLELE_INV)
    # offspring HET. One parent transmitted REF, the other INV.
    if p1_arr == ARRANGEMENT_HOM_REF:
        return (ALLELE_REF, ALLELE_INV)
    if p1_arr == ARRANGEMENT_HOM_INV:
        return (ALLELE_INV, ALLELE_REF)
    if p2_arr == ARRANGEMENT_HOM_REF:
        return (ALLELE_INV, ALLELE_REF)
    if p2_arr == ARRANGEMENT_HOM_INV:
        return (ALLELE_REF, ALLELE_INV)
    # Both parents HET, offspring HET → ambiguous which parent gave which.
    return ("ambiguous", "ambiguous")


def call_transmissions(
    *,
    inversion_id: str,
    karyotype_calls: Sequence[KaryotypeCall],
    dyads: Sequence[DyadPair],
    triads: Sequence[TriadTrio],
    chosen_polarity: str,
) -> List[TransmissionCall]:
    """Emit one TransmissionCall per parent→offspring relationship.

    Triads override dyads: if a (parent, offspring) pair also appears
    in a triad, the triad's transmission inference is used and the
    bare dyad is suppressed.
    """
    band_by = _band_map([c for c in karyotype_calls if c.inversion_id == inversion_id])
    sex_by = {dy.parent_sample_id: dy.parent_sex for dy in dyads if dy.parent_sex}

    out: List[TransmissionCall] = []
    seen_triad_pairs: Set[Tuple[str, str]] = set()

    for tr in triads:
        pb1 = band_by.get(tr.paternal_sample_id)
        pb2 = band_by.get(tr.maternal_sample_id)
        ob = band_by.get(tr.offspring_sample_id)
        if pb1 is None or pb2 is None or ob is None:
            continue
        a1 = _arrangement_for_band(pb1, chosen_polarity)
        a2 = _arrangement_for_band(pb2, chosen_polarity)
        ao = _arrangement_for_band(ob, chosen_polarity)
        t1, t2 = _infer_triad_transmission(a1, a2, ao)
        for parent_id, p_arr, co_arr, t_allele, sex in [
            (tr.paternal_sample_id, a1, a2, t1, "male"),
            (tr.maternal_sample_id, a2, a1, t2, "female"),
        ]:
            informative = (
                p_arr == ARRANGEMENT_HET
                and t_allele in (ALLELE_REF, ALLELE_INV)
            )
            out.append(TransmissionCall(
                inversion_id=inversion_id,
                parent_sample_id=parent_id,
                offspring_sample_id=tr.offspring_sample_id,
                parent_arrangement=p_arr,
                offspring_arrangement=ao,
                transmitted_arrangement=t_allele,
                relationship_type="triad",
                co_parent_arrangement=co_arr,
                parent_sex=sex,
                informative_for_drive=informative,
            ))
            seen_triad_pairs.add((parent_id, tr.offspring_sample_id))

    for dy in dyads:
        if (dy.parent_sample_id, dy.offspring_sample_id) in seen_triad_pairs:
            continue
        pb = band_by.get(dy.parent_sample_id)
        ob = band_by.get(dy.offspring_sample_id)
        if pb is None or ob is None:
            continue
        p_arr = _arrangement_for_band(pb, chosen_polarity)
        o_arr = _arrangement_for_band(ob, chosen_polarity)
        if not _dyad_compatible(p_arr, o_arr):
            t_allele = "contradiction"
            informative = False
        else:
            t_allele = _infer_dyad_transmission(p_arr, o_arr)
            informative = (
                p_arr == ARRANGEMENT_HET
                and t_allele in (ALLELE_REF, ALLELE_INV)
            )
        out.append(TransmissionCall(
            inversion_id=inversion_id,
            parent_sample_id=dy.parent_sample_id,
            offspring_sample_id=dy.offspring_sample_id,
            parent_arrangement=p_arr,
            offspring_arrangement=o_arr,
            transmitted_arrangement=t_allele,
            relationship_type="dyad",
            co_parent_arrangement=None,
            parent_sex=dy.parent_sex or sex_by.get(dy.parent_sample_id),
            informative_for_drive=informative,
        ))
    return out


# ----------------------------------------------------------------------
# Mendelian drive test — two-sided binomial against p=0.5, stdlib.
# ----------------------------------------------------------------------


def binomial_two_sided_pvalue(k: int, n: int, p: float = 0.5) -> float:
    """Exact two-sided binomial p-value against H0: prob = p.

    For p = 0.5 the distribution is symmetric and the p-value is
    2 * min(P(X <= k), P(X >= k)) clipped at 1.0. For p != 0.5 we use
    the "method-of-small-p-values" convention: sum P(X = j) over all
    j with P(X = j) <= P(X = k).
    """
    if not 0 <= k <= n:
        raise ValueError(f"k must be in [0, n]; got k={k}, n={n}")
    if n == 0:
        return 1.0
    if p == 0.5:
        # Symmetric two-sided.
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


def drive_test(transmissions: Sequence[TransmissionCall]) -> DriveStats:
    """Aggregate informative HET-parent transmissions and run the
    binomial drive test against H0: INV transmission rate = 0.5,
    overall and stratified by parent sex."""
    informative = [t for t in transmissions if t.informative_for_drive]
    n = len(informative)
    n_inv = sum(1 for t in informative if t.transmitted_arrangement == ALLELE_INV)
    n_ref = n - n_inv

    overall_rate = (n_inv / n) if n else None
    overall_p = binomial_two_sided_pvalue(n_inv, n, p=0.5) if n else None

    def _strat(sex_label: str) -> Tuple[int, Optional[float], Optional[float]]:
        sub = [t for t in informative if t.parent_sex == sex_label]
        ns = len(sub)
        if ns == 0:
            return (0, None, None)
        ki = sum(1 for t in sub if t.transmitted_arrangement == ALLELE_INV)
        return (ns, ki / ns, binomial_two_sided_pvalue(ki, ns, p=0.5))

    pat_n, pat_rate, pat_p = _strat("male")
    mat_n, mat_rate, mat_p = _strat("female")

    return DriveStats(
        n_informative_transmissions=n,
        n_REF_transmitted=n_ref,
        n_INV_transmitted=n_inv,
        INV_transmission_rate=overall_rate,
        binomial_pvalue=overall_p,
        paternal_n=pat_n,
        paternal_INV_rate=pat_rate,
        paternal_pvalue=pat_p,
        maternal_n=mat_n,
        maternal_INV_rate=mat_rate,
        maternal_pvalue=mat_p,
    )
