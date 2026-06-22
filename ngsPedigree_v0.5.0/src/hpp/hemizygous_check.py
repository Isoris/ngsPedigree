"""
Hemizygous deletion markers — Mendelian transmission scoring for triads
(and dyads, partial).

Why a separate module from inversion_polarization
-------------------------------------------------
The math is the same biallelic Mendelian enumeration, but the
biological interpretation and the genotype-quality caveats are
different enough to keep separate:

  - DEL markers are structural haplotypes, not PC1-cluster bands.
    Genotypes are written as ``0/0`` (no DEL), ``0/1`` (heterozygous /
    hemizygous), ``1/1`` (homozygous DEL); ``./.`` for missing.
  - ``1/1`` is treated with care in short-read data: a homozygous
    DEL call can be confounded with mapping dropout / regional
    coverage loss. ``strict_hom_del`` (default ``True``) trusts
    homozygous DEL calls as-is; setting it ``False`` treats ``1/1``
    as ambiguous and refuses to flag (1/1 parent, 0/0 child) as
    incompatible.
  - SNPs that sit INSIDE common DEL intervals must be excluded from
    SNP-level Mendelian checks (deletion-induced allele dropout will
    masquerade as homozygous SNP transmission). This module operates
    on the DEL markers themselves, not on the SNPs inside them.

Use case — "fake trio" direction inference
------------------------------------------
First-degree θ identifies candidate dyads but does not orient triads
(P1+P2 → C is symmetric with P1+C → P2 and P2+C → P1 at the kinship
level). DEL Mendelian error rates separate the three permutations:
the true child has the lowest DEL error rate against the inferred
parents. ``score_all_three_directions`` returns the per-direction
error counts; ``best_direction`` picks the direction with the lowest
error rate, gated by a minimum-informative-markers floor.

Per-direction scoring is over the set of markers where all three
samples have called genotypes. A marker is informative when the trio
genotype combination has determinable Mendelian outcome (any
combination is allowed under the Punnett table — we count
incompatibilities, not informativeness in the polymorphism sense).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


# Allowed genotype string set. Phased separators are normalised on load.
ALLOWED_GT = {"0/0", "0/1", "1/1"}
MISSING_GT = {"./.", ".", ""}

# Mapping (parent1_gt, parent2_gt) -> set of allowed child genotypes,
# using the Punnett enumeration on biallelic markers.
_TRIAD_ALLOWED_DEL_CHILD: Dict[Tuple[str, str], Set[str]] = {
    ("0/0", "0/0"): {"0/0"},
    ("0/0", "0/1"): {"0/0", "0/1"},
    ("0/1", "0/0"): {"0/0", "0/1"},
    ("0/0", "1/1"): {"0/1"},
    ("1/1", "0/0"): {"0/1"},
    ("0/1", "0/1"): {"0/0", "0/1", "1/1"},
    ("0/1", "1/1"): {"0/1", "1/1"},
    ("1/1", "0/1"): {"0/1", "1/1"},
    ("1/1", "1/1"): {"1/1"},
}


# ----------------------------------------------------------------------
# Records.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class DelMarker:
    marker_id: str
    chrom: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    notes: str = ""


@dataclass(frozen=True)
class DelCall:
    marker_id: str
    sample_id: str
    genotype: str           # "0/0" | "0/1" | "1/1" | "./."
    depth_ratio: Optional[float] = None
    confidence: str = "high"


@dataclass
class TriadDelScore:
    direction: str          # "P1+P2->C"
    p1_sample_id: str
    p2_sample_id: str
    child_sample_id: str
    n_markers_total: int
    n_informative: int      # all three GTs called and non-missing
    n_compatible: int
    n_incompatible: int
    error_rate: Optional[float]   # None when n_informative == 0
    incompatible_markers: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class DyadDelScore:
    parent_sample_id: str
    offspring_sample_id: str
    n_markers_total: int
    n_informative: int
    n_strong_incompatible: int    # (0/0, 1/1) or (1/1, 0/0) under strict mode
    n_indeterminate: int           # cases the lone-dyad cannot resolve
    incompatible_markers: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class TriadVerdict:
    a_sample_id: str
    b_sample_id: str
    c_sample_id: str
    scores: List[TriadDelScore]    # exactly three: A+B->C, A+C->B, B+C->A
    best_direction: Optional[str]  # the lowest-error direction or None
    best_error_rate: Optional[float]
    second_best_error_rate: Optional[float]
    margin: Optional[float]        # second_best - best (gap)
    informative_marker_floor: int

    def to_dict(self) -> dict:
        return {
            "a_sample_id": self.a_sample_id,
            "b_sample_id": self.b_sample_id,
            "c_sample_id": self.c_sample_id,
            "scores": [s.to_dict() for s in self.scores],
            "best_direction": self.best_direction,
            "best_error_rate": self.best_error_rate,
            "second_best_error_rate": self.second_best_error_rate,
            "margin": self.margin,
            "informative_marker_floor": self.informative_marker_floor,
        }


# ----------------------------------------------------------------------
# Genotype helpers.
# ----------------------------------------------------------------------


def normalise_gt(raw: Optional[str]) -> str:
    if raw is None:
        return "./."
    g = raw.strip()
    if not g or g in MISSING_GT:
        return "./."
    g = g.replace("|", "/")
    return g


def is_called(gt: str) -> bool:
    return gt in ALLOWED_GT


# ----------------------------------------------------------------------
# Triad Mendelian-compatibility predicate.
# ----------------------------------------------------------------------


def is_del_mendelian_compatible(
    p1_gt: str, p2_gt: str, child_gt: str,
    *,
    strict_hom_del: bool = True,
) -> Optional[bool]:
    """Return True/False for the (P1, P2 → child) Mendelian check, or
    None when any sample is missing.

    With ``strict_hom_del=False``, treat ``1/1`` as ambiguous: any
    (parent or child) ``1/1`` call is *not* used to flag incompatibility
    (it just makes the marker non-informative). This matches the
    short-read caveat in the module docstring.
    """
    p1 = normalise_gt(p1_gt)
    p2 = normalise_gt(p2_gt)
    c = normalise_gt(child_gt)
    if not (is_called(p1) and is_called(p2) and is_called(c)):
        return None
    if not strict_hom_del and ("1/1" in (p1, p2, c)):
        return None
    return c in _TRIAD_ALLOWED_DEL_CHILD[(p1, p2)]


def dyad_del_strong_incompatible(
    parent_gt: str, offspring_gt: str,
    *,
    strict_hom_del: bool = True,
) -> Optional[bool]:
    """Return True if the (parent, offspring) DEL pair is *strongly*
    incompatible (parent and offspring are opposite homozygotes), False
    if compatible, None if not testable.

    The lone-dyad test catches only (0/0, 1/1) and (1/1, 0/0); a
    heterozygous parent leaves both outcomes possible without a
    co-parent, so most dyad markers are not strongly resolvable.
    """
    p = normalise_gt(parent_gt)
    o = normalise_gt(offspring_gt)
    if not (is_called(p) and is_called(o)):
        return None
    if not strict_hom_del and "1/1" in (p, o):
        return None
    if (p, o) == ("0/0", "1/1") or (p, o) == ("1/1", "0/0"):
        return True
    return False


# ----------------------------------------------------------------------
# Per-direction triad scoring.
# ----------------------------------------------------------------------


def score_triad_direction(
    *,
    p1: str,
    p2: str,
    child: str,
    del_calls_by_marker: Dict[str, Dict[str, str]],
    strict_hom_del: bool = True,
) -> TriadDelScore:
    """Score one assumed direction (p1 + p2 → child) across all markers.

    ``del_calls_by_marker`` is ``{marker_id: {sample_id: genotype_str}}``.
    """
    n_total = len(del_calls_by_marker)
    n_inf = n_compat = n_incompat = 0
    incompat: List[str] = []
    for marker_id, calls in del_calls_by_marker.items():
        verdict = is_del_mendelian_compatible(
            calls.get(p1, "./."), calls.get(p2, "./."), calls.get(child, "./."),
            strict_hom_del=strict_hom_del,
        )
        if verdict is None:
            continue
        n_inf += 1
        if verdict:
            n_compat += 1
        else:
            n_incompat += 1
            incompat.append(marker_id)
    rate = (n_incompat / n_inf) if n_inf else None
    return TriadDelScore(
        direction=f"{p1}+{p2}->{child}",
        p1_sample_id=p1,
        p2_sample_id=p2,
        child_sample_id=child,
        n_markers_total=n_total,
        n_informative=n_inf,
        n_compatible=n_compat,
        n_incompatible=n_incompat,
        error_rate=rate,
        incompatible_markers=incompat,
    )


def score_all_three_directions(
    *,
    a: str, b: str, c: str,
    del_calls_by_marker: Dict[str, Dict[str, str]],
    strict_hom_del: bool = True,
) -> List[TriadDelScore]:
    return [
        score_triad_direction(p1=a, p2=b, child=c,
                              del_calls_by_marker=del_calls_by_marker,
                              strict_hom_del=strict_hom_del),
        score_triad_direction(p1=a, p2=c, child=b,
                              del_calls_by_marker=del_calls_by_marker,
                              strict_hom_del=strict_hom_del),
        score_triad_direction(p1=b, p2=c, child=a,
                              del_calls_by_marker=del_calls_by_marker,
                              strict_hom_del=strict_hom_del),
    ]


def best_direction(
    *,
    a: str, b: str, c: str,
    del_calls_by_marker: Dict[str, Dict[str, str]],
    strict_hom_del: bool = True,
    informative_marker_floor: int = 5,
    min_margin: float = 0.05,
) -> TriadVerdict:
    """Pick the direction with the lowest DEL Mendelian error rate.

    Returns a TriadVerdict where ``best_direction`` is set only when:
      - at least one direction has ``n_informative >= informative_marker_floor``
      - the best direction's error_rate is at least ``min_margin``
        lower than the second-best (so we do not call a winner from noise).
    """
    scores = score_all_three_directions(
        a=a, b=b, c=c,
        del_calls_by_marker=del_calls_by_marker,
        strict_hom_del=strict_hom_del,
    )
    rated = [s for s in scores
             if s.error_rate is not None and s.n_informative >= informative_marker_floor]
    if not rated:
        return TriadVerdict(
            a_sample_id=a, b_sample_id=b, c_sample_id=c,
            scores=scores,
            best_direction=None, best_error_rate=None,
            second_best_error_rate=None, margin=None,
            informative_marker_floor=informative_marker_floor,
        )
    rated.sort(key=lambda s: s.error_rate)
    best = rated[0]
    second_rate = rated[1].error_rate if len(rated) > 1 else None
    margin = (second_rate - best.error_rate) if second_rate is not None else None
    chosen: Optional[str] = best.direction
    if margin is not None and margin < min_margin:
        chosen = None
    return TriadVerdict(
        a_sample_id=a, b_sample_id=b, c_sample_id=c,
        scores=scores,
        best_direction=chosen,
        best_error_rate=best.error_rate,
        second_best_error_rate=second_rate,
        margin=margin,
        informative_marker_floor=informative_marker_floor,
    )


# ----------------------------------------------------------------------
# Dyad scoring (partial; mostly a QC count).
# ----------------------------------------------------------------------


def score_dyad(
    *,
    parent: str,
    offspring: str,
    del_calls_by_marker: Dict[str, Dict[str, str]],
    strict_hom_del: bool = True,
) -> DyadDelScore:
    n_total = len(del_calls_by_marker)
    n_inf = n_strong = n_ind = 0
    incompat: List[str] = []
    for marker_id, calls in del_calls_by_marker.items():
        v = dyad_del_strong_incompatible(
            calls.get(parent, "./."), calls.get(offspring, "./."),
            strict_hom_del=strict_hom_del,
        )
        if v is None:
            continue
        n_inf += 1
        if v:
            n_strong += 1
            incompat.append(marker_id)
        else:
            n_ind += 1
    return DyadDelScore(
        parent_sample_id=parent,
        offspring_sample_id=offspring,
        n_markers_total=n_total,
        n_informative=n_inf,
        n_strong_incompatible=n_strong,
        n_indeterminate=n_ind,
        incompatible_markers=incompat,
    )


# ----------------------------------------------------------------------
# JSON adapter — del_markers.in JSON  (one long-table file per cohort).
# ----------------------------------------------------------------------

DEL_MARKERS_SCHEMA = "ngspedigree_del_markers_v1"


class DelMarkersAdapterError(ValueError):
    pass


def load_del_markers(path) -> Dict[str, Dict[str, str]]:
    """Load a ``ngspedigree_del_markers_v1`` JSON file and return the
    nested map ``{marker_id: {sample_id: genotype_str}}`` consumed by the
    triad/dyad scorers."""
    import json
    from pathlib import Path
    p = Path(path)
    with open(p) as fh:
        doc = json.load(fh)
    if doc.get("schema") != DEL_MARKERS_SCHEMA:
        raise DelMarkersAdapterError(
            f"{p}: schema must be {DEL_MARKERS_SCHEMA!r}; got {doc.get('schema')!r}"
        )
    rows = doc.get("rows")
    if rows is None:
        raise DelMarkersAdapterError(f"{p}: missing 'rows' field")
    out: Dict[str, Dict[str, str]] = {}
    for i, r in enumerate(rows):
        for req in ("marker_id", "sample_id", "genotype"):
            if req not in r:
                raise DelMarkersAdapterError(
                    f"{p}#/rows/{i}: missing required field {req!r}"
                )
        gt = normalise_gt(r["genotype"])
        if gt != "./." and gt not in ALLOWED_GT:
            raise DelMarkersAdapterError(
                f"{p}#/rows/{i}: genotype must be in "
                f"{sorted(ALLOWED_GT | MISSING_GT)}; got {r['genotype']!r}"
            )
        out.setdefault(r["marker_id"], {})[r["sample_id"]] = gt
    return out
