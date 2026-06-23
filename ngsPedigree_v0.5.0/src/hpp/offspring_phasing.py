"""
Offspring-first arrangement phasing — reverse-segregation marker
discovery from variant dosage in sibship cohorts.

Method
------
For each (family × chromosome interval) the input is a children ×
variants dosage matrix. Offspring are split into two segregation
classes by k-means on the filtered variants, and every variant is then
scored by class separation:

    marker_score = | mean(dosage | class A) - mean(dosage | class B) |

Variants are NOT assumed to be SNPs. Any presence/absence marker is
accepted — DEL, DUP, INV, BND, indel, normalised depth, split-read
support, paired-end support, or a soft probability in [0, 1]. The
encoding rule is:

    0.00 - 0.49   → "absent"        (counted as class 0 / dosage 0.0)
    0.50          → "uncertain"     (intermediate)
    0.51 - 1.00   → "present"       (counted as class 1 / dosage 1.0)

For categorical genotype strings the conversion is::

    "0/0" → 0.0
    "0/1" → 0.5
    "1/1" → 1.0
    other → None  (treated as missing)

Filter modes
------------
The variant filter supports two modes that can be invoked side by side
so the user can compare them:

  - ``segregating`` (default): keep variants with carrier frequency in
    [min_freq, max_freq]. The classic 0.2-0.8 segregation band.
  - ``hemizygous_only``: keep variants where the majority of informative
    offspring are heterozygous (intermediate dosage near 0.5),
    consistent with the hom-depleted "situation 1" pattern in which
    both homokaryotypes are absent or extremely rare.

Both modes feed the same downstream scoring, so the two output tables
are directly comparable.

Honesty caveats
---------------
This module does NOT phase SNPs, does NOT estimate arrangement age,
does NOT call lethality, does NOT make any fitness claim. It outputs
arrangement-tagging marker candidates only. The polarity of class A
versus class B is arbitrary (band-relabel symmetry); only the within-
family marker score is meaningful.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ----------------------------------------------------------------------
# Data classes.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class VariantMarker:
    variant_id: str
    chrom: str
    pos: int
    end: int
    variant_type: str        # DEL / DUP / INV / BND / indel / SNP / depth
    notes: str = ""


@dataclass(frozen=True)
class OffspringPhasingRecord:
    family_id: str
    chrom: str
    interval_start: int
    interval_end: int
    variant_id: str
    variant_type: str
    pos: int
    n_class_a: int
    n_class_b: int
    mean_dosage_class_a: Optional[float]
    mean_dosage_class_b: Optional[float]
    marker_score: float
    n_informative: int
    missingness: float
    filter_mode: str
    segregation_pattern: str
    inferred_parental_state: str
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# Dosage encoding.
# ----------------------------------------------------------------------


def encode_dosage(value) -> Optional[float]:
    """Return a dosage in [0, 1] or None for missing.

    Accepts GT strings ("0/0", "0/1", "1/1", "./."), integer dosages
    (0, 1, 2), and floats already in [0, 1] or [0, 2]. Float inputs in
    the (1, 2] range are normalised by dividing by 2.
    """
    if value is None:
        return None
    if isinstance(value, str):
        v = value.strip()
        if v == "0/0":
            return 0.0
        if v == "0/1" or v == "1/0":
            return 0.5
        if v == "1/1":
            return 1.0
        if v in ("./.", ".", ""):
            return None
        try:
            value = float(v)
        except ValueError:
            return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f:                  # NaN
        return None
    if f < 0.0 or f > 2.0:
        return None
    if f > 1.0:
        f = f / 2.0
    return f


# ----------------------------------------------------------------------
# Variant filters.
# ----------------------------------------------------------------------


def _per_variant_dosages(
    variant_id: str,
    dosage_matrix: Dict[str, Dict[str, float]],
    samples: Sequence[str],
) -> List[float]:
    row = dosage_matrix.get(variant_id, {})
    out: List[float] = []
    for s in samples:
        d = encode_dosage(row.get(s))
        if d is not None:
            out.append(d)
    return out


def carrier_frequency(dosages: Sequence[float]) -> Optional[float]:
    """Fraction of offspring with any non-zero dosage."""
    if not dosages:
        return None
    return sum(1 for d in dosages if d > 0.0) / len(dosages)


def het_fraction(
    dosages: Sequence[float],
    low: float = 0.25,
    high: float = 0.75,
) -> Optional[float]:
    """Fraction of offspring whose dosage sits in the intermediate band."""
    if not dosages:
        return None
    return sum(1 for d in dosages if low <= d <= high) / len(dosages)


def filter_variants(
    variant_ids: Sequence[str],
    dosage_matrix: Dict[str, Dict[str, float]],
    samples: Sequence[str],
    *,
    mode: str = "segregating",
    min_freq: float = 0.20,
    max_freq: float = 0.80,
    min_het_fraction: float = 0.60,
    max_missingness: float = 0.40,
) -> List[str]:
    """Return the subset of ``variant_ids`` passing the chosen filter.

    Modes:
      - ``segregating``: carrier frequency in [min_freq, max_freq]
      - ``hemizygous_only``: het fraction >= min_het_fraction (the
        situation-1 pattern: nearly every informative offspring is
        heterozygous because both homokaryotypes are absent)
    """
    if mode not in ("segregating", "hemizygous_only"):
        raise ValueError(f"unknown filter mode: {mode!r}")
    kept: List[str] = []
    n_samples = len(samples)
    for vid in variant_ids:
        dos = _per_variant_dosages(vid, dosage_matrix, samples)
        if not dos:
            continue
        miss = 1.0 - (len(dos) / n_samples) if n_samples else 1.0
        if miss > max_missingness:
            continue
        if mode == "segregating":
            f = carrier_frequency(dos)
            if f is None:
                continue
            if min_freq <= f <= max_freq:
                kept.append(vid)
        else:  # hemizygous_only
            h = het_fraction(dos)
            if h is None:
                continue
            if h >= min_het_fraction:
                kept.append(vid)
    return kept


# ----------------------------------------------------------------------
# Distance + two-class clustering.
# ----------------------------------------------------------------------


def _sample_vector(
    sample: str,
    variant_ids: Sequence[str],
    dosage_matrix: Dict[str, Dict[str, float]],
) -> List[Optional[float]]:
    return [encode_dosage(dosage_matrix.get(v, {}).get(sample))
            for v in variant_ids]


def hamming_like_distance(
    vec_a: Sequence[Optional[float]],
    vec_b: Sequence[Optional[float]],
) -> Tuple[float, int]:
    """L1 distance over jointly observed positions. Returns
    (total_distance, n_pairs_compared). Missing → skip."""
    total = 0.0
    n = 0
    for a, b in zip(vec_a, vec_b):
        if a is None or b is None:
            continue
        total += abs(a - b)
        n += 1
    return (total, n)


def _seed_pair(
    samples: Sequence[str],
    vectors: Dict[str, Sequence[Optional[float]]],
) -> Tuple[str, str]:
    """Pick the two samples furthest apart (deterministic tiebreak by
    sample order)."""
    best_pair = (samples[0], samples[-1])
    best_d = -1.0
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            d, n = hamming_like_distance(vectors[a], vectors[b])
            if n == 0:
                continue
            avg = d / n
            if avg > best_d:
                best_d = avg
                best_pair = (a, b)
    return best_pair


def _centroid(
    members: Sequence[str],
    variant_ids: Sequence[str],
    vectors: Dict[str, Sequence[Optional[float]]],
) -> List[Optional[float]]:
    out: List[Optional[float]] = []
    for j in range(len(variant_ids)):
        vals = [vectors[m][j] for m in members if vectors[m][j] is not None]
        out.append(sum(vals) / len(vals) if vals else None)
    return out


def cluster_offspring_two_class(
    samples: Sequence[str],
    variant_ids: Sequence[str],
    dosage_matrix: Dict[str, Dict[str, float]],
    *,
    max_iter: int = 20,
) -> Dict[str, str]:
    """Two-class k-means-like split of offspring over the filtered
    variants. Returns sample_id → "A" | "B" | "U" (unassigned when no
    informative variants exist)."""
    if not variant_ids or len(samples) < 2:
        return {s: "U" for s in samples}

    vectors = {s: _sample_vector(s, variant_ids, dosage_matrix)
               for s in samples}
    # If a sample has no observed dosages at all → unassigned.
    informative = [s for s in samples
                   if any(v is not None for v in vectors[s])]
    if len(informative) < 2:
        return {s: "U" for s in samples}

    seed_a, seed_b = _seed_pair(informative, vectors)
    cent_a: Sequence[Optional[float]] = vectors[seed_a]
    cent_b: Sequence[Optional[float]] = vectors[seed_b]

    assign: Dict[str, str] = {}
    for _ in range(max_iter):
        changed = False
        new_assign: Dict[str, str] = {}
        for s in informative:
            da, na = hamming_like_distance(vectors[s], cent_a)
            db, nb = hamming_like_distance(vectors[s], cent_b)
            if na == 0 and nb == 0:
                new_assign[s] = "U"
                continue
            if na == 0:
                new_assign[s] = "B"
            elif nb == 0:
                new_assign[s] = "A"
            else:
                avg_a = da / na
                avg_b = db / nb
                new_assign[s] = "A" if avg_a <= avg_b else "B"
            if assign.get(s) != new_assign[s]:
                changed = True
        assign = new_assign
        members_a = [s for s in informative if assign[s] == "A"]
        members_b = [s for s in informative if assign[s] == "B"]
        if not members_a or not members_b:
            # Collapse to single class; nothing to refine.
            break
        cent_a = _centroid(members_a, variant_ids, vectors)
        cent_b = _centroid(members_b, variant_ids, vectors)
        if not changed:
            break

    final: Dict[str, str] = {}
    for s in samples:
        final[s] = assign.get(s, "U")
    return final


# ----------------------------------------------------------------------
# Variant scoring.
# ----------------------------------------------------------------------


def _classify_segregation_pattern(dosages: Sequence[float]) -> str:
    """Tag the observed offspring distribution by Mendelian shape."""
    if not dosages:
        return "no_data"
    n = len(dosages)
    n_absent = sum(1 for d in dosages if d <= 0.25)
    n_het = sum(1 for d in dosages if 0.25 < d < 0.75)
    n_hom = sum(1 for d in dosages if d >= 0.75)

    def _near(a, b, tol=0.15):
        return abs(a - b) <= tol * n

    if n_het / n >= 0.80:
        return "all_het"        # hemizygous-only / situation 1
    if _near(n_absent, n_het + n_hom) and n_hom == 0:
        return "1:1"            # AA x AB
    if _near(n_hom, n_het + n_absent) and n_absent == 0:
        return "1:1_inverted"   # BB x AB
    if _near(n_absent, n_hom) and (n_het / n) >= 0.40:
        return "1:2:1"          # AB x AB
    if n_hom == n and n_het == 0 and n_absent == 0:
        return "fixed_present"
    if n_absent == n and n_het == 0 and n_hom == 0:
        return "fixed_absent"
    return "intermediate"


def _infer_parental_state(
    segregation: str,
    mode: str,
) -> str:
    if mode == "hemizygous_only":
        return "het_x_absent_or_hom_depleted"
    if segregation == "1:1":
        return "AA_x_AB"
    if segregation == "1:1_inverted":
        return "BB_x_AB"
    if segregation == "1:2:1":
        return "AB_x_AB"
    if segregation == "all_het":
        return "het_x_absent_or_hom_depleted"
    if segregation == "fixed_present":
        return "BB_x_BB"
    if segregation == "fixed_absent":
        return "AA_x_AA"
    return "unknown"


def score_variant_class_separation(
    variant_id: str,
    samples: Sequence[str],
    dosage_matrix: Dict[str, Dict[str, float]],
    classes: Dict[str, str],
) -> Tuple[int, int, Optional[float], Optional[float], float, int, float]:
    """Return
    (n_class_a, n_class_b, mean_a, mean_b, marker_score,
     n_informative, missingness)."""
    row = dosage_matrix.get(variant_id, {})
    a_vals: List[float] = []
    b_vals: List[float] = []
    n_obs = 0
    for s in samples:
        d = encode_dosage(row.get(s))
        if d is None:
            continue
        n_obs += 1
        cls = classes.get(s, "U")
        if cls == "A":
            a_vals.append(d)
        elif cls == "B":
            b_vals.append(d)
    mean_a = sum(a_vals) / len(a_vals) if a_vals else None
    mean_b = sum(b_vals) / len(b_vals) if b_vals else None
    if mean_a is None or mean_b is None:
        score = 0.0
    else:
        score = abs(mean_a - mean_b)
    miss = 1.0 - (n_obs / len(samples)) if samples else 1.0
    return (len(a_vals), len(b_vals), mean_a, mean_b,
            score, n_obs, miss)


# ----------------------------------------------------------------------
# Orchestrator: one (family × chromosome × interval) call.
# ----------------------------------------------------------------------


def phase_offspring_interval(
    *,
    family_id: str,
    chrom: str,
    interval_start: int,
    interval_end: int,
    markers: Sequence[VariantMarker],
    dosage_matrix: Dict[str, Dict[str, float]],
    offspring: Sequence[str],
    mode: str = "segregating",
    min_freq: float = 0.20,
    max_freq: float = 0.80,
    min_het_fraction: float = 0.60,
    max_missingness: float = 0.40,
    min_variants_for_clustering: int = 3,
) -> Tuple[Dict[str, str], List[OffspringPhasingRecord]]:
    """Run the offspring-first phasing on a single chromosome interval.

    Returns (offspring → class_label, per-variant phasing records).
    """
    in_interval = [m for m in markers
                   if m.chrom == chrom
                   and interval_start <= m.pos < interval_end]
    variant_ids = [m.variant_id for m in in_interval]
    kept = filter_variants(
        variant_ids, dosage_matrix, offspring,
        mode=mode, min_freq=min_freq, max_freq=max_freq,
        min_het_fraction=min_het_fraction,
        max_missingness=max_missingness,
    )

    if len(kept) < min_variants_for_clustering:
        classes = {s: "U" for s in offspring}
    else:
        classes = cluster_offspring_two_class(
            offspring, kept, dosage_matrix,
        )

    by_id = {m.variant_id: m for m in in_interval}
    records: List[OffspringPhasingRecord] = []
    for vid in kept:
        m = by_id[vid]
        n_a, n_b, mean_a, mean_b, score, n_obs, miss = (
            score_variant_class_separation(
                vid, offspring, dosage_matrix, classes,
            )
        )
        dos = _per_variant_dosages(vid, dosage_matrix, offspring)
        seg = _classify_segregation_pattern(dos)
        parental = _infer_parental_state(seg, mode)
        records.append(OffspringPhasingRecord(
            family_id=family_id,
            chrom=chrom,
            interval_start=interval_start,
            interval_end=interval_end,
            variant_id=vid,
            variant_type=m.variant_type,
            pos=m.pos,
            n_class_a=n_a,
            n_class_b=n_b,
            mean_dosage_class_a=mean_a,
            mean_dosage_class_b=mean_b,
            marker_score=score,
            n_informative=n_obs,
            missingness=miss,
            filter_mode=mode,
            segregation_pattern=seg,
            inferred_parental_state=parental,
        ))

    records.sort(key=lambda r: r.marker_score, reverse=True)
    return classes, records


# ----------------------------------------------------------------------
# TSV writers.
# ----------------------------------------------------------------------


def write_phasing_tsv(
    path,
    records: Sequence[OffspringPhasingRecord],
) -> None:
    cols = [
        "family_id", "chrom", "interval_start", "interval_end",
        "variant_id", "variant_type", "pos",
        "n_class_a", "n_class_b",
        "mean_dosage_class_a", "mean_dosage_class_b",
        "marker_score", "n_informative", "missingness",
        "filter_mode", "segregation_pattern",
        "inferred_parental_state", "notes",
    ]
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
            fh.write("\t".join(_fmt(getattr(r, c)) for c in cols) + "\n")


def write_class_assignments_tsv(
    path,
    family_id: str,
    chrom: str,
    interval_start: int,
    interval_end: int,
    classes: Dict[str, str],
    filter_mode: str,
) -> None:
    cols = ["family_id", "chrom", "interval_start", "interval_end",
            "sample_id", "class_label", "filter_mode"]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sid in sorted(classes):
            fh.write("\t".join([
                family_id, chrom,
                str(interval_start), str(interval_end),
                sid, classes[sid], filter_mode,
            ]) + "\n")
