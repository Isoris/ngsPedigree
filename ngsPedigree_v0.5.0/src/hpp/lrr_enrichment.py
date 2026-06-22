"""
Family-based odds-ratio enrichment for candidate LRR (low-recombination
region) intervals against matched background.

The claim is the user-correct one: this does NOT prove the molecular
inversion breakpoint. It proves that the candidate LRR is enriched for
**transmission-compatible block inheritance** across dyads/triads —
i.e., it behaves like an inherited recombination-suppressed haplotype
block, in contrast to random background regions where transmission
switches more freely.

Per (relationship, region), classify the region as "block-compatible"
when:

  - the region contains ≥ ``min_markers`` informative markers
    (parent het, offspring resolvable);
  - the dominant transmitted allele covers ≥ ``dominance_threshold``
    of those markers (default 0.8);
  - and there are zero Mendelian contradictions within the region.

Then for each candidate LRR draw ``n_background`` matched background
windows (same chromosome, same approximate length, non-overlapping
with any LRR), classify those, and assemble a 2×2:

    block-compat ¬block-compat
    inside LRR        a            b
    outside LRR       c            d

Report OR = (a/b)/(c/d), report log-OR 95% CI via the standard
Woolf SE = sqrt(1/a + 1/b + 1/c + 1/d) with Haldane-Anscombe (+0.5)
when any cell is zero. Also stratify by relationship type so OR_triad
and OR_dyad are reported separately (triads carry more inheritance
information per region than dyads).
"""

from __future__ import annotations

import math
import random
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .del_inheritance import (
    DelMarkerLocus,
    TransmittedAllele,
    transmitted_from_parent,
    transmitted_from_parent_triad,
)


@dataclass(frozen=True)
class LRRInterval:
    lrr_id: str
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


@dataclass(frozen=True)
class RegionClassification:
    relationship_id: str
    relationship_type: str    # "dyad" | "triad"
    parent_sample_id: str
    offspring_sample_id: str
    region_id: str
    region_kind: str          # "LRR" | "background"
    chrom: str
    start: int
    end: int
    n_informative: int
    n_dominant_allele: int
    dominance: float
    mendelian_errors: int
    block_compatible: bool

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class OddsRatio:
    a: int                     # inside, block-compatible
    b: int                     # inside, not
    c: int                     # outside, block-compatible
    d: int                     # outside, not
    odds_ratio: Optional[float]
    log_or: Optional[float]
    se_log_or: Optional[float]
    ci_low: Optional[float]
    ci_high: Optional[float]
    haldane_corrected: bool

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class LRREnrichment:
    lrr_id: str
    chrom: str
    start: int
    end: int
    n_relationships: int
    n_triads: int
    n_dyads: int
    combined: OddsRatio
    triad_only: Optional[OddsRatio]
    dyad_only: Optional[OddsRatio]

    def to_dict(self) -> dict:
        d = asdict(self)
        d["combined"] = self.combined.to_dict()
        d["triad_only"] = self.triad_only.to_dict() if self.triad_only else None
        d["dyad_only"] = self.dyad_only.to_dict() if self.dyad_only else None
        return d


# ----------------------------------------------------------------------
# Region classification.
# ----------------------------------------------------------------------


def classify_region_dyad(
    *,
    parent_sample_id: str,
    offspring_sample_id: str,
    chrom: str,
    start: int,
    end: int,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_markers: int = 3,
    dominance_threshold: float = 0.8,
) -> Tuple[int, int, float, int, bool]:
    """Classify (n_inf, n_dominant, dominance, mendel_errors, block_compat)."""
    in_region = [
        l for l in loci
        if l.chrom == chrom and start <= l.midpoint < end
    ]
    n_inf = 0
    n_ref = n_del = 0
    mend_err = 0
    for loc in in_region:
        calls = genotype_matrix.get(loc.marker_id, {})
        pg = calls.get(parent_sample_id, "./.")
        og = calls.get(offspring_sample_id, "./.")
        t = transmitted_from_parent(pg, og)
        if t == "contradiction":
            mend_err += 1
            continue
        if pg != "0/1":
            continue
        if t == "REF":
            n_ref += 1
            n_inf += 1
        elif t == "DEL":
            n_del += 1
            n_inf += 1
    n_dom = max(n_ref, n_del)
    dom = n_dom / n_inf if n_inf else 0.0
    block_compat = (
        n_inf >= min_markers
        and dom >= dominance_threshold
        and mend_err == 0
    )
    return n_inf, n_dom, dom, mend_err, block_compat


def classify_region_triad(
    *,
    paternal_sample_id: str,
    maternal_sample_id: str,
    offspring_sample_id: str,
    chrom: str,
    start: int,
    end: int,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_markers: int = 3,
    dominance_threshold: float = 0.8,
) -> Tuple[int, int, float, int, bool]:
    """Triad region classification — uses co-parent disambiguation
    to harvest more informative parent-HET markers."""
    in_region = [
        l for l in loci
        if l.chrom == chrom and start <= l.midpoint < end
    ]
    n_inf = 0
    n_ref = n_del = 0
    mend_err = 0
    for loc in in_region:
        calls = genotype_matrix.get(loc.marker_id, {})
        pg_pat = calls.get(paternal_sample_id, "./.")
        pg_mat = calls.get(maternal_sample_id, "./.")
        og = calls.get(offspring_sample_id, "./.")
        t_pat = transmitted_from_parent_triad(pg_pat, pg_mat, og)
        t_mat = transmitted_from_parent_triad(pg_mat, pg_pat, og)
        for t in (t_pat, t_mat):
            if t == "contradiction":
                mend_err += 1
                break
        else:
            # No contradiction. Count parent-HET resolutions.
            if pg_pat == "0/1" and t_pat in ("REF", "DEL"):
                if t_pat == "REF":
                    n_ref += 1
                else:
                    n_del += 1
                n_inf += 1
            if pg_mat == "0/1" and t_mat in ("REF", "DEL"):
                if t_mat == "REF":
                    n_ref += 1
                else:
                    n_del += 1
                n_inf += 1
    n_dom = max(n_ref, n_del)
    dom = n_dom / n_inf if n_inf else 0.0
    block_compat = (
        n_inf >= min_markers
        and dom >= dominance_threshold
        and mend_err == 0
    )
    return n_inf, n_dom, dom, mend_err, block_compat


# ----------------------------------------------------------------------
# Background-region sampling.
# ----------------------------------------------------------------------


def sample_background_regions(
    *,
    chrom: str,
    chrom_length: int,
    region_length: int,
    n_regions: int,
    exclude: Sequence[Tuple[int, int]],
    rng: random.Random,
    max_attempts: int = 100,
) -> List[Tuple[int, int]]:
    """Sample non-overlapping background windows of fixed length that do
    not intersect any excluded interval (typically the LRR list)."""
    out: List[Tuple[int, int]] = []
    excl_sorted = sorted(exclude)
    for _ in range(n_regions):
        for _attempt in range(max_attempts):
            lo = rng.randint(0, max(0, chrom_length - region_length - 1))
            hi = lo + region_length
            if not any(not (hi <= e_lo or lo >= e_hi)
                       for e_lo, e_hi in excl_sorted + out):
                out.append((lo, hi))
                break
    return out


# ----------------------------------------------------------------------
# Odds-ratio + Woolf CI with Haldane-Anscombe correction.
# ----------------------------------------------------------------------


def compute_odds_ratio(a: int, b: int, c: int, d: int) -> OddsRatio:
    haldane = False
    aa, bb, cc, dd = a, b, c, d
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        haldane = True
    try:
        odds_inside = aa / bb
        odds_outside = cc / dd
        or_val = odds_inside / odds_outside if odds_outside else None
    except ZeroDivisionError:
        or_val = None
    log_or = math.log(or_val) if or_val and or_val > 0 else None
    se = (
        math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
        if all(x > 0 for x in (aa, bb, cc, dd)) else None
    )
    if log_or is not None and se is not None:
        ci_low = math.exp(log_or - 1.96 * se)
        ci_high = math.exp(log_or + 1.96 * se)
    else:
        ci_low = ci_high = None
    return OddsRatio(
        a=a, b=b, c=c, d=d,
        odds_ratio=or_val, log_or=log_or, se_log_or=se,
        ci_low=ci_low, ci_high=ci_high,
        haldane_corrected=haldane,
    )


# ----------------------------------------------------------------------
# Aggregator across relationships + background draws.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class Relationship:
    relationship_id: str
    relationship_type: str       # "dyad" | "triad"
    paternal_sample_id: Optional[str]   # triad
    maternal_sample_id: Optional[str]   # triad
    parent_sample_id: Optional[str]     # dyad
    offspring_sample_id: str


def classify_one(
    rel: Relationship,
    chrom: str, start: int, end: int,
    region_id: str, region_kind: str,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_markers: int, dominance_threshold: float,
) -> RegionClassification:
    if rel.relationship_type == "triad":
        n_inf, n_dom, dom, err, ok = classify_region_triad(
            paternal_sample_id=rel.paternal_sample_id,
            maternal_sample_id=rel.maternal_sample_id,
            offspring_sample_id=rel.offspring_sample_id,
            chrom=chrom, start=start, end=end,
            loci=loci, genotype_matrix=genotype_matrix,
            min_markers=min_markers,
            dominance_threshold=dominance_threshold,
        )
        parent_id = f"{rel.paternal_sample_id}+{rel.maternal_sample_id}"
    else:
        n_inf, n_dom, dom, err, ok = classify_region_dyad(
            parent_sample_id=rel.parent_sample_id,
            offspring_sample_id=rel.offspring_sample_id,
            chrom=chrom, start=start, end=end,
            loci=loci, genotype_matrix=genotype_matrix,
            min_markers=min_markers,
            dominance_threshold=dominance_threshold,
        )
        parent_id = rel.parent_sample_id
    return RegionClassification(
        relationship_id=rel.relationship_id,
        relationship_type=rel.relationship_type,
        parent_sample_id=parent_id,
        offspring_sample_id=rel.offspring_sample_id,
        region_id=region_id, region_kind=region_kind,
        chrom=chrom, start=start, end=end,
        n_informative=n_inf, n_dominant_allele=n_dom,
        dominance=dom, mendelian_errors=err,
        block_compatible=ok,
    )


def compute_enrichment(
    *,
    lrrs: Sequence[LRRInterval],
    relationships: Sequence[Relationship],
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    n_background_per_lrr: int = 10,
    chrom_lengths: Optional[Dict[str, int]] = None,
    seed: int = 0,
    min_markers: int = 3,
    dominance_threshold: float = 0.8,
) -> Tuple[List[LRREnrichment], List[RegionClassification]]:
    """Per-LRR odds ratios. Returns (enrichments, raw_classifications)."""
    rng = random.Random(seed)
    # Build chrom_length defaults from the loci themselves if not provided.
    if chrom_lengths is None:
        chrom_lengths = {}
        for l in loci:
            chrom_lengths[l.chrom] = max(chrom_lengths.get(l.chrom, 0), l.end + 1)
    by_chrom_excl: Dict[str, List[Tuple[int, int]]] = {}
    for l in lrrs:
        by_chrom_excl.setdefault(l.chrom, []).append((l.start, l.end))

    classifications: List[RegionClassification] = []
    enrichments: List[LRREnrichment] = []

    for lrr in lrrs:
        # Inside-region classifications.
        for rel in relationships:
            classifications.append(classify_one(
                rel, lrr.chrom, lrr.start, lrr.end,
                region_id=lrr.lrr_id, region_kind="LRR",
                loci=loci, genotype_matrix=genotype_matrix,
                min_markers=min_markers,
                dominance_threshold=dominance_threshold,
            ))
        # Background regions: sample N matched windows on the same chrom.
        bg_regions = sample_background_regions(
            chrom=lrr.chrom,
            chrom_length=chrom_lengths.get(lrr.chrom, lrr.end * 2),
            region_length=lrr.length,
            n_regions=n_background_per_lrr,
            exclude=by_chrom_excl.get(lrr.chrom, []),
            rng=rng,
        )
        for bg_idx, (bg_lo, bg_hi) in enumerate(bg_regions):
            bg_id = f"BG_{lrr.lrr_id}_{bg_idx:03d}"
            for rel in relationships:
                classifications.append(classify_one(
                    rel, lrr.chrom, bg_lo, bg_hi,
                    region_id=bg_id, region_kind="background",
                    loci=loci, genotype_matrix=genotype_matrix,
                    min_markers=min_markers,
                    dominance_threshold=dominance_threshold,
                ))

        # Build 2×2 tables for this LRR.
        lrr_rows = [c for c in classifications if c.region_id == lrr.lrr_id]
        bg_rows = [c for c in classifications
                   if c.region_kind == "background"
                   and c.region_id.startswith(f"BG_{lrr.lrr_id}_")]

        def _table(rows_inside, rows_outside):
            a = sum(1 for r in rows_inside if r.block_compatible)
            b = sum(1 for r in rows_inside if not r.block_compatible)
            c = sum(1 for r in rows_outside if r.block_compatible)
            d = sum(1 for r in rows_outside if not r.block_compatible)
            return a, b, c, d

        a, b, c, d = _table(lrr_rows, bg_rows)
        combined = compute_odds_ratio(a, b, c, d)

        triad_inside = [r for r in lrr_rows if r.relationship_type == "triad"]
        triad_outside = [r for r in bg_rows if r.relationship_type == "triad"]
        dyad_inside = [r for r in lrr_rows if r.relationship_type == "dyad"]
        dyad_outside = [r for r in bg_rows if r.relationship_type == "dyad"]

        triad_or: Optional[OddsRatio] = None
        if triad_inside:
            ta, tb, tc, td = _table(triad_inside, triad_outside)
            triad_or = compute_odds_ratio(ta, tb, tc, td)
        dyad_or: Optional[OddsRatio] = None
        if dyad_inside:
            da, db, dc, dd_ = _table(dyad_inside, dyad_outside)
            dyad_or = compute_odds_ratio(da, db, dc, dd_)

        enrichments.append(LRREnrichment(
            lrr_id=lrr.lrr_id,
            chrom=lrr.chrom, start=lrr.start, end=lrr.end,
            n_relationships=len(relationships),
            n_triads=sum(1 for r in relationships if r.relationship_type == "triad"),
            n_dyads=sum(1 for r in relationships if r.relationship_type == "dyad"),
            combined=combined,
            triad_only=triad_or,
            dyad_only=dyad_or,
        ))

    return enrichments, classifications


# ----------------------------------------------------------------------
# LRR TSV reader.
# ----------------------------------------------------------------------


def load_lrr_list(path: str | Path) -> List[LRRInterval]:
    """TSV with header. Required columns: lrr_id, chrom, start, end."""
    path = Path(path)
    out: List[LRRInterval] = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        required = ("lrr_id", "chrom", "start", "end")
        for c in required:
            if c not in header:
                raise ValueError(f"{path}: missing required column {c!r}")
        idx = {c: header.index(c) for c in required}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0].startswith("#") or not parts[0]:
                continue
            out.append(LRRInterval(
                lrr_id=parts[idx["lrr_id"]],
                chrom=parts[idx["chrom"]],
                start=int(parts[idx["start"]]),
                end=int(parts[idx["end"]]),
            ))
    return out


# ----------------------------------------------------------------------
# TSV writer for enrichment results.
# ----------------------------------------------------------------------


def write_enrichment_tsv(path: str | Path,
                         enrichments: Sequence[LRREnrichment]) -> None:
    cols = [
        "lrr_id", "chrom", "start", "end",
        "n_relationships", "n_triads", "n_dyads",
        "combined_OR", "combined_CI_low", "combined_CI_high",
        "combined_a", "combined_b", "combined_c", "combined_d",
        "combined_haldane_corrected",
        "triad_OR", "triad_CI_low", "triad_CI_high",
        "dyad_OR", "dyad_CI_low", "dyad_CI_high",
    ]

    def _fmt(x):
        if x is None:
            return ""
        if isinstance(x, float):
            return f"{x:.4f}"
        return str(x)

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in enrichments:
            row = [
                e.lrr_id, e.chrom, e.start, e.end,
                e.n_relationships, e.n_triads, e.n_dyads,
                e.combined.odds_ratio, e.combined.ci_low, e.combined.ci_high,
                e.combined.a, e.combined.b, e.combined.c, e.combined.d,
                e.combined.haldane_corrected,
                e.triad_only.odds_ratio if e.triad_only else None,
                e.triad_only.ci_low if e.triad_only else None,
                e.triad_only.ci_high if e.triad_only else None,
                e.dyad_only.odds_ratio if e.dyad_only else None,
                e.dyad_only.ci_low if e.dyad_only else None,
                e.dyad_only.ci_high if e.dyad_only else None,
            ]
            fh.write("\t".join(_fmt(x) for x in row) + "\n")
