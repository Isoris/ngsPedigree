"""
Chromosome inheritance map — independent-assortment level (bloc 19).

The OTHER inheritance map lives in ``del_inheritance.py``; that one is
the **LRR inheritance map** — segment-by-segment transmission within
a recombination-suppressed haplotype block (a candidate inversion /
LRR). It tracks the *flow* of REF vs DEL allele transmissions along
contiguous markers, run-length-encodes into segments, and flags
candidate recombination break points.

THIS module operates one level up the hierarchy. It treats each whole
chromosome as the unit and asks: **is this candidate parent compatible
as a contributor to this whole chromosome of the offspring?** Both
maps are valid; they answer different questions.

```
                        |  level         |  unit        |  Pearson at
  Chromosome inheritance|  whole chrom   |  individual  |  one-pair level
  LRR inheritance       |  segment       |  population  |  ditto, inside LRR
```

The chromosome-level signal
---------------------------
A diploid offspring receives one homolog of each chromosome from each
parent (independent chromosomal assortment). For an autosome:

  - if candidate parent is hom-REF (0/0) at marker M, the homolog that
    parent contributed carries 0. Offspring must therefore carry
    **at least one 0** at M (genotype 0/0 or 0/1).
  - if candidate parent is hom-DEL (1/1), the contributed homolog
    carries 1. Offspring must carry **at least one 1** (0/1 or 1/1).
  - if candidate parent is HET, the contributed homolog carries 0 or 1
    with equal probability — uninformative for hard exclusion at a
    single marker, but the *parent's HET DELs* are nevertheless a
    fingerprint: ~50% of parent's HET DELs should be present in the
    offspring's chromosome.

A true parent has zero (or near-zero) opposite-homozygote
contradictions on every chromosome, *and* the Pearson correlation
between parent and offspring DEL dosage on that chromosome is
markedly positive (~0.5 in expectation under independent assortment
for a single chromosome's worth of markers — higher when the same
homolog runs the whole chromosome). A non-parent shows opposite-hom
contradictions on at least some chromosomes and near-zero Pearson r.

Use cases
---------
  - Triad sanity-check: each of the two confirmed parents should
    score high (compatible) on every chromosome.
  - Distinguishing the true parent from a half-sibling that looks
    first-degree at θ: half-sibs share only ~25% of segregating
    markers, so a per-chromosome score is much lower than a PO score.
  - Cryptic chromosomal events: a uniparental disomy or a sibship
    swap leaves one chromosome with the wrong parent profile while
    others look normal — visible per-chromosome but invisible to
    genome-wide θ.

This is the broke-grad-student stand-in for a SNP-resolution
phased-haplotype map. With ~200k+ DEL markers across the genome we
have ~1k–10k DELs per chromosome — plenty of signal for the Pearson
score to be stable on real data.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .del_inheritance import DelMarkerLocus


# ----------------------------------------------------------------------
# Records.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class ChromosomeInheritanceScore:
    offspring_sample_id: str
    candidate_parent_sample_id: str
    chrom: str
    n_markers_called: int                # both samples called
    n_compatible: int                    # no opposite-hom contradictions
    n_excluding: int                     # opposite-hom contradictions
    compatibility_rate: Optional[float]
    pearson_r: Optional[float]
    # Parent-HET-marker transmission rate (counts how many of parent's
    # HET DELs the offspring carries at least one copy of).
    n_parent_het_markers: int
    n_parent_het_inherited: int
    het_inheritance_rate: Optional[float]
    inheritance_support: str             # see _classify_support
    parent_hom_dels_inherited: int       # n parent-1/1 markers where offspring also carries ≥1 DEL
    parent_hom_dels_total: int

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ChromosomeInheritanceMap:
    """Map for ONE offspring: scores against every candidate parent
    on every chromosome, with the best candidate per chromosome flagged."""
    offspring_sample_id: str
    candidate_parents: List[str]
    scores: List[ChromosomeInheritanceScore] = field(default_factory=list)
    best_per_chrom: Dict[str, str] = field(default_factory=dict)
    best_score_per_chrom: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "offspring_sample_id": self.offspring_sample_id,
            "candidate_parents": list(self.candidate_parents),
            "scores": [s.to_dict() for s in self.scores],
            "best_per_chrom": dict(self.best_per_chrom),
            "best_score_per_chrom": dict(self.best_score_per_chrom),
        }


# ----------------------------------------------------------------------
# Genotype helpers (DEL dosage 0/1/2; missing → None).
# ----------------------------------------------------------------------


def _dosage(gt: str) -> Optional[int]:
    return {"0/0": 0, "0/1": 1, "1/1": 2}.get(gt)


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> Optional[float]:
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    return num / math.sqrt(sxx * syy)


# ----------------------------------------------------------------------
# Single (offspring, parent, chromosome) score.
# ----------------------------------------------------------------------


def _classify_support(
    *,
    n_called: int,
    n_excluding: int,
    pearson_r: Optional[float],
    het_inheritance_rate: Optional[float],
    min_markers: int,
) -> str:
    """Bucket the score:
      'rejected'    — ≥1 opposite-hom contradiction (strict Mendelian)
      'ambiguous'   — too few markers to call
      'compatible'  — Mendelian-clean but no positive signal
      'strong'      — Mendelian-clean with high Pearson r and/or
                       parent-HET inheritance rate near 0.5
    """
    if n_called < min_markers:
        return "ambiguous"
    if n_excluding > 0:
        return "rejected"
    # Mendelian-clean. Look for positive support.
    pos_r = pearson_r is not None and pearson_r >= 0.30
    pos_het = (het_inheritance_rate is not None
                and 0.35 <= het_inheritance_rate <= 0.65)
    if pos_r and pos_het:
        return "strong"
    if pos_r or pos_het:
        return "compatible"
    return "compatible"


def score_chromosome_inheritance(
    *,
    offspring_sample_id: str,
    candidate_parent_sample_id: str,
    chrom: str,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_markers: int = 20,
) -> ChromosomeInheritanceScore:
    par_dos: List[int] = []
    off_dos: List[int] = []
    n_excluding = 0
    n_compat = 0
    n_called = 0
    n_par_het = 0
    n_par_het_inherited = 0
    n_par_hom_total = 0
    n_par_hom_inherited = 0

    for loc in loci:
        if loc.chrom != chrom:
            continue
        gts = genotype_matrix.get(loc.marker_id, {})
        pg = gts.get(candidate_parent_sample_id, "./.")
        og = gts.get(offspring_sample_id, "./.")
        pd = _dosage(pg)
        od = _dosage(og)
        if pd is None or od is None:
            continue
        n_called += 1
        par_dos.append(pd)
        off_dos.append(od)
        # Mendelian-exclusion check (opposite-homozygote).
        if (pg == "0/0" and og == "1/1") or (pg == "1/1" and og == "0/0"):
            n_excluding += 1
        else:
            n_compat += 1
        # Parent-HET inheritance counting.
        if pg == "0/1":
            n_par_het += 1
            if od >= 1:
                n_par_het_inherited += 1
        # Parent-HOM-DEL passage.
        if pg == "1/1":
            n_par_hom_total += 1
            if od >= 1:
                n_par_hom_inherited += 1

    compat_rate = (n_compat / n_called) if n_called else None
    r = _pearson(par_dos, off_dos) if n_called >= 3 else None
    het_rate = (n_par_het_inherited / n_par_het) if n_par_het else None
    support = _classify_support(
        n_called=n_called, n_excluding=n_excluding,
        pearson_r=r, het_inheritance_rate=het_rate,
        min_markers=min_markers,
    )
    return ChromosomeInheritanceScore(
        offspring_sample_id=offspring_sample_id,
        candidate_parent_sample_id=candidate_parent_sample_id,
        chrom=chrom,
        n_markers_called=n_called,
        n_compatible=n_compat,
        n_excluding=n_excluding,
        compatibility_rate=compat_rate,
        pearson_r=r,
        n_parent_het_markers=n_par_het,
        n_parent_het_inherited=n_par_het_inherited,
        het_inheritance_rate=het_rate,
        inheritance_support=support,
        parent_hom_dels_inherited=n_par_hom_inherited,
        parent_hom_dels_total=n_par_hom_total,
    )


# ----------------------------------------------------------------------
# Per-offspring map across all candidate parents and chromosomes.
# ----------------------------------------------------------------------


def build_chromosome_inheritance_map(
    *,
    offspring_sample_id: str,
    candidate_parents: Sequence[str],
    chromosomes: Sequence[str],
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_markers: int = 20,
) -> ChromosomeInheritanceMap:
    """Score every (candidate parent × chromosome) pair for one offspring.

    Returns the full score list plus best-parent-per-chromosome lookup
    (best by Pearson r among non-rejected candidates; ties broken by
    higher compatibility rate then higher n_markers_called).
    """
    m = ChromosomeInheritanceMap(
        offspring_sample_id=offspring_sample_id,
        candidate_parents=list(candidate_parents),
    )
    for chrom in chromosomes:
        for parent in candidate_parents:
            if parent == offspring_sample_id:
                continue
            s = score_chromosome_inheritance(
                offspring_sample_id=offspring_sample_id,
                candidate_parent_sample_id=parent,
                chrom=chrom, loci=loci,
                genotype_matrix=genotype_matrix,
                min_markers=min_markers,
            )
            m.scores.append(s)

    # Pick best candidate per chromosome.
    by_chrom: Dict[str, List[ChromosomeInheritanceScore]] = {}
    for s in m.scores:
        by_chrom.setdefault(s.chrom, []).append(s)
    for chrom, ss in by_chrom.items():
        non_rejected = [s for s in ss if s.inheritance_support != "rejected"]
        if not non_rejected:
            continue
        # Sort by (Pearson r desc, compatibility rate desc, n_called desc).
        non_rejected.sort(key=lambda s: (
            -(s.pearson_r if s.pearson_r is not None else -1),
            -(s.compatibility_rate if s.compatibility_rate is not None else -1),
            -s.n_markers_called,
        ))
        best = non_rejected[0]
        m.best_per_chrom[chrom] = best.candidate_parent_sample_id
        m.best_score_per_chrom[chrom] = best.pearson_r or 0.0
    return m


# ----------------------------------------------------------------------
# TSV writer.
# ----------------------------------------------------------------------


def write_chromosome_inheritance_tsv(
    path, maps: Iterable[ChromosomeInheritanceMap],
) -> None:
    cols = ["offspring_sample_id", "candidate_parent_sample_id", "chrom",
            "n_markers_called", "n_compatible", "n_excluding",
            "compatibility_rate", "pearson_r",
            "n_parent_het_markers", "n_parent_het_inherited",
            "het_inheritance_rate",
            "parent_hom_dels_inherited", "parent_hom_dels_total",
            "inheritance_support"]
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
        for m in maps:
            for s in m.scores:
                fh.write("\t".join(_fmt(getattr(s, c)) for c in cols) + "\n")


def write_chromosome_best_parent_tsv(
    path, maps: Iterable[ChromosomeInheritanceMap],
) -> None:
    """Per (offspring, chromosome) flat table: best candidate parent +
    pearson_r."""
    cols = ["offspring_sample_id", "chrom", "best_parent", "best_pearson_r"]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m in maps:
            for chrom in sorted(m.best_per_chrom):
                fh.write(
                    f"{m.offspring_sample_id}\t{chrom}\t"
                    f"{m.best_per_chrom[chrom]}\t"
                    f"{m.best_score_per_chrom.get(chrom, 0.0):.4f}\n"
                )
