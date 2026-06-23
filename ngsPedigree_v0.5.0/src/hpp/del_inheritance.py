"""
LRR inheritance map (within recombination-suppressed regions).

This is the **inner-level** inheritance map: per-chromosome segment-by-
segment transmission within a contiguous arrangement haplotype block
(a candidate inversion / LRR). It tracks the *flow* of REF vs DEL
allele transmissions along contiguous markers, run-length-encodes
into segments, and flags candidate recombination break points.

The OTHER inheritance map — chromosome-level, independent-assortment
— lives in ``chromosome_inheritance.py`` (bloc 19). That module
operates at the whole-chromosome level and answers: "is this
candidate parent compatible as a contributor to this whole
chromosome of the offspring?" Both modules are valid; they answer
different questions at different scales.

```
                        |  level         |  unit        |  Pearson at
  Chromosome inheritance|  whole chrom   |  individual  |  one-pair level
  LRR inheritance       |  segment       |  population  |  one-pair level inside LRR
```

Chromosome inheritance map from DEL markers.

Given a confirmed dyad (parent → offspring) or triad (P1+P2 → O) and a
DEL-marker genotype catalogue, walk each chromosome and emit per-segment
records of which parental haplotype was transmitted. This is the
broke-grad-student replacement for ngsPedigree Stage 3 — coarser than
BEAGLE-based phasing, but it costs nothing beyond the Delly + Manta
catalogue we already have.

Parent identification by Mendelian exclusion
--------------------------------------------
For a candidate (parent → offspring) pair, "by exclusion" means: if
there exists *any* informative marker where the parent is homozygous
and the offspring is the opposite homozygote, the candidate is NOT the
parent. Concretely:

    parent = 0/0  AND  offspring = 1/1   → impossible (exclusion)
    parent = 1/1  AND  offspring = 0/0   → impossible (exclusion)

A single high-confidence opposite-homozygote excludes parentage. The
test ``exclude_as_parent`` returns the list of excluding markers (so
the caller can audit them).

Transmitted-allele inference per marker
---------------------------------------
At each marker where both parent and offspring are called:

    parent      offspring      transmitted_from_parent
    0/0         any            REF (always)
    1/1         any            DEL (always)
    0/1         0/0            REF  (offspring received parent's REF allele)
    0/1         1/1            DEL  (offspring received parent's DEL allele)
    0/1         0/1            ambiguous (parent transmitted either)

Markers where the parent is homozygous are *not* informative for
detecting recombination breakpoints — both of the parent's
haplotypes carry the same allele. Recombination tracing uses only the
parent-het markers.

Triad version
-------------
With both parents known, we can also discriminate which parent
contributed which offspring allele at parent-HET sites where the
co-parent is HOMOZYGOUS:

    P_het + P_hom_ref  →  HET offspring inherited the DEL from P_het
                          (the co-parent could only give REF)
    P_het + P_hom_del  →  HET offspring inherited the REF from P_het
                          (the co-parent could only give DEL)

So triads turn previously-ambiguous parent-HET-offspring-HET markers
into resolvable ones.

Segments
--------
Once each marker is labeled with its transmitted allele, the
chromosome is reduced to a run-length encoding. A switch from REF to
DEL (or vice versa) is a candidate recombination break. Short isolated
runs (length 1) inside a long opposite-allele tract are smoothed away
as likely noise (controlled by ``min_run_length``).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ----------------------------------------------------------------------
# Marker location helper.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class DelMarkerLocus:
    marker_id: str
    chrom: str
    start: int                  # 1-based midpoint (or breakpoint)
    end: int

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


@dataclass(frozen=True)
class DelExclusionEvidence:
    parent_sample_id: str
    offspring_sample_id: str
    n_informative: int
    n_excluding_markers: int
    excluding_marker_ids: Tuple[str, ...]
    can_be_parent: bool

    def to_dict(self) -> dict:
        d = asdict(self)
        d["excluding_marker_ids"] = list(self.excluding_marker_ids)
        return d


@dataclass(frozen=True)
class TransmittedAllele:
    marker_id: str
    chrom: str
    pos: int                     # marker midpoint
    transmitted: str             # "REF" | "DEL" | "ambiguous" | "contradiction"
    parent_gt: str
    offspring_gt: str
    informative_for_recomb: bool # True only when parent_gt == "0/1"
                                  # AND transmission is resolved


@dataclass(frozen=True)
class InheritanceSegment:
    dyad_or_triad_id: str
    parent_sample_id: str
    offspring_sample_id: str
    chrom: str
    seg_start: int                  # 1-based
    seg_end: int                    # 1-based
    transmitted_allele: str         # "REF" | "DEL" | "ambiguous"
    n_informative_markers: int
    confidence: str                 # "Gold" | "Silver" | "Bronze"
    recomb_event_left: bool
    recomb_event_right: bool
    parental_hap_inherited: str     # "1" | "2" | "ambiguous"
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ----------------------------------------------------------------------
# Parent-by-exclusion test.
# ----------------------------------------------------------------------


def exclude_as_parent(
    parent_sample_id: str,
    offspring_sample_id: str,
    genotype_matrix: Dict[str, Dict[str, str]],
    *,
    min_excluding: int = 1,
) -> DelExclusionEvidence:
    """Identify a putative parent by Mendelian exclusion.

    Returns the count of markers where (parent_gt, offspring_gt) makes
    parentage impossible. The default policy (``min_excluding=1``) is
    strict: a single opposite-homozygote excludes the candidate. Set
    higher for short-read noise tolerance.
    """
    excl: List[str] = []
    n_inf = 0
    for marker_id, calls in genotype_matrix.items():
        pg = calls.get(parent_sample_id, "./.")
        og = calls.get(offspring_sample_id, "./.")
        if pg in ("./.", "") or og in ("./.", ""):
            continue
        n_inf += 1
        if (pg == "0/0" and og == "1/1") or (pg == "1/1" and og == "0/0"):
            excl.append(marker_id)
    return DelExclusionEvidence(
        parent_sample_id=parent_sample_id,
        offspring_sample_id=offspring_sample_id,
        n_informative=n_inf,
        n_excluding_markers=len(excl),
        excluding_marker_ids=tuple(excl),
        can_be_parent=(len(excl) < min_excluding),
    )


# ----------------------------------------------------------------------
# Transmitted-allele inference.
# ----------------------------------------------------------------------


def transmitted_from_parent(parent_gt: str, offspring_gt: str) -> str:
    """Return REF | DEL | ambiguous | contradiction | unknown."""
    if parent_gt in ("./.", "") or offspring_gt in ("./.", ""):
        return "unknown"
    if parent_gt == "0/0":
        if offspring_gt == "1/1":
            return "contradiction"
        return "REF"
    if parent_gt == "1/1":
        if offspring_gt == "0/0":
            return "contradiction"
        return "DEL"
    # parent_gt == "0/1"
    if offspring_gt == "0/0":
        return "REF"
    if offspring_gt == "1/1":
        return "DEL"
    return "ambiguous"


def transmitted_from_parent_triad(
    parent_gt: str, co_parent_gt: str, offspring_gt: str,
) -> str:
    """Triad-disambiguated transmission.

    Same as ``transmitted_from_parent`` for resolvable cases, but at
    parent-HET sites where the dyad call is ambiguous, the co-parent's
    homozygous state lets us pick.
    """
    base = transmitted_from_parent(parent_gt, offspring_gt)
    if base != "ambiguous":
        return base
    # parent_gt == "0/1" and offspring_gt == "0/1"
    if co_parent_gt == "0/0":
        # co-parent could only give REF → parent must have given DEL
        return "DEL"
    if co_parent_gt == "1/1":
        # co-parent could only give DEL → parent must have given REF
        return "REF"
    return "ambiguous"


# ----------------------------------------------------------------------
# Walk markers along one chromosome and yield TransmittedAllele rows.
# ----------------------------------------------------------------------


def walk_chromosome(
    *,
    parent_sample_id: str,
    offspring_sample_id: str,
    chrom: str,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    co_parent_sample_id: Optional[str] = None,
) -> List[TransmittedAllele]:
    out: List[TransmittedAllele] = []
    chrom_loci = [l for l in loci if l.chrom == chrom]
    chrom_loci.sort(key=lambda l: l.midpoint)
    for loc in chrom_loci:
        calls = genotype_matrix.get(loc.marker_id, {})
        pg = calls.get(parent_sample_id, "./.")
        og = calls.get(offspring_sample_id, "./.")
        if co_parent_sample_id:
            cg = calls.get(co_parent_sample_id, "./.")
            t = transmitted_from_parent_triad(pg, cg, og)
        else:
            t = transmitted_from_parent(pg, og)
        if t == "unknown":
            continue
        out.append(TransmittedAllele(
            marker_id=loc.marker_id,
            chrom=loc.chrom,
            pos=loc.midpoint,
            transmitted=t,
            parent_gt=pg,
            offspring_gt=og,
            informative_for_recomb=(pg == "0/1" and t in ("REF", "DEL")),
        ))
    return out


# ----------------------------------------------------------------------
# Run-length encode into inheritance segments.
# ----------------------------------------------------------------------


def _smooth_runs(
    values: Sequence[str], *, min_run_length: int = 2,
) -> List[Tuple[str, int, int]]:
    """Return [(value, start_idx, end_idx_inclusive)] over the input list,
    after collapsing isolated runs of length < ``min_run_length`` that
    sit inside a longer opposing run.

    Conservative smoother — only flips isolated minority entries flanked
    on both sides by the same majority value. Avoids over-smoothing real
    short blocks.
    """
    if not values:
        return []
    vals = list(values)
    # naive run-length encoding
    runs: List[List[int]] = []   # each: [value_index_into_vals_at_start, start, end]
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]:
            j += 1
        runs.append([i, i, j])
        i = j + 1

    # smooth single-marker minority runs flanked by same neighbour
    changed = True
    while changed:
        changed = False
        for k in range(1, len(runs) - 1):
            run_start, _, run_end = runs[k]
            length = run_end - run_start + 1
            if length >= min_run_length:
                continue
            left_val = vals[runs[k - 1][1]]
            right_val = vals[runs[k + 1][1]]
            cur_val = vals[run_start]
            if left_val == right_val and left_val != cur_val:
                for x in range(run_start, run_end + 1):
                    vals[x] = left_val
                changed = True
                break
        if changed:
            # re-encode runs
            runs = []
            i = 0
            while i < len(vals):
                j = i
                while j + 1 < len(vals) and vals[j + 1] == vals[i]:
                    j += 1
                runs.append([i, i, j])
                i = j + 1

    return [(vals[r[1]], r[1], r[2]) for r in runs]


def build_chromosome_segments(
    *,
    dyad_or_triad_id: str,
    parent_sample_id: str,
    offspring_sample_id: str,
    transmissions: Sequence[TransmittedAllele],
    min_run_length: int = 2,
) -> List[InheritanceSegment]:
    """Group transmissions on one chromosome into inheritance segments."""
    if not transmissions:
        return []
    # Use only parent-het informative markers for recombination tracing.
    # All-het / all-hom-parent chromosomes give a single ambiguous segment.
    informative = [t for t in transmissions if t.informative_for_recomb]
    if not informative:
        # The whole chromosome is uninformative for breakpoints. Emit one
        # segment spanning the chrom.
        first, last = transmissions[0], transmissions[-1]
        # Convention for parental_hap when not informative: ambiguous.
        return [InheritanceSegment(
            dyad_or_triad_id=dyad_or_triad_id,
            parent_sample_id=parent_sample_id,
            offspring_sample_id=offspring_sample_id,
            chrom=first.chrom,
            seg_start=first.pos,
            seg_end=last.pos,
            transmitted_allele="ambiguous",
            n_informative_markers=0,
            confidence="Bronze",
            recomb_event_left=False, recomb_event_right=False,
            parental_hap_inherited="ambiguous",
            notes="no parent-het markers on this chromosome",
        )]

    values = [t.transmitted for t in informative]
    runs = _smooth_runs(values, min_run_length=min_run_length)

    segments: List[InheritanceSegment] = []
    n_runs = len(runs)
    for k, (val, a, b) in enumerate(runs):
        seg_start = informative[a].pos
        seg_end = informative[b].pos
        n_markers = b - a + 1
        # Confidence by marker density: Gold ≥ 8, Silver 3-7, Bronze 1-2
        confidence = ("Gold" if n_markers >= 8
                       else "Silver" if n_markers >= 3
                       else "Bronze")
        parental_hap = (
            "1" if val == "REF"
            else "2" if val == "DEL"
            else "ambiguous"
        )
        segments.append(InheritanceSegment(
            dyad_or_triad_id=dyad_or_triad_id,
            parent_sample_id=parent_sample_id,
            offspring_sample_id=offspring_sample_id,
            chrom=informative[0].chrom,
            seg_start=seg_start, seg_end=seg_end,
            transmitted_allele=val,
            n_informative_markers=n_markers,
            confidence=confidence,
            recomb_event_left=(k > 0),
            recomb_event_right=(k < n_runs - 1),
            parental_hap_inherited=parental_hap,
        ))
    return segments


# ----------------------------------------------------------------------
# Top-level driver — emit segments per chromosome for a dyad or triad.
# ----------------------------------------------------------------------


def build_inheritance_map_for_dyad(
    *,
    dyad_id: str,
    parent_sample_id: str,
    offspring_sample_id: str,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_run_length: int = 2,
) -> List[InheritanceSegment]:
    chroms = sorted({l.chrom for l in loci})
    segments: List[InheritanceSegment] = []
    for chrom in chroms:
        t = walk_chromosome(
            parent_sample_id=parent_sample_id,
            offspring_sample_id=offspring_sample_id,
            chrom=chrom,
            loci=loci,
            genotype_matrix=genotype_matrix,
        )
        segments.extend(build_chromosome_segments(
            dyad_or_triad_id=dyad_id,
            parent_sample_id=parent_sample_id,
            offspring_sample_id=offspring_sample_id,
            transmissions=t,
            min_run_length=min_run_length,
        ))
    return segments


def build_inheritance_map_for_triad(
    *,
    triad_id: str,
    paternal_sample_id: str,
    maternal_sample_id: str,
    offspring_sample_id: str,
    loci: Sequence[DelMarkerLocus],
    genotype_matrix: Dict[str, Dict[str, str]],
    min_run_length: int = 2,
) -> Dict[str, List[InheritanceSegment]]:
    """Return ``{"paternal": [...], "maternal": [...]}`` — one inheritance
    trace per parent. Triad mode resolves more parent-HET markers via
    the co-parent's homozygous state."""
    chroms = sorted({l.chrom for l in loci})
    out: Dict[str, List[InheritanceSegment]] = {"paternal": [], "maternal": []}
    for chrom in chroms:
        t_pat = walk_chromosome(
            parent_sample_id=paternal_sample_id,
            offspring_sample_id=offspring_sample_id,
            chrom=chrom,
            loci=loci,
            genotype_matrix=genotype_matrix,
            co_parent_sample_id=maternal_sample_id,
        )
        t_mat = walk_chromosome(
            parent_sample_id=maternal_sample_id,
            offspring_sample_id=offspring_sample_id,
            chrom=chrom,
            loci=loci,
            genotype_matrix=genotype_matrix,
            co_parent_sample_id=paternal_sample_id,
        )
        out["paternal"].extend(build_chromosome_segments(
            dyad_or_triad_id=triad_id,
            parent_sample_id=paternal_sample_id,
            offspring_sample_id=offspring_sample_id,
            transmissions=t_pat,
            min_run_length=min_run_length,
        ))
        out["maternal"].extend(build_chromosome_segments(
            dyad_or_triad_id=triad_id,
            parent_sample_id=maternal_sample_id,
            offspring_sample_id=offspring_sample_id,
            transmissions=t_mat,
            min_run_length=min_run_length,
        ))
    return out
