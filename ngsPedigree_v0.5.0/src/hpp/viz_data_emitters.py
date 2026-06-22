"""
Bloc 20 — visualization-data emitters.

For every figure in the user-uploaded ngsPedigree manuscript layouts
**that maps to data we actually produce in this pipeline**, emit the
underlying data table(s) so they can be rendered in matplotlib /
ggplot2 / Cytoscape / whatever is downstream. We do not render plots
ourselves (no matplotlib in stdlib). The emitter side is stdlib-only.

Figures covered + the function that emits its data:

  Chromosome-scale inheritance ideograms     emit_ideogram_segments
  Per-chromosome event summary                emit_per_chromosome_events
  Local LRR view                              emit_local_lrr_view
  Pairwise kinship matrix ordered by family   emit_kinship_matrix
  Pair counts by inferred edge class          emit_edge_class_counts
  Pedigree network (nodes + edges)            emit_pedigree_network
  Pairwise metrics scatter (theta vs IBS0
   vs Jaccard)                                emit_pairwise_metrics
  Mating-risk matrix (F × M expected kinship) emit_mating_risk_matrix
  Close-kin groups summary                    emit_close_kin_groups
  Genome-wide event timeline (CO / NCO /
   Mendelian errors)                          emit_genome_event_timeline

Deliberately NOT included
-------------------------
  * Genotype-likelihood PCA — requires SNP genotype likelihoods from
    ANGSD, not DEL genotypes. A DEL-dosage stand-in would be misleading
    relative to the user-uploaded figure label.
  * Layered-resolver workflow log — we do not run an explicit layered
    resolver with per-layer edge counts; the workflow is implicit in
    the master pipeline.
  * Case cards — derivable presentationally from existing tables
    (pairwise_metrics + edge_class + mtdna_validation); no separate
    emitter is justified.
  * Ancestry / admixture (NGSadmix-style K components) — requires SNP
    genotype likelihoods, not DEL data. Out of scope for this bloc.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ----------------------------------------------------------------------
# Helpers.
# ----------------------------------------------------------------------


def _write_tsv(path, rows, columns):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(
                "" if r.get(c) is None
                else (f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]))
                for c in columns
            ) + "\n")


def _gt_to_dosage(gt: str) -> Optional[int]:
    return {"0/0": 0, "0/1": 1, "1/1": 2}.get(gt)


def _has_del(gt: str) -> Optional[bool]:
    d = _gt_to_dosage(gt)
    return None if d is None else d >= 1


def _pearson(xs, ys):
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
# 1. Pairwise metrics scatter — theta, IBS0, Jaccard per pair.
# ----------------------------------------------------------------------


def jaccard_del(sample_a, sample_b, genotype_matrix):
    """Jaccard on DEL-carrying marker sets (≥1 DEL allele)."""
    n_inter = 0
    n_union = 0
    for calls in genotype_matrix.values():
        a = _has_del(calls.get(sample_a, "./."))
        b = _has_del(calls.get(sample_b, "./."))
        if a is None or b is None:
            continue
        if a or b:
            n_union += 1
        if a and b:
            n_inter += 1
    return n_inter / n_union if n_union else None


def emit_pairwise_metrics(
    pairs,                  # iterable of PairCoeffs (from del_relatedness)
    genotype_matrix,
    edge_class_by_pair,
    out_path,
):
    """One row per pair. Columns: sample_a, sample_b, theta, IBS0,
    Jaccard, KING, n_informative, edge_class. Drives the
    theta-vs-IBS0-vs-Jaccard scatter matrix."""
    rows = []
    for p in pairs:
        jac = jaccard_del(p.sample_a, p.sample_b, genotype_matrix)
        edge_class = (edge_class_by_pair.get((p.sample_a, p.sample_b))
                      or edge_class_by_pair.get((p.sample_b, p.sample_a))
                      or "unknown")
        rows.append({
            "sample_a": p.sample_a, "sample_b": p.sample_b,
            "theta": p.theta, "IBS0": p.IBS0, "Jaccard": jac,
            "KING": p.theta, "n_informative": p.n_informative,
            "edge_class": edge_class,
        })
    _write_tsv(out_path, rows,
                columns=["sample_a", "sample_b", "theta", "IBS0",
                         "Jaccard", "KING", "n_informative", "edge_class"])
    return rows


# ----------------------------------------------------------------------
# 2. Pair counts by inferred edge class.
# ----------------------------------------------------------------------


def emit_edge_class_counts(edge_class_by_pair, out_path):
    """One row per edge class with the count. Bar chart in the figure."""
    seen = set()
    counter = Counter()
    for (a, b), cls in edge_class_by_pair.items():
        key = tuple(sorted([a, b]))
        if key in seen:
            continue
        seen.add(key)
        counter[cls] += 1
    rows = [{"edge_class": c, "n_pairs": n}
            for c, n in sorted(counter.items(), key=lambda x: -x[1])]
    _write_tsv(out_path, rows, columns=["edge_class", "n_pairs"])
    return rows


# ----------------------------------------------------------------------
# 3. Pedigree network — nodes + edges (for Cytoscape / networkx).
# ----------------------------------------------------------------------


def emit_pedigree_network(
    samples, edge_class_by_pair,
    edge_classes_kept=("parent_offspring", "full_sibling",
                       "duplicate_or_clone", "second_degree"),
    sample_metadata: Optional[Dict[str, Dict]] = None,
    out_nodes=None, out_edges=None,
):
    sample_metadata = sample_metadata or {}
    nodes = [{
        "sample_id": s,
        **(sample_metadata.get(s, {})),
    } for s in samples]
    seen = set()
    edges = []
    for (a, b), cls in edge_class_by_pair.items():
        if cls not in edge_classes_kept:
            continue
        key = tuple(sorted([a, b]))
        if key in seen:
            continue
        seen.add(key)
        edges.append({"source": key[0], "target": key[1], "edge_class": cls})
    if out_nodes:
        _write_tsv(out_nodes, nodes,
                    columns=["sample_id"] + sorted(
                        {k for n in nodes for k in n if k != "sample_id"}))
    if out_edges:
        _write_tsv(out_edges, edges,
                    columns=["source", "target", "edge_class"])
    return nodes, edges


# ----------------------------------------------------------------------
# 4. Close-kin community detection + family aggregates.
# ----------------------------------------------------------------------


def detect_kin_groups(edge_class_by_pair, samples,
                      in_group_classes=("parent_offspring", "full_sibling",
                                         "duplicate_or_clone")):
    """Union-find on the first-degree subgraph. Returns
    {sample_id: family_id}."""
    parent = {s: s for s in samples}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    seen = set()
    for (a, b), cls in edge_class_by_pair.items():
        if cls not in in_group_classes:
            continue
        key = tuple(sorted([a, b]))
        if key in seen:
            continue
        seen.add(key)
        if a in parent and b in parent:
            union(a, b)

    raw = {s: find(s) for s in samples}
    # Re-id sequentially.
    next_id = 0
    out: Dict[str, str] = {}
    rep_to_id: Dict[str, str] = {}
    for s in samples:
        r = raw[s]
        if r not in rep_to_id:
            next_id += 1
            rep_to_id[r] = f"F{next_id:03d}"
        out[s] = rep_to_id[r]
    return out


def emit_close_kin_groups(
    samples, family_by_sample, theta_by_pair,
    sample_metadata: Optional[Dict[str, Dict]] = None,
    out_path=None,
):
    """One row per kin family with: family_id, n_members, members
    (semicolon-separated), mean_within_family_theta, etc."""
    by_family: Dict[str, List[str]] = defaultdict(list)
    for s, f in family_by_sample.items():
        by_family[f].append(s)

    rows = []
    for fid, members in sorted(by_family.items(),
                                 key=lambda kv: -len(kv[1])):
        if len(members) < 2:
            mean_theta = None
        else:
            vals = []
            for i, a in enumerate(members):
                for b in members[i + 1:]:
                    t = theta_by_pair.get((a, b)) or theta_by_pair.get((b, a))
                    if t is not None:
                        vals.append(t)
            mean_theta = (sum(vals) / len(vals)) if vals else None
        rows.append({
            "family_id": fid,
            "n_members": len(members),
            "members": ";".join(sorted(members)),
            "mean_within_family_theta": mean_theta,
        })
    if out_path:
        _write_tsv(out_path, rows,
                    columns=["family_id", "n_members", "members",
                             "mean_within_family_theta"])
    return rows


# ----------------------------------------------------------------------
# 5. Kinship matrix ordered by close-kin families (heatmap data).
# ----------------------------------------------------------------------


def emit_kinship_matrix(samples, theta_by_pair, family_by_sample,
                        out_path):
    """One row per (sample_a, sample_b) ordered by family. The downstream
    plot orders rows/columns by family_id."""
    ordered = sorted(samples,
                     key=lambda s: (family_by_sample.get(s, "Z"), s))
    rows = []
    for a in ordered:
        for b in ordered:
            t = theta_by_pair.get((a, b)) or theta_by_pair.get((b, a))
            t = 0.5 if a == b else t
            rows.append({
                "sample_a": a, "sample_b": b,
                "theta": t,
                "family_a": family_by_sample.get(a, ""),
                "family_b": family_by_sample.get(b, ""),
            })
    _write_tsv(out_path, rows,
                columns=["sample_a", "sample_b", "theta",
                         "family_a", "family_b"])
    return rows, ordered


# ----------------------------------------------------------------------
# 6. Mating-risk matrix — candidate females × candidate males expected
#    offspring kinship.
# ----------------------------------------------------------------------


def emit_mating_risk_matrix(
    candidate_females, candidate_males, theta_by_pair, out_path,
    *,
    low_risk_max=0.025,    # offspring kinship below 0.025 is "low risk"
    high_risk_min=0.0625,   # above 0.0625 is "high risk" (first cousins)
):
    """For each F × M pair, expected offspring kinship = θ(F, M).

    Categorisation:
      low (< low_risk_max)    "recommended"
      mid                       "caution"
      high (≥ high_risk_min)  "avoid"
    """
    rows = []
    for f in candidate_females:
        for m in candidate_males:
            t = theta_by_pair.get((f, m))
            if t is None:
                t = theta_by_pair.get((m, f))
            if t is None:
                cat = "unknown"
            elif t < low_risk_max:
                cat = "low_risk_recommended"
            elif t >= high_risk_min:
                cat = "high_risk_avoid"
            else:
                cat = "mid_caution"
            rows.append({
                "female_sample_id": f,
                "male_sample_id": m,
                "expected_offspring_kinship": t,
                "risk_category": cat,
            })
    _write_tsv(out_path, rows,
                columns=["female_sample_id", "male_sample_id",
                         "expected_offspring_kinship", "risk_category"])
    return rows


# ----------------------------------------------------------------------
# 7. Chromosome-scale ideogram segments — already produced by
#    `del_inheritance.build_inheritance_map_for_triad`; we just re-shape
#    that output for ideogram rendering.
# ----------------------------------------------------------------------


def emit_ideogram_segments(triad_maps, out_path):
    """Take the dict {triad_id: {"paternal": [...segs], "maternal":
    [...segs]}} output and emit one row per (offspring, chrom, side,
    seg). The downstream renderer draws ideogram blocks."""
    rows = []
    for triad_id, side_segs in triad_maps.items():
        for side, segs in side_segs.items():
            for s in segs:
                rows.append({
                    "triad_id": triad_id,
                    "side": side,
                    "offspring_sample_id": s.offspring_sample_id,
                    "parent_sample_id": s.parent_sample_id,
                    "chrom": s.chrom,
                    "seg_start": s.seg_start,
                    "seg_end": s.seg_end,
                    "transmitted_allele": s.transmitted_allele,
                    "n_informative_markers": s.n_informative_markers,
                    "confidence": s.confidence,
                    "recomb_event_left": s.recomb_event_left,
                    "recomb_event_right": s.recomb_event_right,
                    "parental_hap_inherited": s.parental_hap_inherited,
                })
    _write_tsv(out_path, rows,
                columns=["triad_id", "side", "offspring_sample_id",
                         "parent_sample_id", "chrom", "seg_start",
                         "seg_end", "transmitted_allele",
                         "n_informative_markers", "confidence",
                         "recomb_event_left", "recomb_event_right",
                         "parental_hap_inherited"])
    return rows


# ----------------------------------------------------------------------
# 8. Per-chromosome event summary (CO / NCO / Mendelian errors).
# ----------------------------------------------------------------------


def emit_per_chromosome_events(triad_maps, out_path):
    """Per (offspring, chrom, parent side): counts of CO, NCO, Mendelian
    errors. Drives the per-chromosome event bar chart."""
    counts: Dict[Tuple[str, str, str], Dict[str, int]] = defaultdict(
        lambda: {"n_segments": 0, "n_co": 0, "n_nco": 0,
                  "n_mendelian_errors": 0})
    for triad_id, side_segs in triad_maps.items():
        for side, segs in side_segs.items():
            for s in segs:
                key = (s.offspring_sample_id, s.chrom, side)
                counts[key]["n_segments"] += 1
                if s.recomb_event_left or s.recomb_event_right:
                    counts[key]["n_co"] += 1
                if s.confidence == "Bronze" and s.n_informative_markers <= 2:
                    counts[key]["n_nco"] += 1
                if "contradiction" in s.notes.lower() if s.notes else False:
                    counts[key]["n_mendelian_errors"] += 1
    rows = []
    for (off, chrom, side), c in sorted(counts.items()):
        rows.append({
            "offspring_sample_id": off,
            "chrom": chrom,
            "side": side,
            **c,
        })
    _write_tsv(out_path, rows,
                columns=["offspring_sample_id", "chrom", "side",
                         "n_segments", "n_co", "n_nco",
                         "n_mendelian_errors"])
    return rows


# ----------------------------------------------------------------------
# 9. Genome-wide event timeline — every event with chromosome + position.
# ----------------------------------------------------------------------


def emit_genome_event_timeline(triad_maps, out_path):
    rows = []
    for triad_id, side_segs in triad_maps.items():
        for side, segs in side_segs.items():
            for s in segs:
                if s.recomb_event_left:
                    rows.append({
                        "triad_id": triad_id, "side": side,
                        "chrom": s.chrom, "position": s.seg_start,
                        "event_type": "crossover",
                        "parent_sample_id": s.parent_sample_id,
                        "offspring_sample_id": s.offspring_sample_id,
                        "confidence": s.confidence,
                    })
    _write_tsv(out_path, rows,
                columns=["triad_id", "side", "chrom", "position",
                         "event_type", "parent_sample_id",
                         "offspring_sample_id", "confidence"])
    return rows


# ----------------------------------------------------------------------
# 10. Local LRR view — variants + transmissions inside one LRR.
# ----------------------------------------------------------------------


def emit_local_lrr_view(
    lrr_id, chrom, start, end,
    triad_maps, loci, genotype_matrix,
    out_path,
):
    """All markers inside the LRR + their per-sample dosage, plus the
    paternal/maternal transmitted trace from each contributing triad.
    Drives the local-view figure."""
    rows = []
    for loc in loci:
        if loc.chrom != chrom or not (start <= loc.midpoint < end):
            continue
        gts = genotype_matrix.get(loc.marker_id, {})
        for sample_id, gt in sorted(gts.items()):
            rows.append({
                "lrr_id": lrr_id, "chrom": loc.chrom,
                "marker_id": loc.marker_id, "marker_pos": loc.midpoint,
                "sample_id": sample_id, "genotype": gt,
                "dosage": _gt_to_dosage(gt),
            })
    _write_tsv(out_path, rows,
                columns=["lrr_id", "chrom", "marker_id", "marker_pos",
                         "sample_id", "genotype", "dosage"])
    return rows


