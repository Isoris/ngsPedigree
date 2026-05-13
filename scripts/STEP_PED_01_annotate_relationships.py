#!/usr/bin/env python3
"""
ngsPedigree — STEP_PED_01_annotate_relationships.py
====================================================

Stage 1 of ngsPedigree: solve and annotate the first-degree relatedness graph
from a genome-wide ngsRelate .res file.

Philosophy
----------
ngsPedigree is a graph annotator, not a pedigree reconstructor. It takes the
filtered first-degree subgraph (theta >= 0.177 by default) and labels every
edge with a relationship class and every node with its possible role within
its hub. Nothing is filtered out; everything is labeled. Downstream consumers
can then query the annotated graph from either direction (forward: who is
the parent of this hub? reverse: who are this sample's candidate parents?).

Operates in BLIND MODE by default — no metadata required. Optional --sex TSV
is the only metadata input and is used purely to promote `parent_a` / `parent_b`
labels to `mother` / `father` where graph topology has already determined
parental roles. Sex is never used to flip direction or change edge classification.

Inputs
------
  --res PATH          : ngsRelate output file (.res, 23-col standard format)
  --samples PATH      : sample sidecar TSV, one sample_id per line, in
                        ngsRelate's `a`/`b` index order
  --outdir PATH       : output directory for the three artifacts
  [--sex PATH]        : optional 2-col TSV (sample_id\\tsex) where sex in
                        {male, female, unknown}. Only used for label promotion.
  [--run-id STRING]   : provenance string written into the envelope
                        (default: cohort_unspecified_v1)
  [--theta-first PCT] : first-degree threshold (default 0.177, Manichaikul)
  [--theta-second]    : second-degree threshold (default 0.0884)
  [--theta-third]     : third-degree threshold (default 0.0442)
  [--ibs0-po-max]     : max IBS0 for PO classification (default 0.005)
  [--theta-dup-min]   : min theta for duplicate/clone (default 0.45)

Outputs
-------
  pairwise_relationship_classification.tsv : edge-level annotations
  family_hub_roster.tsv                    : node-level annotations
  ngspedigree_run_envelope.json            : self-describing wrapper

See docs/SCHEMA.md for full column documentation.

Usage
-----
  python STEP_PED_01_annotate_relationships.py \\
      --res data/cohort/relatedness/ngsrelate/cohort_226_full_v1/relatedness.res \\
      --samples data/cohort/relatedness/ngsrelate/cohort_226_full_v1/samples.txt \\
      --outdir data/cohort/relatedness/ngsrelate/cohort_226_full_v1/ngspedigree/ \\
      --run-id cohort_226_full_v1

Author: ngsPedigree v0.1, 2026-05-10
"""

import argparse
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Hard dependency: pandas. Soft on numpy (only used for nan handling).
try:
    import pandas as pd
except ImportError:
    sys.stderr.write("ERROR: pandas is required. Install with: pip install pandas\n")
    sys.exit(2)


SCRIPT_VERSION = "v0.1.0"
SCRIPT_NAME    = "STEP_PED_01_annotate_relationships"

# ---------------------------------------------------------------------------
# ngsRelate column conventions
# ---------------------------------------------------------------------------
#
# The standard 23-column ngsRelate output (ANGSD-Korneliussen v2) includes:
#   a, b, nSites, J7, J8, J9, rab, Fa, Fb, theta, inbreed_a, inbreed_b,
#   2of3_IDB, FDiff, loglh, nIter, coverage, IBS0, IBS1, IBS2, R0, R1, KING
#
# The classifier needs at minimum: theta, IBS0, KING. R0, R1, J9 are used
# for confidence scoring and the duplicate-detection branch. The script will
# error cleanly if the minimum set is missing and warn if the optional set is.

REQUIRED_COLS = ["a", "b", "theta", "IBS0"]
RECOMMENDED_COLS = ["KING", "R0", "R1", "J9", "nSites"]

# ---------------------------------------------------------------------------
# Edge classifier
# ---------------------------------------------------------------------------
#
# This is the core function that STEP_PED_02 will reuse for per-chromosome
# annotation. Keep it pure: takes a row dict + thresholds, returns a tuple.
# No I/O, no globals.
#
# Decision order (FIRST MATCH WINS — order matters):
#   1. duplicate_or_clone   — theta >= theta_dup_min AND IBS0 < 0.001
#                              (catch this first; would otherwise hit PO)
#   2. parent_offspring     — theta in [theta_first, theta_dup_min) AND IBS0 < ibs0_po_max
#   3. full_sibling         — theta in [theta_first, theta_dup_min) AND IBS0 >= ibs0_po_max
#   4. ambiguous_first_deg  — theta in [theta_first, theta_dup_min) but IBS0 missing
#   5. second_degree        — theta in [theta_second, theta_first)
#   6. third_degree         — theta in [theta_third, theta_second)
#   7. unrelated            — theta < theta_third

def classify_edge(row, thresholds):
    """Classify a single ngsRelate row.

    Args:
        row: dict-like with at least 'theta' and 'IBS0' (may be NaN).
             Optionally 'KING', 'R0', 'R1', 'J9', 'nSites'.
        thresholds: dict with theta_first, theta_second, theta_third,
                    theta_dup_min, ibs0_po_max.

    Returns:
        dict with:
          edge_class: str (one of the 7 classes)
          confidence: str (high|medium|low)
          reasons: str (semicolon-separated diagnostic notes)
    """
    theta = row.get("theta")
    ibs0  = row.get("IBS0")
    king  = row.get("KING")
    n_sites = row.get("nSites")

    reasons = []
    confidence = "high"

    # Guard: missing theta means we can't classify at all.
    if theta is None or pd.isna(theta):
        return dict(edge_class="undetermined",
                    confidence="low",
                    reasons="theta_missing")

    # Low n_sites is an automatic confidence downgrade.
    if n_sites is not None and not pd.isna(n_sites) and n_sites < 1000:
        confidence = "low"
        reasons.append(f"low_n_sites_{int(n_sites)}")

    # ---- Branch 1: duplicate / clone / monozygotic-equivalent ----
    if theta >= thresholds["theta_dup_min"]:
        if ibs0 is None or pd.isna(ibs0):
            reasons.append("ibs0_missing_dup_inferred_from_theta")
            confidence = "medium" if confidence == "high" else confidence
        elif ibs0 >= 0.001:
            # High theta but non-trivial IBS0 — weird; flag.
            reasons.append("high_theta_unexpected_ibs0")
            confidence = "low"
        return dict(edge_class="duplicate_or_clone",
                    confidence=confidence,
                    reasons=";".join(reasons) if reasons else "")

    # ---- Branch 2/3: first-degree (PO vs FS) ----
    if theta >= thresholds["theta_first"]:
        # Need IBS0 to separate PO from FS.
        if ibs0 is None or pd.isna(ibs0):
            reasons.append("ibs0_missing_cannot_separate_po_fs")
            return dict(edge_class="ambiguous_first_degree",
                        confidence="low",
                        reasons=";".join(reasons))
        if ibs0 < thresholds["ibs0_po_max"]:
            # PO. Confidence high if KING also supports first-degree (>= 0.177).
            if king is not None and not pd.isna(king) and king < 0.15:
                reasons.append(f"king_low_for_first_degree_{king:.3f}")
                confidence = "medium" if confidence == "high" else confidence
            return dict(edge_class="parent_offspring",
                        confidence=confidence,
                        reasons=";".join(reasons) if reasons else "")
        else:
            # FS (or rare other configurations like double-cousin in inbred lines).
            return dict(edge_class="full_sibling",
                        confidence=confidence,
                        reasons=";".join(reasons) if reasons else "")

    # ---- Branch 5: second-degree ----
    if theta >= thresholds["theta_second"]:
        return dict(edge_class="second_degree",
                    confidence=confidence,
                    reasons=";".join(reasons) if reasons else "")

    # ---- Branch 6: third-degree ----
    if theta >= thresholds["theta_third"]:
        return dict(edge_class="third_degree",
                    confidence=confidence,
                    reasons=";".join(reasons) if reasons else "")

    # ---- Branch 7: unrelated ----
    return dict(edge_class="unrelated",
                confidence=confidence,
                reasons=";".join(reasons) if reasons else "")


# ---------------------------------------------------------------------------
# Hub-topology pass (blind-mode role inference)
# ---------------------------------------------------------------------------
#
# After every edge is classified, we look at the structure of each connected
# component (treating only edges with edge_class in {parent_offspring,
# full_sibling, duplicate_or_clone} as in-hub edges) and assign roles to
# each node. Direction is decidable from topology in some cases (the
# forced-parent rule below), uncertain in others (po_dyad_only, mixed).
#
# Forced-parent rule:
#   If a node A has parent_offspring edges to >= 2 other nodes B, C, ...
#   AND those nodes B, C, ... form a full-sibling clique among themselves,
#   then A is the parent. (B, C, ... cannot all be parents of A; an
#   individual has only two parents.)
#
# Hub types:
#   two_parents_with_sibship  : two unrelated parents, shared sibship
#   parent_with_sibship       : one parent, sibship of >= 2 offspring
#   po_dyad_only              : single PO edge, no third sample to triangulate
#   sibship_only              : all FS edges, no PO edges (parents not sampled)
#   duplicate_pair            : just a dup edge
#   mixed_or_complex          : anything that doesn't fit the patterns above

def find_connected_components(edges_in_hub, all_node_ids):
    """Return list of components. Each component is a set of sample_ids."""
    parent = {n: n for n in all_node_ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for a, b in edges_in_hub:
        if a in parent and b in parent:
            union(a, b)

    comps = defaultdict(set)
    for n in all_node_ids:
        comps[find(n)].add(n)
    return list(comps.values())


def infer_hub_topology(component_nodes, edges_by_pair):
    """Given a set of node IDs and a map (a,b)->edge_class, infer hub structure.

    Returns dict with:
      hub_type         : str
      node_roles       : {sample_id: role}
      role_confidence  : {sample_id: high|medium|low}
      reasons          : {sample_id: str}
    """
    nodes = sorted(component_nodes)
    n = len(nodes)

    # Single node: not a hub.
    if n == 1:
        return dict(hub_type="singleton",
                    node_roles={nodes[0]: "isolated"},
                    role_confidence={nodes[0]: "high"},
                    reasons={nodes[0]: "no_first_degree_edges"})

    # Two nodes: classify by edge type.
    if n == 2:
        a, b = nodes
        edge = edges_by_pair.get((a, b)) or edges_by_pair.get((b, a))
        edge_class = edge["edge_class"] if edge else None

        if edge_class == "duplicate_or_clone":
            return dict(hub_type="duplicate_pair",
                        node_roles={a: "duplicate", b: "duplicate"},
                        role_confidence={a: "high", b: "high"},
                        reasons={a: "duplicate_pair", b: "duplicate_pair"})
        if edge_class == "parent_offspring":
            # Direction not decidable from a single edge.
            return dict(hub_type="po_dyad_only",
                        node_roles={a: "ambiguous_first_degree_PO",
                                    b: "ambiguous_first_degree_PO"},
                        role_confidence={a: "low", b: "low"},
                        reasons={a: "dyad_no_triangulation",
                                 b: "dyad_no_triangulation"})
        if edge_class == "full_sibling":
            return dict(hub_type="sibship_only",
                        node_roles={a: "possible_full_sib",
                                    b: "possible_full_sib"},
                        role_confidence={a: "medium", b: "medium"},
                        reasons={a: "fs_dyad", b: "fs_dyad"})
        return dict(hub_type="mixed_or_complex",
                    node_roles={a: "ambiguous", b: "ambiguous"},
                    role_confidence={a: "low", b: "low"},
                    reasons={a: f"unexpected_edge_class_{edge_class}",
                             b: f"unexpected_edge_class_{edge_class}"})

    # Three or more nodes: count PO and FS edges per node within the component.
    po_neighbors = defaultdict(set)
    fs_neighbors = defaultdict(set)
    dup_neighbors = defaultdict(set)

    for a in nodes:
        for b in nodes:
            if a == b:
                continue
            edge = edges_by_pair.get((a, b)) or edges_by_pair.get((b, a))
            if not edge:
                continue
            ec = edge["edge_class"]
            if ec == "parent_offspring":
                po_neighbors[a].add(b)
            elif ec == "full_sibling":
                fs_neighbors[a].add(b)
            elif ec == "duplicate_or_clone":
                dup_neighbors[a].add(b)

    # Find candidate parent nodes: nodes with PO edges to >= 2 others, where
    # those others form an FS clique among themselves.
    candidate_parents = []
    for a in nodes:
        po_n = po_neighbors[a]
        if len(po_n) < 2:
            continue
        # Check that po_n forms an FS clique.
        po_n_list = sorted(po_n)
        is_fs_clique = True
        for i, x in enumerate(po_n_list):
            for y in po_n_list[i+1:]:
                if y not in fs_neighbors[x]:
                    is_fs_clique = False
                    break
            if not is_fs_clique:
                break
        if is_fs_clique:
            candidate_parents.append((a, len(po_n)))

    node_roles = {}
    role_confidence = {}
    reasons = {}

    if len(candidate_parents) >= 2:
        # Two parents — verify they are unrelated to each other.
        candidate_parents.sort(key=lambda x: -x[1])
        p1, p2 = candidate_parents[0][0], candidate_parents[1][0]
        edge = edges_by_pair.get((p1, p2)) or edges_by_pair.get((p2, p1))
        if edge is None or edge["edge_class"] in ("unrelated", "third_degree", "second_degree"):
            # two_parents_with_sibship
            shared_offspring = po_neighbors[p1] & po_neighbors[p2]
            for nd in nodes:
                if nd == p1:
                    node_roles[nd] = "parent_a"
                    role_confidence[nd] = "high" if len(po_neighbors[p1]) >= 3 else "medium"
                    reasons[nd] = f"po_to_{len(po_neighbors[p1])}_fs_clique;two_parent_topology"
                elif nd == p2:
                    node_roles[nd] = "parent_b"
                    role_confidence[nd] = "high" if len(po_neighbors[p2]) >= 3 else "medium"
                    reasons[nd] = f"po_to_{len(po_neighbors[p2])}_fs_clique;two_parent_topology"
                elif nd in shared_offspring:
                    node_roles[nd] = "possible_offspring"
                    role_confidence[nd] = "high"
                    reasons[nd] = f"po_to_both_parents_p1_{p1}_p2_{p2}"
                else:
                    node_roles[nd] = "ambiguous_in_two_parent_hub"
                    role_confidence[nd] = "low"
                    reasons[nd] = "not_po_to_both_inferred_parents"
            return dict(hub_type="two_parents_with_sibship",
                        node_roles=node_roles,
                        role_confidence=role_confidence,
                        reasons=reasons)
        # else: two PO-anchor nodes that are themselves related. Falls through to mixed.

    if len(candidate_parents) == 1:
        # parent_with_sibship.
        p = candidate_parents[0][0]
        offspring = po_neighbors[p]
        for nd in nodes:
            if nd == p:
                node_roles[nd] = "forced_parent" if len(offspring) >= 3 else "likely_parent"
                role_confidence[nd] = "high" if len(offspring) >= 3 else "medium"
                reasons[nd] = f"po_to_{len(offspring)}_fs_clique"
            elif nd in offspring:
                node_roles[nd] = "possible_offspring"
                role_confidence[nd] = "medium"
                reasons[nd] = f"po_to_inferred_parent_{p}"
            else:
                node_roles[nd] = "ambiguous_in_parent_hub"
                role_confidence[nd] = "low"
                reasons[nd] = "in_component_but_not_offspring_of_inferred_parent"
        return dict(hub_type="parent_with_sibship",
                    node_roles=node_roles,
                    role_confidence=role_confidence,
                    reasons=reasons)

    # No PO-anchor found. Could be sibship_only, or mixed.
    has_any_po = any(po_neighbors[n] for n in nodes)
    has_any_fs = any(fs_neighbors[n] for n in nodes)

    if not has_any_po and has_any_fs:
        # Pure sibship.
        for nd in nodes:
            node_roles[nd] = "possible_full_sib"
            role_confidence[nd] = "medium"
            reasons[nd] = "fs_clique_no_parents_sampled"
        return dict(hub_type="sibship_only",
                    node_roles=node_roles,
                    role_confidence=role_confidence,
                    reasons=reasons)

    # Mixed / complex.
    for nd in nodes:
        n_po = len(po_neighbors[nd])
        n_fs = len(fs_neighbors[nd])
        node_roles[nd] = "ambiguous_first_degree"
        role_confidence[nd] = "low"
        reasons[nd] = f"po_deg_{n_po};fs_deg_{n_fs};no_clean_topology"
    return dict(hub_type="mixed_or_complex",
                node_roles=node_roles,
                role_confidence=role_confidence,
                reasons=reasons)


# ---------------------------------------------------------------------------
# Sex-based label promotion (only metadata path)
# ---------------------------------------------------------------------------

def promote_with_sex(roster_df, sex_map):
    """If sex info exists, promote parent_a/parent_b/forced_parent/likely_parent
    to mother/father where unambiguous.

    Mutates roster_df in place. Returns count of promotions.
    """
    if sex_map is None:
        return 0

    promotions = 0
    warnings_per_hub = defaultdict(list)

    # First pass: collect parent_a / parent_b assignments per hub.
    for hub_id, hub_df in roster_df.groupby("hub_id"):
        if hub_id == "_isolated":
            continue
        parents = hub_df[hub_df["possible_role"].isin(
            ["parent_a", "parent_b", "forced_parent", "likely_parent"])]
        sexes = {row["sample_id"]: sex_map.get(row["sample_id"], "unknown")
                 for _, row in parents.iterrows()}

        # two-parent hub with one each → mother/father
        if len(sexes) == 2:
            sids = list(sexes.keys())
            sxs = [sexes[s] for s in sids]
            if set(sxs) == {"male", "female"}:
                for sid, sx in zip(sids, sxs):
                    new_role = "mother" if sx == "female" else "father"
                    roster_df.loc[roster_df["sample_id"] == sid, "possible_role"] = new_role
                    promotions += 1
            elif "male" in sxs and sxs.count("male") == 2:
                warnings_per_hub[hub_id].append("same_sex_male_parent_pair_unusual")
            elif "female" in sxs and sxs.count("female") == 2:
                warnings_per_hub[hub_id].append("same_sex_female_parent_pair_unusual")
            # else: at least one unknown — leave roles as-is

        # single-parent hub → promote if sex known
        elif len(sexes) == 1:
            sid, sx = next(iter(sexes.items()))
            if sx == "female":
                roster_df.loc[roster_df["sample_id"] == sid, "possible_role"] = "mother"
                promotions += 1
            elif sx == "male":
                roster_df.loc[roster_df["sample_id"] == sid, "possible_role"] = "father"
                promotions += 1

    # Append warnings.
    for hub_id, warns in warnings_per_hub.items():
        warn_str = ";".join(warns)
        mask = (roster_df["hub_id"] == hub_id)
        roster_df.loc[mask, "reason"] = roster_df.loc[mask, "reason"].astype(str) + ";" + warn_str

    return promotions


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_samples(path):
    """Read sample sidecar; one sample_id per line. Returns list."""
    samples = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                samples.append(s)
    return samples


def load_res(path):
    """Read ngsRelate .res file (TSV with header). Returns DataFrame."""
    df = pd.read_csv(path, sep="\t")
    # Some ngsRelate versions use 'IBS0' or 'ibs0'; normalize.
    df.columns = [c.strip() for c in df.columns]
    rename_map = {}
    for c in df.columns:
        if c.lower() == "ibs0":
            rename_map[c] = "IBS0"
        elif c.lower() == "ibs1":
            rename_map[c] = "IBS1"
        elif c.lower() == "ibs2":
            rename_map[c] = "IBS2"
        elif c.lower() == "king":
            rename_map[c] = "KING"
        elif c.lower() == "nsites":
            rename_map[c] = "nSites"
    df = df.rename(columns=rename_map)
    return df


def load_sex(path):
    """Read 2-col sex TSV. Returns {sample_id: sex} or None."""
    if path is None:
        return None
    sex_map = {}
    df = pd.read_csv(path, sep="\t")
    if not {"sample_id", "sex"}.issubset(df.columns):
        sys.stderr.write(f"WARNING: --sex file {path} missing required columns; ignored\n")
        return None
    for _, row in df.iterrows():
        sx = str(row["sex"]).strip().lower()
        if sx in ("male", "female"):
            sex_map[str(row["sample_id"]).strip()] = sx
    return sex_map if sex_map else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--res", required=True, help="ngsRelate .res file")
    p.add_argument("--samples", required=True, help="sample sidecar (one ID per line)")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--sex", default=None, help="optional sex TSV (sample_id\\tsex)")
    p.add_argument("--run-id", default="cohort_unspecified_v1",
                   help="provenance string for the envelope")
    p.add_argument("--theta-first", type=float, default=0.177,
                   help="first-degree threshold (default 0.177)")
    p.add_argument("--theta-second", type=float, default=0.0884,
                   help="second-degree threshold (default 0.0884)")
    p.add_argument("--theta-third", type=float, default=0.0442,
                   help="third-degree threshold (default 0.0442)")
    p.add_argument("--ibs0-po-max", type=float, default=0.005,
                   help="max IBS0 for PO classification (default 0.005)")
    p.add_argument("--theta-dup-min", type=float, default=0.45,
                   help="min theta for duplicate (default 0.45)")
    args = p.parse_args()

    thresholds = dict(
        theta_first=args.theta_first,
        theta_second=args.theta_second,
        theta_third=args.theta_third,
        theta_dup_min=args.theta_dup_min,
        ibs0_po_max=args.ibs0_po_max,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load inputs ----
    sys.stderr.write(f"[ngsPedigree] Loading samples from {args.samples}\n")
    samples = load_samples(args.samples)
    sys.stderr.write(f"[ngsPedigree]   {len(samples)} samples\n")

    sys.stderr.write(f"[ngsPedigree] Loading ngsRelate output from {args.res}\n")
    df = load_res(args.res)
    sys.stderr.write(f"[ngsPedigree]   {len(df)} pairs, {len(df.columns)} columns\n")

    # ---- Validate columns ----
    missing_required = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing_required:
        sys.stderr.write(f"FATAL: required columns missing from .res: {missing_required}\n")
        sys.stderr.write(f"  found columns: {list(df.columns)}\n")
        sys.stderr.write(f"  re-run ngsRelate with the standard 23-column output.\n")
        sys.exit(3)
    missing_recommended = [c for c in RECOMMENDED_COLS if c not in df.columns]
    if missing_recommended:
        sys.stderr.write(f"[ngsPedigree] WARN: recommended columns missing: {missing_recommended}\n")
        sys.stderr.write(f"[ngsPedigree]   classifier will fall back to theta+IBS0 only; confidence may be lower\n")

    # ---- Load optional sex map ----
    sex_map = load_sex(args.sex)
    if sex_map:
        sys.stderr.write(f"[ngsPedigree] Loaded sex info for {len(sex_map)} samples\n")

    # ---- Map ngsRelate indices to sample IDs ----
    if df["a"].max() >= len(samples) or df["b"].max() >= len(samples):
        sys.stderr.write(f"FATAL: .res indices out of bounds for samples sidecar\n")
        sys.stderr.write(f"  max a={df['a'].max()}, max b={df['b'].max()}, n_samples={len(samples)}\n")
        sys.exit(4)
    df["sample_a"] = df["a"].map(lambda i: samples[i])
    df["sample_b"] = df["b"].map(lambda i: samples[i])

    # ---- Classify every edge ----
    sys.stderr.write(f"[ngsPedigree] Classifying {len(df)} edges...\n")
    classifications = df.apply(lambda r: classify_edge(r.to_dict(), thresholds), axis=1)
    df["edge_class"] = classifications.map(lambda c: c["edge_class"])
    df["confidence"] = classifications.map(lambda c: c["confidence"])
    df["reasons"]    = classifications.map(lambda c: c["reasons"])

    # ---- Build edges_by_pair lookup for hub topology ----
    edges_by_pair = {}
    for _, row in df.iterrows():
        edges_by_pair[(row["sample_a"], row["sample_b"])] = {
            "edge_class": row["edge_class"],
            "confidence": row["confidence"],
        }

    # ---- Build first-degree subgraph and find connected components ----
    in_hub_classes = {"parent_offspring", "full_sibling", "duplicate_or_clone",
                      "ambiguous_first_degree"}
    in_hub_edges = [(r["sample_a"], r["sample_b"]) for _, r in df.iterrows()
                    if r["edge_class"] in in_hub_classes]

    components = find_connected_components(in_hub_edges, samples)
    # Sort components: largest first; isolated singletons last.
    multi_comps = [c for c in components if len(c) > 1]
    singletons = [c for c in components if len(c) == 1]
    multi_comps.sort(key=lambda c: -len(c))
    sys.stderr.write(f"[ngsPedigree] First-degree graph: {len(multi_comps)} multi-node components, "
                     f"{len(singletons)} singletons\n")
    if multi_comps:
        sys.stderr.write(f"[ngsPedigree]   largest component: n={len(multi_comps[0])}\n")

    # ---- Infer hub topology + node roles ----
    sys.stderr.write(f"[ngsPedigree] Inferring hub topology...\n")
    roster_rows = []
    hub_types = defaultdict(int)

    for hub_idx, comp in enumerate(multi_comps):
        hub_id = f"H{hub_idx+1:03d}"
        topo = infer_hub_topology(comp, edges_by_pair)
        hub_types[topo["hub_type"]] += 1
        for sid in comp:
            roster_rows.append({
                "sample_id": sid,
                "hub_id": hub_id,
                "hub_type": topo["hub_type"],
                "hub_size": len(comp),
                "possible_role": topo["node_roles"][sid],
                "role_confidence": topo["role_confidence"][sid],
                "reason": topo["reasons"][sid],
            })

    for comp in singletons:
        sid = next(iter(comp))
        roster_rows.append({
            "sample_id": sid,
            "hub_id": "_isolated",
            "hub_type": "isolated",
            "hub_size": 1,
            "possible_role": "isolated",
            "role_confidence": "high",
            "reason": "no_first_degree_edges",
        })

    roster_df = pd.DataFrame(roster_rows)

    # ---- Sex-based label promotion (optional) ----
    n_promotions = promote_with_sex(roster_df, sex_map)
    if n_promotions > 0:
        sys.stderr.write(f"[ngsPedigree] Promoted {n_promotions} parent labels to mother/father using --sex\n")

    # ---- Write outputs ----
    edge_out_cols = ["sample_a", "sample_b", "a", "b", "nSites", "theta", "IBS0",
                     "KING", "R0", "R1", "J9",
                     "edge_class", "confidence", "reasons"]
    edge_out_cols = [c for c in edge_out_cols if c in df.columns]

    pairwise_path = outdir / "pairwise_relationship_classification.tsv"
    df[edge_out_cols].to_csv(pairwise_path, sep="\t", index=False)
    sys.stderr.write(f"[ngsPedigree] Wrote {pairwise_path}\n")

    roster_path = outdir / "family_hub_roster.tsv"
    roster_df.to_csv(roster_path, sep="\t", index=False)
    sys.stderr.write(f"[ngsPedigree] Wrote {roster_path}\n")

    # ---- Build envelope ----
    n_pairs_by_class = df["edge_class"].value_counts().to_dict()
    envelope = {
        "schema": "ngspedigree_run_envelope_v1",
        "produced_by": {
            "tool": SCRIPT_NAME,
            "version": SCRIPT_VERSION,
            "params": {
                "thresholds": thresholds,
                "sex_provided": sex_map is not None,
                "n_sex_samples_provided": len(sex_map) if sex_map else 0,
            },
        },
        "inputs": {
            "res_path": str(args.res),
            "samples_path": str(args.samples),
            "sex_path": str(args.sex) if args.sex else None,
        },
        "run_id": args.run_id,
        "mode": "sex_assisted" if sex_map else "blind",
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "n_samples": len(samples),
        "n_pairs": int(len(df)),
        "n_pairs_by_class": {k: int(v) for k, v in n_pairs_by_class.items()},
        "n_hubs_by_type": {k: int(v) for k, v in hub_types.items()},
        "n_multi_node_components": len(multi_comps),
        "n_singletons": len(singletons),
        "n_label_promotions_via_sex": n_promotions,
        "artifacts": {
            "pairwise_relationship_classification": str(pairwise_path.name),
            "family_hub_roster": str(roster_path.name),
        },
        "schema_status": "v1_draft",
        "downstream": {
            "next_step": "STEP_PED_02_per_chromosome_annotation.py",
            "consumers": ["family_segregation_gate (Slice 3)",
                          "STEP_PED_03_inheritance_map.py",
                          "ngsTracts (NCO/DCO calling)"],
        },
    }

    envelope_path = outdir / "ngspedigree_run_envelope.json"
    with open(envelope_path, "w") as fh:
        json.dump(envelope, fh, indent=2)
    sys.stderr.write(f"[ngsPedigree] Wrote {envelope_path}\n")

    # ---- Summary to stderr ----
    sys.stderr.write(f"\n[ngsPedigree] ===== SUMMARY =====\n")
    sys.stderr.write(f"  mode: {envelope['mode']}\n")
    sys.stderr.write(f"  n_samples: {len(samples)}\n")
    sys.stderr.write(f"  n_pairs: {len(df)}\n")
    sys.stderr.write(f"  edge classes:\n")
    for k, v in sorted(n_pairs_by_class.items(), key=lambda x: -x[1]):
        sys.stderr.write(f"    {k:30s} {v}\n")
    sys.stderr.write(f"  hub types:\n")
    for k, v in sorted(hub_types.items(), key=lambda x: -x[1]):
        sys.stderr.write(f"    {k:30s} {v}\n")
    sys.stderr.write(f"  outputs in: {outdir}\n")
    sys.stderr.write(f"\n")


if __name__ == "__main__":
    main()
