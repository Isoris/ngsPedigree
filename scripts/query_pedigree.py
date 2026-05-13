#!/usr/bin/env python3
"""
query_pedigree.py — query the annotated pedigree graph from either direction.

Usage examples:
  # Forward: who's the parent of hub H001?
  python query_pedigree.py hub H001 --outdir <ngspedigree_output_dir>

  # Reverse: who are the candidate parents of CGA187?
  python query_pedigree.py sample CGA187 --outdir <ngspedigree_output_dir>

  # List all PO edges
  python query_pedigree.py edges --class parent_offspring --outdir <...>

  # List all hubs of a given type
  python query_pedigree.py hubs --type two_parents_with_sibship --outdir <...>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_artifacts(outdir):
    outdir = Path(outdir)
    edges = pd.read_csv(outdir / "pairwise_relationship_classification.tsv", sep="\t")
    roster = pd.read_csv(outdir / "family_hub_roster.tsv", sep="\t")
    return edges, roster


def query_sample(sample_id, edges, roster):
    """Reverse query: everything about this sample."""
    print(f"\n=== Sample: {sample_id} ===")

    in_roster = roster[roster["sample_id"] == sample_id]
    if in_roster.empty:
        print(f"  Sample not found in roster.")
        return
    r = in_roster.iloc[0]
    print(f"  Hub:        {r['hub_id']} ({r['hub_type']}, size {r['hub_size']})")
    print(f"  Role:       {r['possible_role']} (confidence: {r['role_confidence']})")
    print(f"  Reason:     {r['reason']}")

    # First-degree neighbors
    pairs_a = edges[(edges["sample_a"] == sample_id) &
                    (edges["edge_class"].isin(["parent_offspring", "full_sibling",
                                                "duplicate_or_clone", "ambiguous_first_degree"]))]
    pairs_b = edges[(edges["sample_b"] == sample_id) &
                    (edges["edge_class"].isin(["parent_offspring", "full_sibling",
                                                "duplicate_or_clone", "ambiguous_first_degree"]))]
    neighbors = []
    for _, row in pairs_a.iterrows():
        neighbors.append((row["sample_b"], row["edge_class"], row["confidence"], row["theta"], row["IBS0"]))
    for _, row in pairs_b.iterrows():
        neighbors.append((row["sample_a"], row["edge_class"], row["confidence"], row["theta"], row["IBS0"]))

    if not neighbors:
        print(f"  No first-degree neighbors.")
        return

    print(f"\n  First-degree neighbors ({len(neighbors)}):")
    print(f"    {'sample':12s}  {'class':22s}  {'confidence':10s}  {'theta':>8s}  {'IBS0':>8s}")
    for nb_id, ec, conf, theta, ibs0 in sorted(neighbors, key=lambda x: x[1]):
        print(f"    {nb_id:12s}  {ec:22s}  {conf:10s}  {theta:>8.4f}  {ibs0:>8.4f}")

    # Specifically: candidate parents
    po_neighbors = [n for n in neighbors if n[1] == "parent_offspring"]
    if po_neighbors:
        print(f"\n  Candidate parents (PO neighbors): {[n[0] for n in po_neighbors]}")
        print(f"    Direction NOT determined here; check hub roster for parent role.")


def query_hub(hub_id, edges, roster):
    """Forward query: everything about this hub."""
    print(f"\n=== Hub: {hub_id} ===")
    members = roster[roster["hub_id"] == hub_id]
    if members.empty:
        print(f"  Hub not found.")
        return
    hub_type = members.iloc[0]["hub_type"]
    print(f"  Type: {hub_type}")
    print(f"  Members ({len(members)}):")
    print(f"    {'sample':12s}  {'role':30s}  {'confidence':10s}  reason")
    for _, r in members.iterrows():
        print(f"    {r['sample_id']:12s}  {r['possible_role']:30s}  {r['role_confidence']:10s}  {r['reason']}")

    # Edges within hub
    hub_samples = set(members["sample_id"])
    in_hub = edges[edges["sample_a"].isin(hub_samples) & edges["sample_b"].isin(hub_samples)]
    print(f"\n  Edges within hub ({len(in_hub)}):")
    for _, e in in_hub.iterrows():
        print(f"    {e['sample_a']:12s} -- {e['sample_b']:12s}  {e['edge_class']:22s}  ({e['confidence']})")


def list_edges(edges, edge_class=None):
    if edge_class:
        sub = edges[edges["edge_class"] == edge_class]
        print(f"\n=== Edges of class {edge_class} ({len(sub)}) ===")
    else:
        sub = edges
        print(f"\n=== All edges ({len(sub)}) ===")
    for _, e in sub.iterrows():
        print(f"  {e['sample_a']:12s} -- {e['sample_b']:12s}  {e['edge_class']:22s}  "
              f"theta={e['theta']:.4f}  IBS0={e['IBS0']:.4f}  ({e['confidence']})")


def list_hubs(roster, hub_type=None):
    if hub_type:
        sub = roster[roster["hub_type"] == hub_type]
    else:
        sub = roster
    hub_ids = sub["hub_id"].unique()
    print(f"\n=== Hubs ({len(hub_ids)}) ===")
    for hid in sorted(hub_ids):
        if hid == "_isolated":
            continue
        members = sub[sub["hub_id"] == hid]
        ht = members.iloc[0]["hub_type"]
        size = members.iloc[0]["hub_size"]
        sample_list = ", ".join(sorted(members["sample_id"]))
        print(f"  {hid}  ({ht}, n={size})  {sample_list}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mode", choices=["sample", "hub", "edges", "hubs"],
                   help="query mode")
    p.add_argument("target", nargs="?", default=None,
                   help="sample_id (for 'sample') or hub_id (for 'hub'); ignored for 'edges'/'hubs'")
    p.add_argument("--outdir", required=True, help="ngsPedigree output directory")
    p.add_argument("--class", dest="edge_class", default=None,
                   help="for 'edges' mode: filter to this edge_class")
    p.add_argument("--type", dest="hub_type", default=None,
                   help="for 'hubs' mode: filter to this hub_type")
    args = p.parse_args()

    edges, roster = load_artifacts(args.outdir)

    if args.mode == "sample":
        if not args.target:
            sys.stderr.write("ERROR: 'sample' mode requires sample_id\n")
            sys.exit(1)
        query_sample(args.target, edges, roster)
    elif args.mode == "hub":
        if not args.target:
            sys.stderr.write("ERROR: 'hub' mode requires hub_id\n")
            sys.exit(1)
        query_hub(args.target, edges, roster)
    elif args.mode == "edges":
        list_edges(edges, args.edge_class)
    elif args.mode == "hubs":
        list_hubs(roster, args.hub_type)


if __name__ == "__main__":
    main()
