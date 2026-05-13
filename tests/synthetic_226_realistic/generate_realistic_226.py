#!/usr/bin/env python3
"""
Generate a realistic synthetic 226-sample ngsRelate .res that matches the
hub-size distribution observed in the real catfish_226 cohort, with the
FULL 23-column output (theta + IBS0 + KING + R0 + R1 + J9 + ...).

Hub design (matches catfish_226 hub_id_1st distribution exactly):
  2 hubs of size 16
  1 hub of size 12
  1 hub of size 9
  2 hubs of size 8
  3 hubs of size 6
  5 hubs of size 5
  9 hubs of size 4
  9 hubs of size 3
  15 hubs of size 2
  21 singletons
  Total: 47 multi-node components + 21 singletons = 68 components, 226 samples

Hub-topology assignment (designed to exercise every classifier branch):
  size 16 (×2): two_parents_with_sibship  (2 parents + 14 offspring)
  size 12: parent_with_sibship             (1 parent + 11 offspring)
  size  9: parent_with_sibship             (1 parent + 8 offspring)
  size  8 (×2): sibship_only               (parents not sampled, 8 full sibs)
  size  6 (×3): mix:
              [0] two_parents_with_sibship (2 parents + 4 offspring)
              [1] parent_with_sibship      (1 parent + 5 offspring)
              [2] sibship_only
  size  5 (×5): mix:
              [0,1] parent_with_sibship    (1 parent + 4 offspring)
              [2,3,4] sibship_only
  size  4 (×9): mix:
              [0..3] parent_with_sibship   (1 parent + 3 offspring)
              [4..8] sibship_only
  size  3 (×9): mix:
              [0..3] parent_with_sibship   (1 parent + 2 offspring)
              [4..8] sibship_only          (3 full sibs)
  size  2 (×15): mix:
              [0..7] po_dyad_only          (parent + offspring, but only 1 edge)
              [8..14] sibship_only         (2 full sibs)

Plus duplicate hubs:
  Insert 2 duplicate pairs as small extra components — but we already have
  size-2 hubs, so we'll demote 2 of the size-2 sibship_only hubs to duplicate_pair.

Result counts after this design:
  duplicate_pair: 2
  two_parents_with_sibship: 3 (sizes 16,16,6)
  parent_with_sibship: 1+1+1+2+4+4 = 13
  sibship_only: 2+1+3+5+5+(15-8-2) = 21  (the size-2 sibship_only ones)
  po_dyad_only: 8

Output:
  <outdir>/relatedness.res (full 23-col)
  <outdir>/samples.txt
  <outdir>/expected_truth.tsv (per-sample expected hub_type and role)
  <outdir>/hub_design.tsv     (per-hub design summary for the test runner)
"""

import argparse
import math
import random
from pathlib import Path

# Hub design — list of (n_members, designed_topology) tuples.
# This is the spec; the test runner reads it and asserts against script output.
HUB_DESIGN = []
def add(n, topo, count=1):
    for _ in range(count):
        HUB_DESIGN.append((n, topo))

add(16, "two_parents_with_sibship", 2)
add(12, "parent_with_sibship")
add( 9, "parent_with_sibship")
add( 8, "sibship_only", 2)
add( 6, "two_parents_with_sibship")
add( 6, "parent_with_sibship")
add( 6, "sibship_only")
add( 5, "parent_with_sibship", 2)
add( 5, "sibship_only", 3)
add( 4, "parent_with_sibship", 4)
add( 4, "sibship_only", 5)
add( 3, "parent_with_sibship", 4)
add( 3, "sibship_only", 5)
add( 2, "po_dyad_only", 8)
add( 2, "sibship_only", 5)
add( 2, "duplicate_pair", 2)
# Total: 47 multi-node hubs.

# Plus 21 singletons (no in-hub edges).
N_SINGLETONS = 21


def synth_pair_values(rel):
    """Generate a row of ngsRelate columns consistent with `rel`.

    Adds small random jitter so values aren't identical across pairs.
    """
    j = lambda lo, hi: random.uniform(lo, hi)

    if rel == "PO":
        return dict(
            theta=j(0.235, 0.260), IBS0=j(0.0001, 0.0040),
            KING=j(0.230, 0.255), R0=j(0.005, 0.045), R1=j(0.78, 0.92),
            J7=j(0.02, 0.06), J8=j(0.02, 0.06), J9=j(0.02, 0.08),
            rab=j(0.46, 0.52), nSites=int(j(80000, 120000)))
    if rel == "FS":
        return dict(
            theta=j(0.235, 0.275), IBS0=j(0.020, 0.080),
            KING=j(0.230, 0.265), R0=j(0.13, 0.24), R1=j(0.55, 0.70),
            J7=j(0.04, 0.10), J8=j(0.04, 0.10), J9=j(0.18, 0.30),
            rab=j(0.46, 0.54), nSites=int(j(80000, 120000)))
    if rel == "DUP":
        return dict(
            theta=j(0.480, 0.499), IBS0=j(0.00001, 0.00050),
            KING=j(0.480, 0.499), R0=j(0.005, 0.020), R1=j(0.005, 0.020),
            J7=j(0.005, 0.020), J8=j(0.005, 0.020), J9=j(0.93, 0.99),
            rab=j(0.96, 1.00), nSites=int(j(80000, 120000)))
    if rel == "2nd":
        return dict(
            theta=j(0.10, 0.16), IBS0=j(0.06, 0.13),
            KING=j(0.09, 0.16), R0=j(0.28, 0.40), R1=j(0.40, 0.55),
            J7=j(0.04, 0.10), J8=j(0.04, 0.10), J9=j(0.06, 0.15),
            rab=j(0.18, 0.32), nSites=int(j(80000, 120000)))
    if rel == "3rd":
        return dict(
            theta=j(0.05, 0.085), IBS0=j(0.10, 0.16),
            KING=j(0.04, 0.08), R0=j(0.38, 0.48), R1=j(0.28, 0.40),
            J7=j(0.02, 0.06), J8=j(0.02, 0.06), J9=j(0.02, 0.08),
            rab=j(0.10, 0.18), nSites=int(j(80000, 120000)))
    # unrelated
    return dict(
        theta=j(0.000, 0.040), IBS0=j(0.13, 0.22),
        KING=j(-0.02, 0.04), R0=j(0.45, 0.55), R1=j(0.45, 0.55),
        J7=j(0.01, 0.04), J8=j(0.01, 0.04), J9=j(0.01, 0.06),
        rab=j(0.00, 0.10), nSites=int(j(80000, 120000)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Build sample IDs and assign them to hubs ----
    samples = []        # list of sample_ids in order
    hub_assign = []     # parallel list: which hub index (0..46) or -1 for singleton
    hub_role_design = []  # parallel list: designed role within hub ('parent_a','parent_b','offspring','fs','dup', or '_isolated')

    next_id = 1
    def mkid():
        nonlocal next_id
        s = f"CGA{next_id:03d}"
        next_id += 1
        return s

    # Multi-node hubs
    for hub_idx, (n, topo) in enumerate(HUB_DESIGN):
        if topo == "two_parents_with_sibship":
            n_off = n - 2
            assert n_off >= 2
            ids = [mkid() for _ in range(n)]
            for k, sid in enumerate(ids):
                samples.append(sid)
                hub_assign.append(hub_idx)
                if k == 0:    hub_role_design.append("parent_a")
                elif k == 1:  hub_role_design.append("parent_b")
                else:         hub_role_design.append("offspring")
        elif topo == "parent_with_sibship":
            ids = [mkid() for _ in range(n)]
            for k, sid in enumerate(ids):
                samples.append(sid)
                hub_assign.append(hub_idx)
                if k == 0:    hub_role_design.append("parent")
                else:         hub_role_design.append("offspring")
        elif topo == "sibship_only":
            ids = [mkid() for _ in range(n)]
            for sid in ids:
                samples.append(sid)
                hub_assign.append(hub_idx)
                hub_role_design.append("fs")
        elif topo == "po_dyad_only":
            assert n == 2
            ids = [mkid() for _ in range(2)]
            for sid in ids:
                samples.append(sid)
                hub_assign.append(hub_idx)
                hub_role_design.append("po_dyad")
        elif topo == "duplicate_pair":
            assert n == 2
            ids = [mkid() for _ in range(2)]
            for sid in ids:
                samples.append(sid)
                hub_assign.append(hub_idx)
                hub_role_design.append("dup")
        else:
            raise RuntimeError(f"unknown topo {topo}")

    # Singletons
    for _ in range(N_SINGLETONS):
        samples.append(mkid())
        hub_assign.append(-1)
        hub_role_design.append("_isolated")

    n_total = len(samples)
    print(f"[gen] total samples: {n_total}")
    print(f"[gen] multi-node hubs: {len(HUB_DESIGN)}, singletons: {N_SINGLETONS}")

    # ---- Build per-pair relationship matrix ----
    # Default: unrelated. Then overwrite for in-hub pairs based on roles.
    rels = {}  # (a_idx, b_idx) -> 'PO'|'FS'|'DUP'|'unrelated' (a < b)

    # Sprinkle in a few 2nd/3rd-degree edges as background structure
    # (to confirm they classify correctly and don't pollute hubs).
    # We'll add 50 random 2nd-degree edges and 100 random 3rd-degree edges
    # between samples in DIFFERENT hubs.
    all_pairs_list = [(i, j) for i in range(n_total) for j in range(i+1, n_total)]
    random.shuffle(all_pairs_list)

    # First mark in-hub edges.
    for hub_idx in range(len(HUB_DESIGN)):
        members = [i for i, h in enumerate(hub_assign) if h == hub_idx]
        topo = HUB_DESIGN[hub_idx][1]
        roles = {i: hub_role_design[i] for i in members}

        for i, mi in enumerate(members):
            for mj in members[i+1:]:
                key = (min(mi, mj), max(mi, mj))
                ri, rj = roles[mi], roles[mj]

                if topo == "two_parents_with_sibship":
                    # parent_a-parent_b: unrelated
                    # parent-offspring: PO
                    # offspring-offspring: FS
                    if {ri, rj} == {"parent_a", "parent_b"}:
                        rels[key] = "unrelated"
                    elif "offspring" in {ri, rj} and ("parent_a" in {ri, rj} or "parent_b" in {ri, rj}):
                        rels[key] = "PO"
                    elif ri == "offspring" and rj == "offspring":
                        rels[key] = "FS"
                elif topo == "parent_with_sibship":
                    if {ri, rj} == {"parent", "offspring"}:
                        rels[key] = "PO"
                    elif ri == "offspring" and rj == "offspring":
                        rels[key] = "FS"
                elif topo == "sibship_only":
                    rels[key] = "FS"
                elif topo == "po_dyad_only":
                    rels[key] = "PO"
                elif topo == "duplicate_pair":
                    rels[key] = "DUP"

    # Add background 2nd / 3rd-degree edges between unmarked pairs in different hubs
    n_2nd = 50
    n_3rd = 100
    placed_2nd = 0
    placed_3rd = 0
    for (i, j) in all_pairs_list:
        if (i, j) in rels:
            continue
        if hub_assign[i] != -1 and hub_assign[j] != -1 and hub_assign[i] == hub_assign[j]:
            continue
        if placed_2nd < n_2nd:
            rels[(i, j)] = "2nd"
            placed_2nd += 1
        elif placed_3rd < n_3rd:
            rels[(i, j)] = "3rd"
            placed_3rd += 1
        else:
            break

    # Default everything else to unrelated.
    # (we don't store unrelated explicitly; the writer falls through.)

    # ---- Write samples.txt ----
    with open(outdir / "samples.txt", "w") as fh:
        for s in samples:
            fh.write(s + "\n")

    # ---- Write the .res with the full 23-col header ----
    cols = ["a", "b", "nSites", "J7", "J8", "J9", "rab", "Fa", "Fb",
            "theta", "inbreed_a", "inbreed_b", "2of3_IDB", "FDiff",
            "loglh", "nIter", "coverage", "IBS0", "IBS1", "IBS2",
            "R0", "R1", "KING"]
    with open(outdir / "relatedness.res", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for (i, j) in all_pairs_list:
            a, b = (i, j) if i < j else (j, i)
            rel = rels.get((a, b), "unrelated")
            v = synth_pair_values(rel)
            row = {
                "a": a, "b": b, "nSites": v["nSites"],
                "J7": v["J7"], "J8": v["J8"], "J9": v["J9"],
                "rab": v["rab"], "Fa": 0.01, "Fb": 0.01,
                "theta": v["theta"], "inbreed_a": 0.01, "inbreed_b": 0.01,
                "2of3_IDB": 0.5, "FDiff": 0.0,
                "loglh": -1000.0, "nIter": 20, "coverage": 9.0,
                "IBS0": v["IBS0"], "IBS1": 0.5, "IBS2": max(0.0, 0.5 - v["IBS0"]),
                "R0": v["R0"], "R1": v["R1"], "KING": v["KING"],
            }
            fh.write("\t".join(f"{row[c]:.6f}" if isinstance(row[c], float) else str(row[c])
                               for c in cols) + "\n")

    # ---- Write expected_truth.tsv ----
    # Translate per-role design into expected (hub_type, role) for the script's output.
    role_translate_blind = {
        # parent_with_sibship of size>=4 (3+ offspring) → forced_parent
        ("parent_with_sibship_3plus_off", "parent"): "forced_parent",
        ("parent_with_sibship_2_off",     "parent"): "likely_parent",  # size-3 hubs
        ("parent_with_sibship_3plus_off", "offspring"): "possible_offspring",
        ("parent_with_sibship_2_off",     "offspring"): "possible_offspring",
        ("two_parents_with_sibship_3plus_off", "parent_a"): "parent_a",
        ("two_parents_with_sibship_3plus_off", "parent_b"): "parent_b",
        ("two_parents_with_sibship_2plus_off", "parent_a"): "parent_a",
        ("two_parents_with_sibship_2plus_off", "parent_b"): "parent_b",
        ("two_parents_with_sibship_3plus_off", "offspring"): "possible_offspring",
        ("two_parents_with_sibship_2plus_off", "offspring"): "possible_offspring",
        ("sibship_only", "fs"): "possible_full_sib",
        ("po_dyad_only", "po_dyad"): "ambiguous_first_degree_PO",
        ("duplicate_pair", "dup"): "duplicate",
        ("isolated", "_isolated"): "isolated",
    }

    with open(outdir / "expected_truth.tsv", "w") as fh:
        fh.write("sample_id\texpected_hub_type\texpected_role\n")
        for sid_idx, sid in enumerate(samples):
            h = hub_assign[sid_idx]
            r = hub_role_design[sid_idx]
            if h == -1:
                fh.write(f"{sid}\tisolated\tisolated\n")
                continue
            n_h, topo = HUB_DESIGN[h]
            # Resolve sub-topology variant for parent rules
            if topo == "parent_with_sibship":
                n_off = n_h - 1
                key = ("parent_with_sibship_3plus_off" if n_off >= 3 else "parent_with_sibship_2_off", r)
                expected_role = role_translate_blind[key]
                expected_hub = "parent_with_sibship"
            elif topo == "two_parents_with_sibship":
                n_off = n_h - 2
                key = ("two_parents_with_sibship_3plus_off" if n_off >= 3 else "two_parents_with_sibship_2plus_off", r)
                expected_role = role_translate_blind[key]
                expected_hub = "two_parents_with_sibship"
            else:
                expected_role = role_translate_blind[(topo, r)]
                expected_hub = topo
            fh.write(f"{sid}\t{expected_hub}\t{expected_role}\n")

    # ---- Write hub_design.tsv summary ----
    with open(outdir / "hub_design.tsv", "w") as fh:
        fh.write("designed_hub_idx\tn_members\tdesigned_topology\n")
        for hub_idx, (n, topo) in enumerate(HUB_DESIGN):
            fh.write(f"{hub_idx}\t{n}\t{topo}\n")

    # Summary of design
    from collections import Counter
    topo_counts = Counter(t for _, t in HUB_DESIGN)
    print(f"[gen] designed hub-type distribution:")
    for k, v in sorted(topo_counts.items()):
        print(f"  {k:30s} {v}")
    print(f"[gen] wrote {outdir}/relatedness.res ({len(all_pairs_list)} pairs)")
    print(f"[gen] wrote {outdir}/samples.txt ({n_total} samples)")
    print(f"[gen] wrote {outdir}/expected_truth.tsv")
    print(f"[gen] wrote {outdir}/hub_design.tsv")


if __name__ == "__main__":
    main()
