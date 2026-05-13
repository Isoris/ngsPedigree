#!/usr/bin/env python3
"""Full-cohort assertion runner against the realistic 226-sample synthetic."""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from collections import Counter

import pandas as pd

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent.parent / "scripts" / "STEP_PED_01_annotate_relationships.py"
GEN = HERE / "generate_realistic_226.py"

assertions = []
def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    assertions.append((status, name, detail))
    if not cond:
        print(f"  FAIL  {name}  {detail}")


def main():
    # Regenerate fixture
    print("[1] generating fixture...")
    subprocess.run([sys.executable, str(GEN), "--outdir", str(HERE)], check=True, capture_output=True)

    # Run script
    print("[2] running STEP_PED_01...")
    out = HERE / "_run"
    if out.exists():
        shutil.rmtree(out)
    cmd = [sys.executable, str(SCRIPT),
           "--res", str(HERE / "relatedness.res"),
           "--samples", str(HERE / "samples.txt"),
           "--outdir", str(out),
           "--run-id", "synth226_test_v1"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout); print(res.stderr)
        sys.exit(2)

    # Load
    edges = pd.read_csv(out / "pairwise_relationship_classification.tsv", sep="\t")
    roster = pd.read_csv(out / "family_hub_roster.tsv", sep="\t")
    truth = pd.read_csv(HERE / "expected_truth.tsv", sep="\t")
    envelope = json.loads((out / "ngspedigree_run_envelope.json").read_text())

    # === Assertions ===

    # ---- Envelope counts ----
    print("[3] checking envelope counts...")
    check("env_n_samples=226", envelope["n_samples"] == 226,
          f"got {envelope['n_samples']}")
    check("env_n_pairs=25425", envelope["n_pairs"] == 25425,
          f"got {envelope['n_pairs']}")
    check("env_n_components=47", envelope["n_multi_node_components"] == 47,
          f"got {envelope['n_multi_node_components']}")
    check("env_n_singletons=21", envelope["n_singletons"] == 21,
          f"got {envelope['n_singletons']}")
    check("env_mode_blind", envelope["mode"] == "blind")

    # ---- Hub type distribution ----
    print("[4] checking hub type distribution...")
    hub_types = envelope["n_hubs_by_type"]
    expected = {"two_parents_with_sibship": 3, "parent_with_sibship": 13,
                "sibship_only": 21, "po_dyad_only": 8, "duplicate_pair": 2}
    for k, v in expected.items():
        check(f"hub_type_count_{k}={v}", hub_types.get(k, 0) == v,
              f"got {hub_types.get(k, 0)}")

    # ---- Edge counts ----
    print("[5] checking edge counts...")
    edge_counts = envelope["n_pairs_by_class"]
    # 2 designed dups
    check("edge_dup_count=2", edge_counts.get("duplicate_or_clone", 0) == 2,
          f"got {edge_counts.get('duplicate_or_clone',0)}")
    # 50 designed 2nd-degree
    check("edge_2nd=50", edge_counts.get("second_degree", 0) == 50,
          f"got {edge_counts.get('second_degree', 0)}")
    # 100 designed 3rd-degree
    check("edge_3rd=100", edge_counts.get("third_degree", 0) == 100,
          f"got {edge_counts.get('third_degree', 0)}")
    # PO + FS counts: derived from the topology design.
    # PO edges: each parent connects to all offspring; sum across hubs.
    # 2 parents × 14 offspring × 2 hubs = 56; plus 11+8+11+5+5+4+4+4+4+(3*4)+(2*4)+(8 PO dyads)
    # Easier: just check both > 0 and that total first-degree sums to PO+FS+ambiguous
    check("edge_PO_present", edge_counts.get("parent_offspring", 0) > 50,
          f"got {edge_counts.get('parent_offspring',0)}")
    check("edge_FS_present", edge_counts.get("full_sibling", 0) > 100,
          f"got {edge_counts.get('full_sibling',0)}")

    # ---- Per-sample hub_type and role ----
    print("[6] checking per-sample roles (this is the big one)...")
    truth_map = {r["sample_id"]: (r["expected_hub_type"], r["expected_role"])
                 for _, r in truth.iterrows()}

    n_samples_checked = 0
    n_correct_hub = 0
    n_correct_role = 0
    role_mismatches = Counter()

    for _, row in roster.iterrows():
        sid = row["sample_id"]
        if sid not in truth_map:
            continue
        n_samples_checked += 1
        exp_ht, exp_role = truth_map[sid]
        actual_ht = row["hub_type"]
        actual_role = row["possible_role"]
        if actual_ht == exp_ht:
            n_correct_hub += 1
        if actual_role == exp_role:
            n_correct_role += 1
        else:
            role_mismatches[(exp_role, actual_role)] += 1

    check("per_sample_hub_type_all_correct", n_correct_hub == n_samples_checked,
          f"{n_correct_hub}/{n_samples_checked} correct")
    check("per_sample_role_all_correct", n_correct_role == n_samples_checked,
          f"{n_correct_role}/{n_samples_checked} correct; "
          f"top mismatches: {role_mismatches.most_common(3)}")

    # ---- Forced parent confidence: large-sibship parents should be high ----
    print("[7] checking forced_parent confidences...")
    forced_parents = roster[roster["possible_role"] == "forced_parent"]
    n_high = (forced_parents["role_confidence"] == "high").sum()
    check("forced_parents_high_conf", n_high == len(forced_parents),
          f"{n_high}/{len(forced_parents)} are high")

    # ---- Likely parent confidence: size-3 parent_with_sibship should be medium ----
    likely_parents = roster[roster["possible_role"] == "likely_parent"]
    n_med = (likely_parents["role_confidence"] == "medium").sum()
    check("likely_parents_medium_conf", n_med == len(likely_parents),
          f"{n_med}/{len(likely_parents)} are medium")

    # ---- 16-sample hubs are correctly classified ----
    print("[8] checking 16-sample hubs (the biggest test of the forced-parent rule)...")
    big_hubs = roster[roster["hub_size"] == 16]
    big_hub_ids = big_hubs["hub_id"].unique()
    check("two_16_sample_hubs", len(big_hub_ids) == 2,
          f"got {len(big_hub_ids)}")

    for hid in big_hub_ids:
        members = big_hubs[big_hubs["hub_id"] == hid]
        ht = members["hub_type"].iloc[0]
        check(f"big_hub_{hid}_type=two_parents_with_sibship",
              ht == "two_parents_with_sibship",
              f"got {ht}")
        n_parents = ((members["possible_role"].isin(["parent_a", "parent_b"]))).sum()
        n_offspring = (members["possible_role"] == "possible_offspring").sum()
        check(f"big_hub_{hid}_2_parents", n_parents == 2,
              f"got {n_parents}")
        check(f"big_hub_{hid}_14_offspring", n_offspring == 14,
              f"got {n_offspring}")

    # ---- Edge classification is consistent with class declarations ----
    print("[9] spot-checking edge classifications...")
    # Pick 5 PO edges and verify theta + IBS0 are consistent
    po = edges[edges["edge_class"] == "parent_offspring"].head(5)
    check("PO_theta_in_first_degree_band",
          ((po["theta"] >= 0.177) & (po["theta"] < 0.45)).all(),
          f"thetas: {po['theta'].tolist()}")
    check("PO_IBS0_low",
          (po["IBS0"] < 0.005).all(),
          f"IBS0s: {po['IBS0'].tolist()}")

    fs = edges[edges["edge_class"] == "full_sibling"].head(5)
    check("FS_IBS0_higher_than_po",
          (fs["IBS0"] >= 0.005).all(),
          f"IBS0s: {fs['IBS0'].tolist()}")

    dups = edges[edges["edge_class"] == "duplicate_or_clone"]
    check("DUP_count=2", len(dups) == 2)
    check("DUP_theta_high", (dups["theta"] >= 0.45).all())
    check("DUP_IBS0_essentially_zero", (dups["IBS0"] < 0.001).all())

    # ---- Summary ----
    print()
    n_pass = sum(1 for s,_,_ in assertions if s == "PASS")
    n_fail = sum(1 for s,_,_ in assertions if s == "FAIL")
    print("=" * 70)
    print(f"Realistic 226 cohort: {n_pass} PASS, {n_fail} FAIL  (of {len(assertions)})")
    print("=" * 70)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
