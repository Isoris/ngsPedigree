#!/usr/bin/env python3
"""
Test runner for STEP_PED_01_annotate_relationships.py.

Generates the synthetic fixture, runs the script, then verifies:
  1. Edge classifications match TRUTH_EDGES.
  2. Hub roster matches expected_truth.tsv (blind mode).
  3. Sex-assisted run promotes the 3 expected labels.

Usage:
  cd ngsPedigree/tests/synthetic_12sample_fixture
  python run_tests.py

Exit code 0 = all pass; non-zero = at least one assertion failed.
"""

import json
import subprocess
import sys
from pathlib import Path
import shutil

import pandas as pd

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent.parent / "scripts" / "STEP_PED_01_annotate_relationships.py"

assertions = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    assertions.append((status, name, detail))
    if cond:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")


def run_step1(outdir, sex_path=None):
    cmd = [sys.executable, str(SCRIPT),
           "--res", str(HERE / "synth.res"),
           "--samples", str(HERE / "samples.txt"),
           "--outdir", str(outdir),
           "--run-id", "synth_test_v1"]
    if sex_path:
        cmd += ["--sex", str(sex_path)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("STDOUT:", res.stdout)
        print("STDERR:", res.stderr)
        raise RuntimeError(f"STEP_PED_01 failed with exit code {res.returncode}")
    return res.stderr


def main():
    print("=" * 70)
    print("ngsPedigree STEP_PED_01 — synthetic fixture test")
    print("=" * 70)

    # 0. Generate fixture
    print("\n[0/4] Generating synthetic fixture...")
    gen_script = HERE / "generate_fixture.py"
    subprocess.run([sys.executable, str(gen_script)], check=True)

    # 1. Blind-mode run
    blind_out = HERE / "_run_blind"
    if blind_out.exists():
        shutil.rmtree(blind_out)
    print("\n[1/4] Running blind-mode pipeline...")
    log = run_step1(blind_out)

    # Load outputs
    edges = pd.read_csv(blind_out / "pairwise_relationship_classification.tsv", sep="\t")
    roster = pd.read_csv(blind_out / "family_hub_roster.tsv", sep="\t")
    with open(blind_out / "ngspedigree_run_envelope.json") as fh:
        envelope = json.load(fh)

    # 2. Verify edge classifications
    print("\n[2/4] Verifying edge classifications...")
    from generate_fixture import TRUTH_EDGES
    for (sa, sb), expected_class in TRUTH_EDGES.items():
        row = edges[((edges["sample_a"] == sa) & (edges["sample_b"] == sb)) |
                    ((edges["sample_a"] == sb) & (edges["sample_b"] == sa))]
        if row.empty:
            check(f"edge_{sa}_{sb}_present", False, "edge missing from output")
            continue
        actual = row.iloc[0]["edge_class"]
        check(f"edge_{sa}_{sb}_class={expected_class}",
              actual == expected_class,
              f"got {actual}")

    # Spot-check a few unrelated edges
    unrelated_edges = edges[edges["edge_class"] == "unrelated"]
    check("edges_unrelated_count_reasonable",
          len(unrelated_edges) > 30,
          f"got {len(unrelated_edges)} (78 pairs total, ~50 expected unrelated)")

    # 3. Verify hub roster against truth
    print("\n[3/4] Verifying hub roster (blind mode)...")
    truth = pd.read_csv(HERE / "expected_truth.tsv", sep="\t")
    truth_map = {r["sample_id"]: (r["expected_hub_type"], r["expected_role"])
                 for _, r in truth.iterrows()}
    for _, row in roster.iterrows():
        sid = row["sample_id"]
        if sid not in truth_map:
            continue
        exp_ht, exp_role = truth_map[sid]
        check(f"hub_type_{sid}={exp_ht}",
              row["hub_type"] == exp_ht,
              f"got {row['hub_type']}")
        check(f"role_{sid}={exp_role}",
              row["possible_role"] == exp_role,
              f"got {row['possible_role']}")

    # Verify forced_parent in H002 (CGA006 has 3 PO edges → forced_parent, high confidence)
    cga006 = roster[roster["sample_id"] == "CGA006"].iloc[0]
    check("CGA006_forced_parent_high_confidence",
          cga006["role_confidence"] == "high",
          f"got {cga006['role_confidence']}")

    # Verify envelope sanity
    check("envelope_mode_blind", envelope["mode"] == "blind")
    check("envelope_n_samples=13", envelope["n_samples"] == 13)
    check("envelope_n_pairs=78", envelope["n_pairs"] == 78)
    check("envelope_has_PO_edges",
          envelope["n_pairs_by_class"].get("parent_offspring", 0) == 9)
    check("envelope_has_FS_edges",
          envelope["n_pairs_by_class"].get("full_sibling", 0) == 7)
    check("envelope_has_dup_edge",
          envelope["n_pairs_by_class"].get("duplicate_or_clone", 0) == 1)

    # 4. Sex-assisted run
    print("\n[4/4] Running sex-assisted pipeline...")
    sex_out = HERE / "_run_sex"
    if sex_out.exists():
        shutil.rmtree(sex_out)
    run_step1(sex_out, sex_path=HERE / "sex.tsv")

    sex_roster = pd.read_csv(sex_out / "family_hub_roster.tsv", sep="\t")
    with open(sex_out / "ngspedigree_run_envelope.json") as fh:
        sex_envelope = json.load(fh)

    cga001 = sex_roster[sex_roster["sample_id"] == "CGA001"].iloc[0]
    cga002 = sex_roster[sex_roster["sample_id"] == "CGA002"].iloc[0]
    cga006 = sex_roster[sex_roster["sample_id"] == "CGA006"].iloc[0]

    check("sex_promote_CGA001_father",
          cga001["possible_role"] == "father",
          f"got {cga001['possible_role']}")
    check("sex_promote_CGA002_mother",
          cga002["possible_role"] == "mother",
          f"got {cga002['possible_role']}")
    check("sex_promote_CGA006_mother",
          cga006["possible_role"] == "mother",
          f"got {cga006['possible_role']}")
    check("sex_envelope_mode_sex_assisted",
          sex_envelope["mode"] == "sex_assisted")
    check("sex_envelope_n_promotions=3",
          sex_envelope["n_label_promotions_via_sex"] == 3,
          f"got {sex_envelope['n_label_promotions_via_sex']}")

    # Summary
    print("\n" + "=" * 70)
    n_pass = sum(1 for s, _, _ in assertions if s == "PASS")
    n_fail = sum(1 for s, _, _ in assertions if s == "FAIL")
    print(f"Results: {n_pass} PASS, {n_fail} FAIL  (of {len(assertions)} assertions)")
    print("=" * 70)

    if n_fail > 0:
        print("\nFailures:")
        for s, name, detail in assertions:
            if s == "FAIL":
                print(f"  {name}  {detail}")
        sys.exit(1)
    else:
        print("\nAll assertions passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
