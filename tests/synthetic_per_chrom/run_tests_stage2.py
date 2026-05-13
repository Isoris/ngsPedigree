#!/usr/bin/env python3
"""Test runner for Stage 2 — synthetic 30-chromosome fixture."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S2_SCRIPT = HERE.parent.parent / "scripts" / "STEP_PED_02_per_chromosome_qc.py"
GEN = HERE / "generate_per_chrom_fixture.py"

assertions = []
def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    assertions.append((status, name, detail))
    if not cond:
        print(f"  FAIL  {name}  {detail}")


def main():
    fixture_dir = HERE / "_fixture"
    if fixture_dir.exists():
        shutil.rmtree(fixture_dir)
    fixture_dir.mkdir()

    print("[1] generating per-chromosome fixture...")
    subprocess.run([sys.executable, str(GEN), "--outdir", str(fixture_dir)],
                   check=True, capture_output=True)

    print("[2] running Stage 2...")
    s2_out = fixture_dir / "stage2_out"
    if s2_out.exists():
        shutil.rmtree(s2_out)
    res = subprocess.run([sys.executable, str(S2_SCRIPT),
                          "--per-chrom-dir", str(fixture_dir / "per_chrom"),
                          "--stage1-outdir", str(fixture_dir / "stage1"),
                          "--samples", str(fixture_dir / "samples.txt"),
                          "--outdir", str(s2_out)],
                         capture_output=True, text=True)
    if res.returncode != 0:
        print("STDOUT:", res.stdout)
        print("STDERR:", res.stderr)
        sys.exit(2)

    # Load outputs
    extended = pd.read_csv(s2_out / "pairwise_relationship_classification.tsv", sep="\t")
    qc = pd.read_csv(s2_out / "per_chromosome_qc_flags.tsv", sep="\t")
    summary = pd.read_csv(s2_out / "per_chromosome_summary.tsv", sep="\t")
    envelope = json.loads((s2_out / "ngspedigree_run_envelope.json").read_text())
    truth = pd.read_csv(fixture_dir / "expected_disagreements.tsv", sep="\t")

    print("\n[3] checking envelope...")
    check("env_30_chroms", envelope["n_chromosomes_processed"] == 30,
          f"got {envelope['n_chromosomes_processed']}")
    check("env_stage=stage2", envelope.get("stage") == "stage2")

    # ---- Per-chromosome summary checks ----
    print("[4] checking per-chromosome summary...")
    check("summary_has_30_rows", len(summary) == 30, f"got {len(summary)}")
    # LG29, LG30 should have low mean_n_sites
    for low_chrom in ("LG29", "LG30"):
        row = summary[summary["chrom"] == low_chrom]
        check(f"{low_chrom}_low_mean_sites",
              not row.empty and row["mean_n_sites"].iloc[0] < 1000,
              f"got mean_n_sites={row['mean_n_sites'].iloc[0] if not row.empty else 'missing'}")

    # ---- Disagreement detection ----
    print("[5] checking disagreement detection on LG14...")
    expected_pairs = set((r["sample_a"], r["sample_b"]) for _, r in truth.iterrows()
                         if r["chrom"] == "LG14")

    # The QC table should have a row for each expected disagreement.
    lg14_flags = qc[(qc["chrom"] == "LG14") & (qc["flag"] == "disagreement")]
    detected_pairs = set((r["sample_a"], r["sample_b"]) for _, r in lg14_flags.iterrows())

    check("LG14_disagreements_detected",
          expected_pairs.issubset(detected_pairs),
          f"expected {expected_pairs}, got {detected_pairs}")

    # ---- low_data flag detection ----
    print("[6] checking low_data flag for LG29, LG30...")
    for low_chrom in ("LG29", "LG30"):
        low_flags = qc[(qc["chrom"] == low_chrom) & (qc["flag"] == "low_data")]
        # Many pairs should be flagged low_data on these chromosomes
        check(f"{low_chrom}_low_data_flagged",
              len(low_flags) > 100,  # we have >25k pairs; expect lots flagged
              f"got {len(low_flags)} flags on {low_chrom}")

    # ---- Pair-level review flag ----
    print("[7] checking pair-level review flags...")
    # 5 PO pairs had their LG14 flipped. n_chrom_disagreements for those 5 is 1.
    # With default --disagreement-threshold=3, none should be flagged for review.
    # Verify the disagreement count was 1 for the affected pairs.
    for sa, sb in expected_pairs:
        row = extended[((extended["sample_a"] == sa) & (extended["sample_b"] == sb)) |
                       ((extended["sample_a"] == sb) & (extended["sample_b"] == sa))]
        if row.empty:
            check(f"pair_{sa}_{sb}_in_table", False, "missing")
            continue
        n_dis = row["n_chrom_disagreements"].iloc[0]
        check(f"pair_{sa}_{sb}_n_disagreements=1", n_dis == 1,
              f"got {n_dis}")

    # ---- Genome-wide table preserved ----
    print("[8] checking Stage 1 columns preserved...")
    for col in ["edge_class", "confidence", "theta", "IBS0", "KING"]:
        check(f"stage1_col_{col}_present", col in extended.columns)

    # ---- Per-chromosome columns present ----
    chrom_cols = [c for c in extended.columns if c.startswith("edge_class_LG")]
    check("30_chrom_class_cols", len(chrom_cols) == 30, f"got {len(chrom_cols)}")

    # ---- Summary ----
    print()
    n_pass = sum(1 for s,_,_ in assertions if s == "PASS")
    n_fail = sum(1 for s,_,_ in assertions if s == "FAIL")
    print("=" * 70)
    print(f"Stage 2 per-chromosome QC: {n_pass} PASS, {n_fail} FAIL  (of {len(assertions)})")
    print("=" * 70)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
