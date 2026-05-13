#!/usr/bin/env python3
"""
Generate a synthetic per-chromosome fixture for Stage 2 testing.

Design:
  - Reuse the 226-sample synthetic from synthetic_226_realistic/ as Stage 1 input.
  - Generate 30 per-chromosome .res files. For each pair:
      * 27 chromosomes agree with genome-wide
      * 2 chromosomes are "low_data" (n_sites < threshold) for some pairs
      * 1 chromosome carries a deliberate disagreement for a subset of pairs
        (simulating a paralog/mappability artifact on that chromosome)

  - 5 PO pairs get LG14 flipped to "unrelated" locally, simulating a swap
    or a real biological anomaly worth flagging.
  - 2 chromosomes (LG29, LG30) have low n_sites (<1000) for ~all pairs,
    simulating chromosomes too short/thinned to classify reliably.
"""

import argparse
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
S226_DIR = HERE.parent / "synthetic_226_realistic"
S226_GEN = S226_DIR / "generate_realistic_226.py"
S1_SCRIPT = HERE.parent.parent / "scripts" / "STEP_PED_01_annotate_relationships.py"


def synth_pair_values(rel, jitter_seed=None, low_data=False):
    """Same recipe as the 226-sample generator, with optional low-data flag."""
    if jitter_seed is not None:
        random.seed(jitter_seed)
    j = lambda lo, hi: random.uniform(lo, hi)

    base_n = int(j(80000, 120000)) if not low_data else int(j(100, 800))

    if rel == "PO":
        return dict(theta=j(0.235, 0.260), IBS0=j(0.0001, 0.0040),
                    KING=j(0.230, 0.255), R0=j(0.005, 0.045), R1=j(0.78, 0.92),
                    J7=j(0.02, 0.06), J8=j(0.02, 0.06), J9=j(0.02, 0.08),
                    rab=j(0.46, 0.52), nSites=base_n)
    if rel == "FS":
        return dict(theta=j(0.235, 0.275), IBS0=j(0.020, 0.080),
                    KING=j(0.230, 0.265), R0=j(0.13, 0.24), R1=j(0.55, 0.70),
                    J7=j(0.04, 0.10), J8=j(0.04, 0.10), J9=j(0.18, 0.30),
                    rab=j(0.46, 0.54), nSites=base_n)
    if rel == "DUP":
        return dict(theta=j(0.480, 0.499), IBS0=j(0.00001, 0.00050),
                    KING=j(0.480, 0.499), R0=j(0.005, 0.020), R1=j(0.005, 0.020),
                    J7=j(0.005, 0.020), J8=j(0.005, 0.020), J9=j(0.93, 0.99),
                    rab=j(0.96, 1.00), nSites=base_n)
    if rel == "2nd":
        return dict(theta=j(0.10, 0.16), IBS0=j(0.06, 0.13),
                    KING=j(0.09, 0.16), R0=j(0.28, 0.40), R1=j(0.40, 0.55),
                    J7=j(0.04, 0.10), J8=j(0.04, 0.10), J9=j(0.06, 0.15),
                    rab=j(0.18, 0.32), nSites=base_n)
    if rel == "3rd":
        return dict(theta=j(0.05, 0.085), IBS0=j(0.10, 0.16),
                    KING=j(0.04, 0.08), R0=j(0.38, 0.48), R1=j(0.28, 0.40),
                    J7=j(0.02, 0.06), J8=j(0.02, 0.06), J9=j(0.02, 0.08),
                    rab=j(0.10, 0.18), nSites=base_n)
    return dict(theta=j(0.000, 0.040), IBS0=j(0.13, 0.22),
                KING=j(-0.02, 0.04), R0=j(0.45, 0.55), R1=j(0.45, 0.55),
                J7=j(0.01, 0.04), J8=j(0.01, 0.04), J9=j(0.01, 0.06),
                rab=j(0.00, 0.10), nSites=base_n)


def edge_class_to_rel(ec):
    """Translate Stage 1's edge_class output → the synthetic generator's rel string."""
    return {"parent_offspring": "PO", "full_sibling": "FS",
            "duplicate_or_clone": "DUP", "second_degree": "2nd",
            "third_degree": "3rd"}.get(ec, "unrelated")


def write_chrom_res(chrom_path, edges_df, samples, rel_overrides=None,
                    low_data_for_all=False, seed=42):
    """Write a single per-chromosome .res file.

    edges_df: from Stage 1's pairwise_relationship_classification.tsv
    rel_overrides: dict {(sample_a, sample_b): override_rel} for synthetic
                   disagreements. Both directions checked.
    low_data_for_all: if True, every pair gets low n_sites.
    """
    cols = ["a", "b", "nSites", "J7", "J8", "J9", "rab", "Fa", "Fb",
            "theta", "inbreed_a", "inbreed_b", "2of3_IDB", "FDiff",
            "loglh", "nIter", "coverage", "IBS0", "IBS1", "IBS2",
            "R0", "R1", "KING"]

    sample_to_idx = {s: i for i, s in enumerate(samples)}
    rel_overrides = rel_overrides or {}

    with open(chrom_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, row in edges_df.iterrows():
            sa, sb = row["sample_a"], row["sample_b"]
            ec = row["edge_class"]
            # Translate to rel
            rel = edge_class_to_rel(ec)
            # Apply override if present
            override = rel_overrides.get((sa, sb)) or rel_overrides.get((sb, sa))
            if override:
                rel = override
            v = synth_pair_values(rel, jitter_seed=seed + i, low_data=low_data_for_all)
            ai = sample_to_idx[sa]
            bi = sample_to_idx[sb]
            r = {
                "a": ai, "b": bi, "nSites": v["nSites"],
                "J7": v["J7"], "J8": v["J8"], "J9": v["J9"],
                "rab": v["rab"], "Fa": 0.01, "Fb": 0.01,
                "theta": v["theta"], "inbreed_a": 0.01, "inbreed_b": 0.01,
                "2of3_IDB": 0.5, "FDiff": 0.0,
                "loglh": -1000.0, "nIter": 20, "coverage": 9.0,
                "IBS0": v["IBS0"], "IBS1": 0.5, "IBS2": max(0.0, 0.5 - v["IBS0"]),
                "R0": v["R0"], "R1": v["R1"], "KING": v["KING"],
            }
            fh.write("\t".join(f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--n-chroms", type=int, default=30)
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    s1_dir = out / "stage1"
    s1_dir.mkdir(exist_ok=True)
    per_chrom_dir = out / "per_chrom"
    per_chrom_dir.mkdir(exist_ok=True)

    # 1. Make sure 226-sample fixture exists, then run Stage 1 to get the
    #    edge table that per-chromosome generation uses as ground truth.
    print("[1] generating 226-sample fixture...")
    subprocess.run([sys.executable, str(S226_GEN), "--outdir", str(S226_DIR)],
                   check=True, capture_output=True)

    print("[2] running Stage 1 on 226-sample fixture...")
    if (s1_dir / "pairwise_relationship_classification.tsv").exists():
        shutil.rmtree(s1_dir)
        s1_dir.mkdir()
    subprocess.run([sys.executable, str(S1_SCRIPT),
                    "--res", str(S226_DIR / "relatedness.res"),
                    "--samples", str(S226_DIR / "samples.txt"),
                    "--outdir", str(s1_dir),
                    "--run-id", "synth_per_chrom"],
                   check=True, capture_output=True)

    # Load Stage 1's edge table — this is the ground truth genome-wide.
    s1_edges = pd.read_csv(s1_dir / "pairwise_relationship_classification.tsv", sep="\t")
    samples = [s.strip() for s in (S226_DIR / "samples.txt").read_text().splitlines() if s.strip()]
    print(f"[2]   {len(s1_edges)} edges in Stage 1 table")

    # Copy the samples sidecar to outdir
    shutil.copy(S226_DIR / "samples.txt", out / "samples.txt")

    # 3. Generate 30 per-chromosome .res files.
    chrom_names = [f"LG{i:02d}" for i in range(1, args.n_chroms + 1)]
    print(f"[3] generating {len(chrom_names)} per-chromosome .res files...")

    # Pick 5 PO pairs to flip to 'unrelated' on LG14 specifically
    po_edges = s1_edges[s1_edges["edge_class"] == "parent_offspring"].head(5)
    lg14_overrides = {(r["sample_a"], r["sample_b"]): "unrelated"
                      for _, r in po_edges.iterrows()}

    expected_disagreements = {}  # (sa, sb) -> set of chroms where disagreement injected
    for sa, sb in lg14_overrides:
        expected_disagreements.setdefault((sa, sb), set()).add("LG14")

    for chrom in chrom_names:
        chrom_path = per_chrom_dir / f"{chrom}.res"
        if chrom == "LG14":
            write_chrom_res(chrom_path, s1_edges, samples,
                            rel_overrides=lg14_overrides, seed=14)
        elif chrom in ("LG29", "LG30"):
            # Low-data chromosomes
            write_chrom_res(chrom_path, s1_edges, samples,
                            low_data_for_all=True, seed=int(chrom[-2:]))
        else:
            write_chrom_res(chrom_path, s1_edges, samples,
                            seed=int(chrom[-2:]))

    # Write expected disagreements truth file
    truth_rows = []
    for (sa, sb), chroms in expected_disagreements.items():
        for c in chroms:
            truth_rows.append({"sample_a": sa, "sample_b": sb,
                               "chrom": c, "expected_flag": "disagreement"})
    pd.DataFrame(truth_rows).to_csv(out / "expected_disagreements.tsv",
                                    sep="\t", index=False)

    print(f"[gen] Stage 1 dir: {s1_dir}")
    print(f"[gen] per-chrom dir: {per_chrom_dir}")
    print(f"[gen] {args.n_chroms} chromosome .res files")
    print(f"[gen]   LG14: 5 PO pairs flipped to unrelated (expected disagreement flag)")
    print(f"[gen]   LG29, LG30: low n_sites (expected low_data flag)")


if __name__ == "__main__":
    main()
