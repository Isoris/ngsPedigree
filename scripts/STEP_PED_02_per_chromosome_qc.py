#!/usr/bin/env python3
"""
ngsPedigree — STEP_PED_02_per_chromosome_qc.py
================================================

Stage 2 of ngsPedigree: per-chromosome QC of the relationship classifications
produced by Stage 1.

What this does
--------------
For each per-chromosome ngsRelate .res file produced by
`STEP_A07b_relatedness_per_chrom.sh`, run the SAME edge classifier as Stage 1
(imported from STEP_PED_01) on that chromosome's data alone. Then add one
column per chromosome to the existing edge table from Stage 1, and flag
pairs where local class disagrees with genome-wide class.

What this does NOT do
---------------------
- Does NOT re-classify hubs. Genome-wide hub topology from Stage 1 is the
  authoritative classification; Stage 2 only adds per-chromosome QC columns.
- Does NOT compute haplotype phase. That's Stage 3.
- Does NOT call recombination, gene conversion, or double crossovers.
  That's ngsTracts.

Inputs
------
  --per-chrom-dir PATH    : directory containing per-chromosome .res files,
                             one per chromosome. File naming convention:
                             <chrom>.res (e.g. LG01.res, LG02.res, ...)
  --stage1-outdir PATH    : Stage 1 output directory (contains
                             pairwise_relationship_classification.tsv,
                             family_hub_roster.tsv, ngspedigree_run_envelope.json)
  --samples PATH          : sample sidecar (same one Stage 1 used)
  --outdir PATH           : where to write Stage 2 outputs
  [--low-data-threshold]  : min n_sites per chromosome to attempt classification
                             (default: 1000)
  [--ibs0-po-max-perchrom]: per-chromosome IBS0 threshold for PO call
                             (default: 0.008 — slightly looser than the
                             genome-wide 0.005 because per-chromosome IBS0
                             has more sampling noise)
  [--disagreement-threshold]: min number of chromosomes where local class
                             must disagree with genome-wide for the pair to
                             be flagged (default: 3)

Outputs
-------
  pairwise_relationship_classification.tsv  : Stage 1's table extended with
                                              one column per chromosome
                                              (edge_class_LG01, ...) plus
                                              n_chrom_disagreements
  per_chromosome_qc_flags.tsv               : per-(pair, chromosome) flags
                                              for chromosomes flagged as
                                              suspicious or low_data
  per_chromosome_summary.tsv                : per-chromosome summary stats
                                              (n_pairs_classified, mean
                                              n_sites, etc.)
  ngspedigree_run_envelope.json             : updated envelope including
                                              Stage 2 metadata

Usage
-----
  python STEP_PED_02_per_chromosome_qc.py \\
      --per-chrom-dir data/cohort/relatedness/ngsrelate/cohort_226_full_v1/per_chrom/ \\
      --stage1-outdir data/cohort/relatedness/ngsrelate/cohort_226_full_v1/ngspedigree/ \\
      --samples data/cohort/relatedness/ngsrelate/cohort_226_full_v1/samples.txt \\
      --outdir data/cohort/relatedness/ngsrelate/cohort_226_full_v1/ngspedigree_stage2/

Author: ngsPedigree v0.2, 2026-05-10
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.stderr.write("ERROR: pandas is required. Install with: pip install pandas\n")
    sys.exit(2)

# Import Stage 1's classifier (this is the load-bearing reuse).
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from STEP_PED_01_annotate_relationships import (
    classify_edge,
    load_samples,
    load_res,
    REQUIRED_COLS,
)

SCRIPT_VERSION = "v0.2.0"
SCRIPT_NAME    = "STEP_PED_02_per_chromosome_qc"


def discover_per_chrom_files(per_chrom_dir):
    """Find all .res files in the per-chrom dir. Returns dict {chrom_name: path}."""
    per_chrom_dir = Path(per_chrom_dir)
    if not per_chrom_dir.is_dir():
        sys.stderr.write(f"FATAL: --per-chrom-dir {per_chrom_dir} is not a directory\n")
        sys.exit(3)
    files = sorted(per_chrom_dir.glob("*.res"))
    if not files:
        sys.stderr.write(f"FATAL: no .res files found in {per_chrom_dir}\n")
        sys.exit(3)
    chrom_files = {}
    for f in files:
        chrom_name = f.stem  # filename without extension, e.g. "LG01"
        chrom_files[chrom_name] = f
    return chrom_files


def classify_chrom(chrom_name, res_path, samples, thresholds, low_data_threshold):
    """Run the Stage 1 classifier on one chromosome's .res.

    Returns DataFrame with columns: a, b, sample_a, sample_b,
    edge_class_<chrom>, n_sites_<chrom>, IBS0_<chrom>.
    """
    df = load_res(res_path)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.stderr.write(f"  [{chrom_name}] WARN: missing required columns {missing}; skipping\n")
        return None

    # Map index to sample id
    df["sample_a"] = df["a"].map(lambda i: samples[i])
    df["sample_b"] = df["b"].map(lambda i: samples[i])

    # Run the classifier per-row using Stage 1's function
    classifications = df.apply(lambda r: classify_edge(r.to_dict(), thresholds), axis=1)
    df[f"edge_class_{chrom_name}"] = classifications.map(lambda c: c["edge_class"])

    # Mark low-data rows explicitly
    if "nSites" in df.columns:
        low_mask = df["nSites"].fillna(0) < low_data_threshold
        df.loc[low_mask, f"edge_class_{chrom_name}"] = "low_data"

    # Keep only the columns we need to merge back
    cols = ["sample_a", "sample_b", f"edge_class_{chrom_name}"]
    if "nSites" in df.columns:
        cols.append("nSites")
        df = df.rename(columns={"nSites": f"n_sites_{chrom_name}"})
        cols[-1] = f"n_sites_{chrom_name}"
    if "IBS0" in df.columns:
        df = df.rename(columns={"IBS0": f"IBS0_{chrom_name}"})
        cols.append(f"IBS0_{chrom_name}")
    return df[cols]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--per-chrom-dir", required=True)
    p.add_argument("--stage1-outdir", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--low-data-threshold", type=int, default=1000)
    p.add_argument("--ibs0-po-max-perchrom", type=float, default=0.008)
    p.add_argument("--disagreement-threshold", type=int, default=3)
    p.add_argument("--theta-first", type=float, default=0.177)
    p.add_argument("--theta-second", type=float, default=0.0884)
    p.add_argument("--theta-third", type=float, default=0.0442)
    p.add_argument("--theta-dup-min", type=float, default=0.45)
    args = p.parse_args()

    # Per-chromosome thresholds. Same as Stage 1 except IBS0_PO_MAX is looser.
    thresholds_perchrom = dict(
        theta_first=args.theta_first,
        theta_second=args.theta_second,
        theta_third=args.theta_third,
        theta_dup_min=args.theta_dup_min,
        ibs0_po_max=args.ibs0_po_max_perchrom,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load Stage 1 outputs ----
    s1_dir = Path(args.stage1_outdir)
    s1_edges_path = s1_dir / "pairwise_relationship_classification.tsv"
    s1_envelope_path = s1_dir / "ngspedigree_run_envelope.json"
    if not s1_edges_path.exists():
        sys.stderr.write(f"FATAL: Stage 1 edges file missing: {s1_edges_path}\n")
        sys.exit(4)
    if not s1_envelope_path.exists():
        sys.stderr.write(f"FATAL: Stage 1 envelope missing: {s1_envelope_path}\n")
        sys.exit(4)

    sys.stderr.write(f"[ngsPedigree-S2] Loading Stage 1 edges from {s1_edges_path}\n")
    s1_edges = pd.read_csv(s1_edges_path, sep="\t")
    sys.stderr.write(f"[ngsPedigree-S2]   {len(s1_edges)} pairs from Stage 1\n")

    with open(s1_envelope_path) as fh:
        s1_envelope = json.load(fh)

    # ---- Load samples ----
    samples = load_samples(args.samples)
    sys.stderr.write(f"[ngsPedigree-S2] {len(samples)} samples\n")

    # ---- Discover per-chromosome .res files ----
    chrom_files = discover_per_chrom_files(args.per_chrom_dir)
    sys.stderr.write(f"[ngsPedigree-S2] Found {len(chrom_files)} per-chromosome .res files\n")
    for c in sorted(chrom_files):
        sys.stderr.write(f"  {c}\n")

    # ---- Classify each chromosome and collect results ----
    sys.stderr.write(f"\n[ngsPedigree-S2] Classifying each chromosome...\n")
    per_chrom_dfs = []
    per_chrom_summary = []
    for chrom_name in sorted(chrom_files):
        sys.stderr.write(f"[ngsPedigree-S2] {chrom_name}...\n")
        chrom_df = classify_chrom(
            chrom_name, chrom_files[chrom_name], samples,
            thresholds_perchrom, args.low_data_threshold)
        if chrom_df is None:
            continue
        per_chrom_dfs.append(chrom_df)
        # summary
        col = f"edge_class_{chrom_name}"
        n_sites_col = f"n_sites_{chrom_name}"
        summary = {
            "chrom": chrom_name,
            "n_pairs": len(chrom_df),
            "n_low_data": (chrom_df[col] == "low_data").sum(),
            "n_PO": (chrom_df[col] == "parent_offspring").sum(),
            "n_FS": (chrom_df[col] == "full_sibling").sum(),
            "n_dup": (chrom_df[col] == "duplicate_or_clone").sum(),
            "n_unrelated": (chrom_df[col] == "unrelated").sum(),
        }
        if n_sites_col in chrom_df.columns:
            summary["mean_n_sites"] = float(chrom_df[n_sites_col].mean())
            summary["min_n_sites"] = int(chrom_df[n_sites_col].min())
        per_chrom_summary.append(summary)

    # ---- Merge all per-chromosome columns into the Stage 1 table ----
    sys.stderr.write(f"\n[ngsPedigree-S2] Merging per-chromosome columns into Stage 1 table...\n")
    extended = s1_edges.copy()
    for chrom_df in per_chrom_dfs:
        extended = extended.merge(chrom_df, on=["sample_a", "sample_b"], how="left")

    # ---- Compute disagreement count per pair ----
    chrom_class_cols = [f"edge_class_{c}" for c in sorted(chrom_files)
                        if f"edge_class_{c}" in extended.columns]
    sys.stderr.write(f"[ngsPedigree-S2] Computing disagreement counts across {len(chrom_class_cols)} chromosomes\n")

    def count_disagreements(row):
        gw = row["edge_class"]
        n_disagree = 0
        n_compared = 0
        for c in chrom_class_cols:
            local = row.get(c)
            if pd.isna(local) or local == "low_data":
                continue
            n_compared += 1
            if local != gw:
                n_disagree += 1
        return pd.Series({"n_chrom_compared": n_compared,
                          "n_chrom_disagreements": n_disagree,
                          "frac_disagreement": (n_disagree / n_compared) if n_compared else 0.0})

    disagree_df = extended.apply(count_disagreements, axis=1)
    extended = pd.concat([extended, disagree_df], axis=1)

    # ---- Per-chromosome QC flags ----
    qc_rows = []
    for chrom_name in sorted(chrom_files):
        col = f"edge_class_{chrom_name}"
        if col not in extended.columns:
            continue
        for _, row in extended.iterrows():
            local = row[col]
            if pd.isna(local):
                continue
            if local == "low_data":
                qc_rows.append({
                    "sample_a": row["sample_a"],
                    "sample_b": row["sample_b"],
                    "chrom": chrom_name,
                    "edge_class_genome_wide": row["edge_class"],
                    "edge_class_local": local,
                    "flag": "low_data",
                })
            elif local != row["edge_class"] and row["edge_class"] != "unrelated":
                qc_rows.append({
                    "sample_a": row["sample_a"],
                    "sample_b": row["sample_b"],
                    "chrom": chrom_name,
                    "edge_class_genome_wide": row["edge_class"],
                    "edge_class_local": local,
                    "flag": "disagreement",
                })
    qc_df = pd.DataFrame(qc_rows) if qc_rows else pd.DataFrame(
        columns=["sample_a", "sample_b", "chrom",
                 "edge_class_genome_wide", "edge_class_local", "flag"])

    # ---- Pair-level review flag ----
    extended["pair_review_flag"] = "OK"
    extended.loc[extended["n_chrom_disagreements"] >= args.disagreement_threshold,
                 "pair_review_flag"] = "REVIEW_disagreement_above_threshold"
    n_review = (extended["pair_review_flag"] != "OK").sum()
    sys.stderr.write(f"[ngsPedigree-S2] {n_review} pairs flagged for review "
                     f"(>= {args.disagreement_threshold} chromosomes disagree)\n")

    # ---- Write outputs ----
    out_edges = outdir / "pairwise_relationship_classification.tsv"
    extended.to_csv(out_edges, sep="\t", index=False)
    sys.stderr.write(f"[ngsPedigree-S2] Wrote {out_edges}\n")

    out_qc = outdir / "per_chromosome_qc_flags.tsv"
    qc_df.to_csv(out_qc, sep="\t", index=False)
    sys.stderr.write(f"[ngsPedigree-S2] Wrote {out_qc} ({len(qc_df)} flag rows)\n")

    out_summary = outdir / "per_chromosome_summary.tsv"
    pd.DataFrame(per_chrom_summary).to_csv(out_summary, sep="\t", index=False)
    sys.stderr.write(f"[ngsPedigree-S2] Wrote {out_summary}\n")

    # ---- Updated envelope ----
    envelope = {
        "schema": "ngspedigree_run_envelope_v1",
        "stage": "stage2",
        "produced_by": {
            "tool": SCRIPT_NAME,
            "version": SCRIPT_VERSION,
            "params": {
                "thresholds_perchrom": thresholds_perchrom,
                "low_data_threshold": args.low_data_threshold,
                "disagreement_threshold": args.disagreement_threshold,
            },
        },
        "inputs": {
            "stage1_outdir": str(args.stage1_outdir),
            "per_chrom_dir": str(args.per_chrom_dir),
            "samples_path": str(args.samples),
        },
        "computed_at": datetime.now(timezone.utc).isoformat(),
        "stage1_envelope": s1_envelope,
        "n_chromosomes_processed": len(per_chrom_dfs),
        "n_pairs": int(len(extended)),
        "n_pairs_flagged_for_review": int(n_review),
        "n_qc_flag_rows": int(len(qc_df)),
        "artifacts": {
            "pairwise_relationship_classification": "pairwise_relationship_classification.tsv",
            "per_chromosome_qc_flags": "per_chromosome_qc_flags.tsv",
            "per_chromosome_summary": "per_chromosome_summary.tsv",
        },
        "downstream": {
            "next_step": "STEP_PED_03_inheritance_map.py",
            "consumers": ["STEP_PED_03_inheritance_map.py", "ngsTracts"],
        },
    }
    out_env = outdir / "ngspedigree_run_envelope.json"
    with open(out_env, "w") as fh:
        json.dump(envelope, fh, indent=2, default=str)
    sys.stderr.write(f"[ngsPedigree-S2] Wrote {out_env}\n")

    # Summary
    sys.stderr.write(f"\n[ngsPedigree-S2] ===== SUMMARY =====\n")
    sys.stderr.write(f"  chromosomes processed: {len(per_chrom_dfs)}\n")
    sys.stderr.write(f"  pairs in table: {len(extended)}\n")
    sys.stderr.write(f"  pairs flagged for review: {n_review}\n")
    sys.stderr.write(f"  QC flag rows: {len(qc_df)}\n")
    sys.stderr.write(f"  outputs in: {outdir}\n")


if __name__ == "__main__":
    main()
