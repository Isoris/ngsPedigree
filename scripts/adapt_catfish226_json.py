#!/usr/bin/env python3
"""
Adapt catfish_226_relatedness.json (schema_v2, theta-only) into the
.res + samples.txt shape that STEP_PED_01_annotate_relationships.py expects.

The source JSON carries only theta — no IBS0, no KING, no R0/R1, no J9.
This adapter writes IBS0 as NaN so the classifier triggers its
'ambiguous_first_degree' branch (instead of guessing PO vs FS), per the
May 10 spec correction.

Outputs:
  <outdir>/relatedness.res     — TSV, columns: a, b, nSites, theta, IBS0
  <outdir>/samples.txt         — one sample ID per line
"""

import argparse
import json
import math
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--json", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.json) as fh:
        d = json.load(fh)
    r = d["relatedness"]
    samples = r["samples"]
    n = r["n_samples"]
    a_arr = r["pairs"]["a"]
    b_arr = r["pairs"]["b"]
    t_arr = r["pairs"]["theta"]

    # samples.txt
    with open(outdir / "samples.txt", "w") as fh:
        for s in samples:
            fh.write(s + "\n")

    # relatedness.res — minimal column set; classifier handles missing IBS0
    # by emitting 'ambiguous_first_degree' for first-degree pairs.
    with open(outdir / "relatedness.res", "w") as fh:
        fh.write("a\tb\tnSites\ttheta\tIBS0\n")
        for a, b, t in zip(a_arr, b_arr, t_arr):
            # Write IBS0 as 'NA' so pandas reads it as NaN.
            fh.write(f"{a}\t{b}\t100000\t{t}\tNA\n")

    print(f"Wrote {outdir}/samples.txt ({n} samples)")
    print(f"Wrote {outdir}/relatedness.res ({len(t_arr)} pairs)")
    print()
    print(f"Note: source JSON is theta-only. IBS0 written as NA.")
    print(f"      First-degree pairs will be classified 'ambiguous_first_degree'.")
    print(f"      To get PO vs FS, re-run ngsRelate to get the full 23-column .res.")


if __name__ == "__main__":
    main()
