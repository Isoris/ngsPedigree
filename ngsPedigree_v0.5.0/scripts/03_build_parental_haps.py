#!/usr/bin/env python3
"""
03_build_parental_haps.py — build per-parent hap-1 / hap-2 variant lists
from a (synthetic) joint VCF + parent_phase placeholder, and print a
summary to stdout.

Usage:
  python 03_build_parental_haps.py \
      --vcf          PATH \
      --parent-phase PATH \
      --parents      PID [PID ...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.parental_haps import build_for_parents  # noqa: E402
from hpp.stage3_placeholder import load_parent_phase  # noqa: E402
from hpp.vcf_lite import read_vcf  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--parent-phase", required=True)
    ap.add_argument("--parents", required=True, nargs="+")
    args = ap.parse_args()

    variants = list(read_vcf(args.vcf))
    phase = load_parent_phase(args.parent_phase)
    per_parent = build_for_parents(args.parents, variants, phase)

    for pid, rec in per_parent.items():
        print(f"[{pid}] n_hap1={len(rec.hap1)} "
              f"n_hap2={len(rec.hap2)} "
              f"n_unphased={len(rec.unphased)} "
              f"n_total={rec.n_total()}")
        for tag, lst in (("  hap1", rec.hap1), ("  hap2", rec.hap2),
                         ("  unphased", rec.unphased)):
            for v in lst:
                print(f"{tag}: {v.variant_id} (GT={v.genotypes[pid]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
