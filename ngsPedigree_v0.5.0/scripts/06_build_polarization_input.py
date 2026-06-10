#!/usr/bin/env python3
"""
06_build_polarization_input.py — assemble the polarization IN JSON
from ngsPedigree Stage 1/2 outputs + PCAngsd karyotype calls.

This closes the gap between the pedigree-classification half of
ngsPedigree (Stages 1-2, working from ngsRelate coefficients) and the
inversion-polarization half (which feeds ngsTracts via
scripts/05_polarize_inversion.py).

Usage:
  python 06_build_polarization_input.py \\
      --stage1-edges  PATH \\
      --stage1-roster PATH \\
      [--stage2-edges PATH] \\
      --karyotype     PATH \\
      --inversion-id  STR \\
      --polarity-hint {band_0_is_REF,band_0_is_INV} \\
      --out           PATH

Karyotype TSV format:
  Required columns: sample_id, band  (band ∈ {0, 1, 2}).
  Optional columns: confidence, inversion_id  (when inversion_id is
  present rows are filtered to --inversion-id, otherwise every row
  is taken to belong to --inversion-id).

Downstream:
  python 05_polarize_inversion.py --in <OUT> [--mtdna ...] --out ...
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.pedigree_extract import build_polarization_input  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage1-edges", required=True,
                    help="Stage 1 pairwise_relationship_classification.tsv")
    ap.add_argument("--stage1-roster", required=True,
                    help="Stage 1 family_hub_roster.tsv")
    ap.add_argument("--stage2-edges", default=None,
                    help=("optional Stage 2 extended pairwise table; if "
                          "supplied, pairs whose pair_review_flag != 'OK' "
                          "are excluded from dyad/triad construction"))
    ap.add_argument("--karyotype", required=True,
                    help="karyotype calls TSV (sample_id, band, [confidence, inversion_id])")
    ap.add_argument("--inversion-id", required=True)
    ap.add_argument("--polarity-hint", required=True,
                    choices=["band_0_is_REF", "band_0_is_INV"])
    ap.add_argument("--out", required=True,
                    help="polarization IN JSON (ngspedigree_karyotype_calls_in_v1)")
    args = ap.parse_args()

    bundle = build_polarization_input(
        stage1_edges_path=args.stage1_edges,
        stage1_roster_path=args.stage1_roster,
        karyotype_path=args.karyotype,
        inversion_id=args.inversion_id,
        polarity_hint=args.polarity_hint,
        stage2_edges_path=args.stage2_edges,
    )

    for w in bundle.warnings:
        print(f"[pedigree-extract] warn: {w}", file=sys.stderr)

    print(f"[pedigree-extract] inversion={bundle.inversion_id}  "
          f"karyotype_calls={len(bundle.karyotype_calls)}  "
          f"dyads={len(bundle.dyads)}  triads={len(bundle.triads)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(bundle.to_in_json(), fh, indent=2)
    print(f"[pedigree-extract] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
