#!/usr/bin/env python3
"""
02_ingest_stage3.py — load the placeholder Stage 3 inputs and emit a
summary to stdout. Useful as a smoke test against a fixture; not a
production pipeline step.

Usage:
  python 02_ingest_stage3.py \
      --dyad-map  PATH \
      [--triad-map PATH] \
      [--parent-phase PATH]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.stage3_placeholder import load_stage3  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dyad-map", required=False)
    ap.add_argument("--triad-map", required=False)
    ap.add_argument("--parent-phase", required=False)
    args = ap.parse_args()

    if not any([args.dyad_map, args.triad_map, args.parent_phase]):
        ap.error("at least one of --dyad-map / --triad-map / --parent-phase required")

    inputs = load_stage3(
        dyad_map_path=args.dyad_map,
        triad_map_path=args.triad_map,
        parent_phase_path=args.parent_phase,
    )

    print(f"dyad_segments     : {len(inputs.dyad_segments)}")
    print(f"triad_segments    : {len(inputs.triad_segments)}")
    print(f"parent_phase rows : {len(inputs.parent_phase)}")
    print(f"dyad_ids          : {inputs.dyad_ids()}")
    print(f"triad_ids         : {inputs.triad_ids()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
