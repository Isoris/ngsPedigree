#!/usr/bin/env python3
"""
07_run_synthetic_panels.py — run every synthetic test panel through the
full ngsPedigree → ngsTracts pipeline (Stage 1 shadow classifier →
roster construction → pedigree extract → mtDNA check → polarization →
transmission calling) and print a recovery report card.

The report shows, per panel, the edge-classification accuracy, the
per-class recall, the dyad/triad recovery rates against ground truth,
the mtDNA swap-detection rates (where injected), and the polarization
output. A PhD is the floor we are insuring; this lets you see at a
glance whether any panel is below an acceptable recovery rate.

Usage:
  python 07_run_synthetic_panels.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.recovery_harness import format_report, run_panel  # noqa: E402
from hpp.synthetic_panels import PANEL_BUILDERS  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", default=None,
                    help="run only this panel (default: all)")
    args = ap.parse_args()

    panels = (
        [(args.panel, *PANEL_BUILDERS[args.panel])]
        if args.panel and args.panel in PANEL_BUILDERS
        else [(k, *v) for k, v in PANEL_BUILDERS.items()]
    )
    if args.panel and args.panel not in PANEL_BUILDERS:
        sys.stderr.write(
            f"unknown panel {args.panel!r}; available: {list(PANEL_BUILDERS)}\n"
        )
        return 2

    print("=" * 72)
    print("ngsPedigree — synthetic-panel recovery report")
    print("=" * 72)
    print()

    for key, builder, cfg in panels:
        pedigree = builder(cfg)
        rep = run_panel(pedigree, cfg)
        print(format_report(rep))
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
