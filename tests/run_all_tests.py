#!/usr/bin/env python3
"""Master test runner. Runs all three test suites and aggregates results."""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

SUITES = [
    ("Stage 1 — 13-sample fixture",
     HERE / "synthetic_12sample_fixture" / "run_tests.py"),
    ("Stage 1 — realistic 226-sample fixture",
     HERE / "synthetic_226_realistic" / "run_tests_226.py"),
    ("Stage 2 — per-chromosome QC (30-chrom fixture)",
     HERE / "synthetic_per_chrom" / "run_tests_stage2.py"),
    ("Stage 4 — HPP v0.5.0 (MVP 1 + adapters)",
     HERE.parent / "ngsPedigree_v0.5.0" / "tests" / "run_tests.py"),
]

print("=" * 70)
print("ngsPedigree — full test suite")
print("=" * 70)

all_passed = True
for name, script in SUITES:
    print(f"\n>>> {name}")
    res = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
    # Capture and print just the final summary lines
    lines = res.stdout.splitlines()
    for line in lines[-6:]:
        if "PASS" in line or "FAIL" in line or "===" in line:
            print(f"    {line}")
    if res.returncode != 0:
        all_passed = False
        print(f"    !!! suite FAILED (exit code {res.returncode})")
        print("    --- stdout tail ---")
        for line in lines[-15:]:
            print(f"    {line}")
        print("    --- stderr tail ---")
        for line in res.stderr.splitlines()[-15:]:
            print(f"    {line}")

print()
print("=" * 70)
if all_passed:
    print("ALL SUITES PASS")
else:
    print("AT LEAST ONE SUITE FAILED")
    sys.exit(1)
print("=" * 70)
