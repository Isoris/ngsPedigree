#!/usr/bin/env python3
"""
HPP MVP 1 master test runner. Runs all v0.5.0 unit-test modules and
prints a summary line at the end. Stdlib only — no external deps.

Usage:
  python ngsPedigree_v0.5.0/tests/run_tests.py
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))
sys.path.insert(0, str(THIS_DIR))


def main() -> int:
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir=str(THIS_DIR), pattern="test_*.py")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    print()
    print(f"tests run: {result.testsRun}")
    print(f"failures : {len(result.failures)}")
    print(f"errors   : {len(result.errors)}")
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
