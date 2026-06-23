"""
test_hemizygous_check — Mendelian DEL transmission rules + triad
direction inference (the "fake trio" discriminator).
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.hemizygous_check import (  # noqa: E402
    ALLOWED_GT,
    DEL_MARKERS_SCHEMA,
    DelMarkersAdapterError,
    best_direction,
    dyad_del_strong_incompatible,
    is_del_mendelian_compatible,
    load_del_markers,
    normalise_gt,
    score_all_three_directions,
    score_dyad,
    score_triad_direction,
)

FIX = THIS_DIR / "fixtures" / "synthetic_hemizygous"


# ----------------------------------------------------------------------
# Punnett table — exhaustive Mendelian checks.
# ----------------------------------------------------------------------


class TestMendelianTable(unittest.TestCase):
    def test_hom_x_hom_opposite_forces_het(self):
        self.assertTrue(is_del_mendelian_compatible("0/0", "1/1", "0/1"))
        self.assertFalse(is_del_mendelian_compatible("0/0", "1/1", "0/0"))
        self.assertFalse(is_del_mendelian_compatible("0/0", "1/1", "1/1"))

    def test_hom_ref_x_het(self):
        self.assertTrue(is_del_mendelian_compatible("0/0", "0/1", "0/0"))
        self.assertTrue(is_del_mendelian_compatible("0/0", "0/1", "0/1"))
        self.assertFalse(is_del_mendelian_compatible("0/0", "0/1", "1/1"))

    def test_hom_inv_x_het(self):
        self.assertFalse(is_del_mendelian_compatible("1/1", "0/1", "0/0"))
        self.assertTrue(is_del_mendelian_compatible("1/1", "0/1", "0/1"))
        self.assertTrue(is_del_mendelian_compatible("1/1", "0/1", "1/1"))

    def test_het_x_het_allows_all(self):
        for c in ("0/0", "0/1", "1/1"):
            self.assertTrue(is_del_mendelian_compatible("0/1", "0/1", c))

    def test_hom_ref_x_hom_ref(self):
        self.assertTrue(is_del_mendelian_compatible("0/0", "0/0", "0/0"))
        self.assertFalse(is_del_mendelian_compatible("0/0", "0/0", "0/1"))
        self.assertFalse(is_del_mendelian_compatible("0/0", "0/0", "1/1"))

    def test_missing_returns_none(self):
        self.assertIsNone(is_del_mendelian_compatible("./.", "0/0", "0/1"))
        self.assertIsNone(is_del_mendelian_compatible("0/0", "./.", "0/1"))
        self.assertIsNone(is_del_mendelian_compatible("0/0", "0/0", "./."))

    def test_relaxed_hom_del_makes_marker_uninformative(self):
        # 1/1 parent under strict mode: incompatible if child=0/0.
        self.assertFalse(
            is_del_mendelian_compatible("1/1", "0/1", "0/0", strict_hom_del=True))
        # Relaxed: returns None (uninformative) rather than False.
        self.assertIsNone(
            is_del_mendelian_compatible("1/1", "0/1", "0/0", strict_hom_del=False))


# ----------------------------------------------------------------------
# Dyad helper.
# ----------------------------------------------------------------------


class TestDyadHelper(unittest.TestCase):
    def test_opposite_homozygotes_strong_incompatible(self):
        self.assertTrue(dyad_del_strong_incompatible("0/0", "1/1"))
        self.assertTrue(dyad_del_strong_incompatible("1/1", "0/0"))

    def test_other_cases_compatible(self):
        for p, o in [("0/0", "0/0"), ("0/0", "0/1"),
                     ("0/1", "0/0"), ("0/1", "0/1"), ("0/1", "1/1"),
                     ("1/1", "0/1"), ("1/1", "1/1")]:
            self.assertFalse(dyad_del_strong_incompatible(p, o))

    def test_relaxed_hom_del_yields_none_on_opposite_hom(self):
        # Strict mode flags the opposite-hom DEL dyad as incompatible.
        self.assertTrue(
            dyad_del_strong_incompatible("0/0", "1/1", strict_hom_del=True))
        # Relaxed treats it as untestable (short-read 1/1 caveat).
        self.assertIsNone(
            dyad_del_strong_incompatible("0/0", "1/1", strict_hom_del=False))


# ----------------------------------------------------------------------
# Genotype normalisation.
# ----------------------------------------------------------------------


class TestNormalise(unittest.TestCase):
    def test_phased_to_unphased(self):
        self.assertEqual(normalise_gt("0|1"), "0/1")
        self.assertEqual(normalise_gt("1|1"), "1/1")

    def test_missing_variants(self):
        for x in (None, "", "./.", "."):
            self.assertEqual(normalise_gt(x), "./.")


# ----------------------------------------------------------------------
# Triad direction scoring + best-direction picker.
# ----------------------------------------------------------------------


class TestTriadDirection(unittest.TestCase):
    def setUp(self):
        # Mendelian-clean (P_F + P_M -> C) across 10 markers.
        self.calls = load_del_markers(FIX / "del_markers.json")

    def test_correct_direction_has_zero_errors(self):
        s = score_triad_direction(
            p1="P_F", p2="P_M", child="C",
            del_calls_by_marker=self.calls,
        )
        self.assertEqual(s.n_incompatible, 0)
        self.assertEqual(s.n_informative, 10)
        self.assertEqual(s.error_rate, 0.0)

    def test_wrong_direction_picks_up_errors(self):
        # Try (P_F + C -> P_M). At several markers the math fails:
        #   DEL_002: (P_F=0/0, C=0/0) -> P_M expected 0/0, observed 0/1 → error
        s = score_triad_direction(
            p1="P_F", p2="C", child="P_M",
            del_calls_by_marker=self.calls,
        )
        self.assertGreater(s.n_incompatible, 0)

    def test_best_direction_picks_correct_trio(self):
        v = best_direction(
            a="P_F", b="P_M", c="C",
            del_calls_by_marker=self.calls,
            informative_marker_floor=5,
            min_margin=0.05,
        )
        # 10 markers all consistent with P_F + P_M -> C and inconsistent
        # with the other two permutations.
        self.assertEqual(v.best_direction, "P_F+P_M->C")
        self.assertEqual(v.best_error_rate, 0.0)
        self.assertIsNotNone(v.margin)
        self.assertGreater(v.margin, 0.0)

    def test_best_direction_refuses_with_no_margin(self):
        # Symmetric trio: every marker is HET in all three.
        sym = {
            f"DEL_{i:03d}": {"X": "0/1", "Y": "0/1", "Z": "0/1"}
            for i in range(10)
        }
        v = best_direction(
            a="X", b="Y", c="Z",
            del_calls_by_marker=sym,
            informative_marker_floor=5,
            min_margin=0.05,
        )
        # All three directions show zero errors → no margin → no winner.
        self.assertIsNone(v.best_direction)
        self.assertEqual(v.best_error_rate, 0.0)
        self.assertEqual(v.margin, 0.0)

    def test_best_direction_below_informative_floor(self):
        sparse = {"DEL_001": {"X": "0/0", "Y": "0/1", "Z": "0/1"}}
        v = best_direction(
            a="X", b="Y", c="Z",
            del_calls_by_marker=sparse,
            informative_marker_floor=5,
            min_margin=0.05,
        )
        self.assertIsNone(v.best_direction)


# ----------------------------------------------------------------------
# Dyad scoring (QC count).
# ----------------------------------------------------------------------


class TestDyadScoring(unittest.TestCase):
    def test_compatible_dyad_zero_strong_incompat(self):
        calls = load_del_markers(FIX / "del_markers.json")
        s = score_dyad(parent="P_F", offspring="C", del_calls_by_marker=calls)
        self.assertEqual(s.n_strong_incompatible, 0)

    def test_opposite_hom_dyad_flagged(self):
        bad = {
            "DEL_X": {"P": "0/0", "O": "1/1"},
            "DEL_Y": {"P": "0/0", "O": "0/0"},
        }
        s = score_dyad(parent="P", offspring="O", del_calls_by_marker=bad)
        self.assertEqual(s.n_strong_incompatible, 1)
        self.assertEqual(s.incompatible_markers, ["DEL_X"])


# ----------------------------------------------------------------------
# IN JSON adapter.
# ----------------------------------------------------------------------


class TestLoader(unittest.TestCase):
    def test_loads_fixture(self):
        calls = load_del_markers(FIX / "del_markers.json")
        self.assertEqual(len(calls), 10)
        self.assertEqual(calls["DEL_001"]["P_F"], "0/1")
        self.assertEqual(calls["DEL_005"]["P_F"], "1/1")

    def test_normalises_phased(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({
                "schema": DEL_MARKERS_SCHEMA,
                "rows": [
                    {"marker_id": "M", "sample_id": "X", "genotype": "0|1"},
                ],
            }, fh)
            p = Path(fh.name)
        try:
            self.assertEqual(load_del_markers(p)["M"]["X"], "0/1")
        finally:
            p.unlink()

    def test_bad_schema_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({"schema": "wrong", "rows": []}, fh)
            p = Path(fh.name)
        try:
            with self.assertRaises(DelMarkersAdapterError):
                load_del_markers(p)
        finally:
            p.unlink()

    def test_missing_required_field_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({
                "schema": DEL_MARKERS_SCHEMA,
                "rows": [{"marker_id": "M", "sample_id": "X"}],
            }, fh)
            p = Path(fh.name)
        try:
            with self.assertRaises(DelMarkersAdapterError):
                load_del_markers(p)
        finally:
            p.unlink()
