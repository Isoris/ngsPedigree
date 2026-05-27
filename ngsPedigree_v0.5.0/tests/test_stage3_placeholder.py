"""
test_stage3_placeholder — Stage 3 placeholder loader assertions.

Run from the repo root or via ``run_tests.py``.
"""

from __future__ import annotations

import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.stage3_placeholder import (  # noqa: E402
    Stage3SchemaError,
    load_dyad_map,
    load_parent_phase,
    load_stage3,
    load_triad_map,
)

DYAD_FIX = THIS_DIR / "fixtures" / "synthetic_dyad"
TRIAD_FIX = THIS_DIR / "fixtures" / "synthetic_triad"


class TestDyadLoader(unittest.TestCase):
    def test_loads_two_segments(self):
        segs = load_dyad_map(DYAD_FIX / "inheritance_map_dyad.tsv")
        self.assertEqual(len(segs), 2)

    def test_segment_fields(self):
        segs = load_dyad_map(DYAD_FIX / "inheritance_map_dyad.tsv")
        s1, s2 = segs
        self.assertEqual(s1.dyad_id, "dyad_PA_OA")
        self.assertEqual(s1.parent_sample_id, "P_A")
        self.assertEqual(s1.offspring_sample_id, "O_A")
        self.assertEqual(s1.chrom, "LG01")
        self.assertEqual(s1.seg_start, 0)
        self.assertEqual(s1.seg_end, 1000)
        self.assertEqual(s1.parental_hap_inherited, "1")
        self.assertEqual(s1.segment_confidence, "Gold")
        self.assertFalse(s1.recomb_event_left)
        self.assertTrue(s1.recomb_event_right)
        self.assertEqual(s1.n_informative_markers, 42)
        self.assertEqual(s2.segment_confidence, "Silver")
        self.assertEqual(s2.parental_hap_inherited, "2")

    def test_bad_confidence_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write(textwrap.dedent("""\
                dyad_id\tparent_sample_id\toffspring_sample_id\tchrom\tseg_start\tseg_end\tparental_hap_inherited\tsegment_confidence\trecomb_event_left\trecomb_event_right\tn_informative_markers
                d\tp\to\tL\t0\t10\t1\tPlatinum\tfalse\tfalse\t3
            """))
            tmp = Path(fh.name)
        try:
            with self.assertRaises(Stage3SchemaError):
                load_dyad_map(tmp)
        finally:
            tmp.unlink()

    def test_bad_segment_range_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write(textwrap.dedent("""\
                dyad_id\tparent_sample_id\toffspring_sample_id\tchrom\tseg_start\tseg_end\tparental_hap_inherited\tsegment_confidence\trecomb_event_left\trecomb_event_right\tn_informative_markers
                d\tp\to\tL\t100\t50\t1\tGold\tfalse\tfalse\t3
            """))
            tmp = Path(fh.name)
        try:
            with self.assertRaises(Stage3SchemaError):
                load_dyad_map(tmp)
        finally:
            tmp.unlink()

    def test_missing_column_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("dyad_id\tparent_sample_id\n")
            fh.write("d\tp\n")
            tmp = Path(fh.name)
        try:
            with self.assertRaises(Stage3SchemaError):
                load_dyad_map(tmp)
        finally:
            tmp.unlink()


class TestTriadLoader(unittest.TestCase):
    def test_loads_two_segments(self):
        segs = load_triad_map(TRIAD_FIX / "inheritance_map_triad.tsv")
        self.assertEqual(len(segs), 2)

    def test_ambiguous_maternal_hap(self):
        segs = load_triad_map(TRIAD_FIX / "inheritance_map_triad.tsv")
        seg_bronze = [s for s in segs if s.segment_confidence == "Bronze"][0]
        self.assertEqual(seg_bronze.maternal_hap_inherited, "ambiguous")
        self.assertEqual(seg_bronze.paternal_hap_inherited, "2")


class TestParentPhaseLoader(unittest.TestCase):
    def test_dyad_phase_loads(self):
        phase = load_parent_phase(DYAD_FIX / "parent_phase.tsv")
        self.assertEqual(len(phase), 3)
        self.assertEqual(phase[("P_A", "LG01:100:A:T")], "1")
        self.assertEqual(phase[("P_A", "LG01:200:G:C")], "2")
        self.assertEqual(phase[("P_A", "LG01:400:T:G")], "1")
        self.assertNotIn(("P_A", "LG01:600:G:T"), phase)

    def test_triad_phase_loads(self):
        phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        self.assertEqual(phase[("P_pat", "LG01:100:A:T")], "1")
        self.assertEqual(phase[("P_mat", "LG01:100:A:T")], "2")
        self.assertEqual(phase[("P_pat", "LG01:500:C:G")], "2")
        self.assertEqual(phase[("P_mat", "LG01:800:T:A")], "1")

    def test_contradictory_entries_raise(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("parent_sample_id\tchrom\tpos\tref\talt\tparental_hap\n")
            fh.write("P\tL\t1\tA\tT\t1\n")
            fh.write("P\tL\t1\tA\tT\t2\n")
            tmp = Path(fh.name)
        try:
            with self.assertRaises(Stage3SchemaError):
                load_parent_phase(tmp)
        finally:
            tmp.unlink()


class TestComposite(unittest.TestCase):
    def test_load_stage3_aggregates(self):
        inputs = load_stage3(
            dyad_map_path=DYAD_FIX / "inheritance_map_dyad.tsv",
            parent_phase_path=DYAD_FIX / "parent_phase.tsv",
        )
        self.assertEqual(len(inputs.dyad_segments), 2)
        self.assertEqual(len(inputs.parent_phase), 3)
        self.assertEqual(inputs.dyad_ids(), ["dyad_PA_OA"])
        self.assertEqual(inputs.triad_ids(), [])

    def test_load_stage3_triad(self):
        inputs = load_stage3(
            triad_map_path=TRIAD_FIX / "inheritance_map_triad.tsv",
            parent_phase_path=TRIAD_FIX / "parent_phase.tsv",
        )
        self.assertEqual(len(inputs.triad_segments), 2)
        self.assertEqual(len(inputs.parent_phase), 4)
        self.assertEqual(inputs.triad_ids(), ["triad_T1"])


if __name__ == "__main__":
    unittest.main()
