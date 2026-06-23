"""
test_viz_data_emitters — every figure's underlying data table.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.viz_data_emitters import (  # noqa: E402
    _pearson,
    detect_kin_groups, emit_close_kin_groups,
    emit_edge_class_counts, emit_genome_event_timeline,
    emit_ideogram_segments, emit_kinship_matrix, emit_mating_risk_matrix,
    emit_pairwise_metrics, emit_pedigree_network,
    emit_per_chromosome_events,
    jaccard_del,
)


def _mock_pair(a, b, theta=0.25, IBS0=0.001, n=100):
    class P:
        pass
    p = P()
    p.sample_a = a
    p.sample_b = b
    p.theta = theta
    p.IBS0 = IBS0
    p.n_informative = n
    return p


# ----------------------------------------------------------------------
# Jaccard.
# ----------------------------------------------------------------------


class TestJaccard(unittest.TestCase):
    def test_identical_samples(self):
        # Both samples are DEL-carriers at the same 5 markers.
        g = {f"M{i}": {"A": "0/1", "B": "0/1"} for i in range(5)}
        self.assertAlmostEqual(jaccard_del("A", "B", g), 1.0)

    def test_disjoint_samples(self):
        g = {f"M{i}": {"A": "0/1", "B": "0/0"} for i in range(5)}
        self.assertAlmostEqual(jaccard_del("A", "B", g), 0.0)

    def test_partial(self):
        g = {
            "M0": {"A": "0/1", "B": "0/1"},
            "M1": {"A": "0/1", "B": "0/0"},
            "M2": {"A": "0/0", "B": "0/1"},
            "M3": {"A": "0/1", "B": "0/1"},
        }
        # union = 4, intersect = 2 → 0.5
        self.assertAlmostEqual(jaccard_del("A", "B", g), 0.5)


# ----------------------------------------------------------------------
# Pairwise metrics emitter.
# ----------------------------------------------------------------------


class TestPairwiseMetrics(unittest.TestCase):
    def test_emit(self):
        pairs = [_mock_pair("A", "B"), _mock_pair("A", "C", theta=0.0,
                                                    IBS0=0.1)]
        g = {f"M{i}": {"A": "0/1", "B": "0/1", "C": "0/0"} for i in range(5)}
        ecbp = {("A", "B"): "parent_offspring",
                ("A", "C"): "unrelated"}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            rows = emit_pairwise_metrics(pairs, g, ecbp, p)
            self.assertEqual(len(rows), 2)
            data = p.read_text().splitlines()
            self.assertIn("Jaccard", data[0])
            self.assertIn("parent_offspring", data[1])
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Edge-class counts.
# ----------------------------------------------------------------------


class TestEdgeClassCounts(unittest.TestCase):
    def test_dedupe_pairs(self):
        ecbp = {
            ("A", "B"): "parent_offspring",
            ("B", "A"): "parent_offspring",   # mirror
            ("A", "C"): "unrelated",
            ("B", "C"): "unrelated",
        }
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            rows = emit_edge_class_counts(ecbp, p)
            counts = {r["edge_class"]: r["n_pairs"] for r in rows}
            self.assertEqual(counts["parent_offspring"], 1)
            self.assertEqual(counts["unrelated"], 2)
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Pedigree network nodes + edges.
# ----------------------------------------------------------------------


class TestPedigreeNetwork(unittest.TestCase):
    def test_emit_nodes_edges(self):
        samples = ["A", "B", "C"]
        ecbp = {("A", "B"): "parent_offspring",
                ("A", "C"): "unrelated"}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fhn:
            pn = Path(fhn.name)
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fhe:
            pe = Path(fhe.name)
        try:
            nodes, edges = emit_pedigree_network(
                samples, ecbp,
                edge_classes_kept=("parent_offspring",),
                out_nodes=pn, out_edges=pe,
            )
            self.assertEqual(len(nodes), 3)
            self.assertEqual(len(edges), 1)
            self.assertEqual(edges[0]["edge_class"], "parent_offspring")
            self.assertEqual({edges[0]["source"], edges[0]["target"]}, {"A", "B"})
        finally:
            pn.unlink()
            pe.unlink()


# ----------------------------------------------------------------------
# Close-kin groups + family aggregates.
# ----------------------------------------------------------------------


class TestKinGroups(unittest.TestCase):
    def test_detect(self):
        samples = ["A", "B", "C", "D", "E"]
        ecbp = {
            ("A", "B"): "parent_offspring",
            ("B", "C"): "full_sibling",
            ("D", "E"): "parent_offspring",
        }
        fam = detect_kin_groups(ecbp, samples)
        # A, B, C should all be in the same family; D + E in another;
        # E if not connected → its own.
        self.assertEqual(fam["A"], fam["B"])
        self.assertEqual(fam["B"], fam["C"])
        self.assertEqual(fam["D"], fam["E"])
        self.assertNotEqual(fam["A"], fam["D"])

    def test_aggregates(self):
        samples = ["A", "B", "C"]
        fam = {"A": "F001", "B": "F001", "C": "F002"}
        theta = {("A", "B"): 0.25}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            rows = emit_close_kin_groups(samples, fam, theta, out_path=p)
            f1 = next(r for r in rows if r["family_id"] == "F001")
            self.assertEqual(f1["n_members"], 2)
            self.assertAlmostEqual(f1["mean_within_family_theta"], 0.25, places=4)
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Kinship-matrix emitter.
# ----------------------------------------------------------------------


class TestKinshipMatrix(unittest.TestCase):
    def test_self_kinship_05(self):
        samples = ["A", "B"]
        theta = {("A", "B"): 0.0625}
        fam = {"A": "F1", "B": "F1"}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            rows, ordered = emit_kinship_matrix(samples, theta, fam, p)
            # 2×2 matrix → 4 rows
            self.assertEqual(len(rows), 4)
            aa = next(r for r in rows if r["sample_a"] == "A" and r["sample_b"] == "A")
            self.assertEqual(aa["theta"], 0.5)
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Mating-risk matrix.
# ----------------------------------------------------------------------


class TestMatingRisk(unittest.TestCase):
    def test_categories(self):
        females = ["F1", "F2"]
        males = ["M1", "M2"]
        theta = {
            ("F1", "M1"): 0.0,       # low risk
            ("F1", "M2"): 0.10,      # high risk
            ("F2", "M1"): 0.03,      # mid caution
            ("F2", "M2"): 0.001,     # low
        }
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            rows = emit_mating_risk_matrix(females, males, theta, p)
            cats = {(r["female_sample_id"], r["male_sample_id"]): r["risk_category"]
                     for r in rows}
            self.assertEqual(cats[("F1", "M1")], "low_risk_recommended")
            self.assertEqual(cats[("F1", "M2")], "high_risk_avoid")
            self.assertEqual(cats[("F2", "M1")], "mid_caution")
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Pearson, power-iteration, PCA emitter.
# ----------------------------------------------------------------------


class TestPearson(unittest.TestCase):
    def test_perfect_positive(self):
        self.assertAlmostEqual(_pearson([1, 2, 3], [2, 4, 6]), 1.0)

    def test_perfect_negative(self):
        self.assertAlmostEqual(_pearson([1, 2, 3], [3, 2, 1]), -1.0)


# ----------------------------------------------------------------------
# Ideogram + event emitters smoke-test on minimal fixtures.
# ----------------------------------------------------------------------


class _MockSeg:
    def __init__(self, chrom, start, end, allele="REF",
                 confidence="Gold", recomb_left=False, recomb_right=False,
                 hap="1"):
        self.dyad_or_triad_id = "T1"
        self.offspring_sample_id = "C"
        self.parent_sample_id = "P"
        self.chrom = chrom
        self.seg_start = start
        self.seg_end = end
        self.transmitted_allele = allele
        self.n_informative_markers = 10
        self.confidence = confidence
        self.recomb_event_left = recomb_left
        self.recomb_event_right = recomb_right
        self.parental_hap_inherited = hap
        self.notes = ""


class TestIdeogramEmitter(unittest.TestCase):
    def test_emit_segments_and_event_summary(self):
        triad_maps = {
            "T1": {
                "paternal": [
                    _MockSeg("Chr1", 0, 50000),
                    _MockSeg("Chr1", 50001, 100000, allele="DEL",
                              recomb_left=True),
                ],
                "maternal": [
                    _MockSeg("Chr1", 0, 100000, allele="DEL"),
                ],
            }
        }
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fhs:
            ps = Path(fhs.name)
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fhe:
            pe = Path(fhe.name)
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fht:
            pt = Path(fht.name)
        try:
            rows = emit_ideogram_segments(triad_maps, ps)
            self.assertEqual(len(rows), 3)
            events = emit_per_chromosome_events(triad_maps, pe)
            pat = next(e for e in events if e["side"] == "paternal")
            self.assertEqual(pat["n_segments"], 2)
            self.assertGreaterEqual(pat["n_co"], 1)
            timeline = emit_genome_event_timeline(triad_maps, pt)
            self.assertGreaterEqual(len(timeline), 1)
            self.assertEqual(timeline[0]["event_type"], "crossover")
        finally:
            ps.unlink()
            pe.unlink()
            pt.unlink()
