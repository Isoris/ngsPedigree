"""
test_chromosome_inheritance — independent-assortment-level
inheritance map (bloc 19). Distinct from the LRR-level inheritance
map in del_inheritance.py.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.chromosome_inheritance import (  # noqa: E402
    build_chromosome_inheritance_map,
    score_chromosome_inheritance,
)
from hpp.del_inheritance import DelMarkerLocus  # noqa: E402


def _make_loci(chrom, n, step=1000):
    return [
        DelMarkerLocus(
            marker_id=f"DEL_{chrom}_{i:03d}",
            chrom=chrom, start=i * step, end=i * step + 500,
        )
        for i in range(n)
    ]


# ----------------------------------------------------------------------
# Mendelian-exclusion logic (the "rejected" path).
# ----------------------------------------------------------------------


class TestExclusion(unittest.TestCase):
    def test_opposite_homozygote_rejects(self):
        loci = _make_loci("Chr1", 30)
        gmatrix = {}
        for i, loc in enumerate(loci):
            # parent 0/0, offspring 1/1 at every marker → fully impossible
            gmatrix[loc.marker_id] = {"P": "0/0", "C": "1/1"}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        self.assertEqual(s.n_excluding, 30)
        self.assertEqual(s.inheritance_support, "rejected")

    def test_clean_inheritance_no_exclusion(self):
        # Parent 0/0 at every marker → offspring 0/0 or 0/1 always
        # compatible. No contradictions.
        loci = _make_loci("Chr1", 30)
        gmatrix = {}
        for i, loc in enumerate(loci):
            og = "0/1" if i % 2 == 0 else "0/0"
            gmatrix[loc.marker_id] = {"P": "0/0", "C": og}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        self.assertEqual(s.n_excluding, 0)
        self.assertEqual(s.compatibility_rate, 1.0)
        self.assertNotEqual(s.inheritance_support, "rejected")


# ----------------------------------------------------------------------
# Parent-HET-DEL transmission rate ≈ 0.5 for true parent.
# ----------------------------------------------------------------------


class TestHetTransmission(unittest.TestCase):
    def test_half_inheritance_rate_recovered(self):
        # Parent het at every marker, offspring inherits ~50% as
        # 0/1 and ~50% as 0/0 — independent of which marker.
        loci = _make_loci("Chr1", 40)
        gmatrix = {}
        for i, loc in enumerate(loci):
            og = "0/1" if i % 2 == 0 else "0/0"
            gmatrix[loc.marker_id] = {"P": "0/1", "C": og}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        self.assertEqual(s.n_parent_het_markers, 40)
        self.assertEqual(s.n_parent_het_inherited, 20)
        self.assertAlmostEqual(s.het_inheritance_rate, 0.5, places=4)
        # Mendelian-clean → not rejected
        self.assertNotEqual(s.inheritance_support, "rejected")


# ----------------------------------------------------------------------
# Pearson r positive for true parent, near zero for stranger.
# ----------------------------------------------------------------------


class TestPearsonScore(unittest.TestCase):
    def test_true_parent_positive_pearson(self):
        # Parent dosage co-segregates with offspring dosage strongly.
        loci = _make_loci("Chr1", 30)
        gmatrix = {}
        # parent dosage [0, 0, 1, 1, 1, 2, 2, ...]; offspring same/+1
        for i, loc in enumerate(loci):
            pg = ["0/0", "0/1", "1/1"][i % 3]
            og = pg if i % 5 != 0 else "0/1"   # mostly track parent
            gmatrix[loc.marker_id] = {"P": pg, "C": og}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        # Some opposite-hom contradictions exist (P=0/0,C=0/1 OK but
        # the cycle in line above doesn't produce 0/0↔1/1). Confirm
        # Pearson r is strongly positive.
        self.assertIsNotNone(s.pearson_r)
        self.assertGreater(s.pearson_r, 0.5)

    def test_stranger_low_pearson(self):
        loci = _make_loci("Chr1", 30)
        # Parent and offspring uncorrelated.
        gmatrix = {}
        for i, loc in enumerate(loci):
            pg = ["0/0", "0/1", "1/1"][i % 3]
            og = ["0/1", "1/1", "0/0"][i % 3]   # cyclic shift
            gmatrix[loc.marker_id] = {"P": pg, "C": og}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        # If opposite-hom contradictions exist → rejected. Otherwise
        # support is at most 'compatible'.
        self.assertNotEqual(s.inheritance_support, "strong")


# ----------------------------------------------------------------------
# Map across chromosomes + best-parent picker.
# ----------------------------------------------------------------------


class TestMap(unittest.TestCase):
    def test_best_parent_per_chrom(self):
        loci = _make_loci("Chr1", 30) + _make_loci("Chr2", 30)
        gmatrix = {}
        # Parent A: matches offspring exactly on Chr1
        # Parent B: matches offspring exactly on Chr2
        for i, loc in enumerate(loci):
            if loc.chrom == "Chr1":
                pa = ["0/0", "0/1", "1/1"][i % 3]
                pb = "0/0"
                cg = pa if i % 5 != 0 else "0/1"
            else:
                pa = "0/0"
                pb = ["0/0", "0/1", "1/1"][i % 3]
                cg = pb if i % 5 != 0 else "0/1"
            gmatrix[loc.marker_id] = {"A": pa, "B": pb, "C": cg}
        m = build_chromosome_inheritance_map(
            offspring_sample_id="C",
            candidate_parents=["A", "B"],
            chromosomes=["Chr1", "Chr2"],
            loci=loci, genotype_matrix=gmatrix,
            min_markers=10,
        )
        # Should produce 2 chromosomes × 2 candidate parents = 4 scores
        self.assertEqual(len(m.scores), 4)
        # Best parent per chromosome should match the construction.
        self.assertEqual(m.best_per_chrom.get("Chr1"), "A")
        self.assertEqual(m.best_per_chrom.get("Chr2"), "B")

    def test_min_markers_threshold_marks_ambiguous(self):
        loci = _make_loci("Chr1", 5)
        gmatrix = {loc.marker_id: {"P": "0/0", "C": "0/0"} for loc in loci}
        s = score_chromosome_inheritance(
            offspring_sample_id="C", candidate_parent_sample_id="P",
            chrom="Chr1", loci=loci, genotype_matrix=gmatrix,
            min_markers=20,
        )
        self.assertEqual(s.inheritance_support, "ambiguous")
