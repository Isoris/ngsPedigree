"""
test_offspring_phasing — reverse-segregation phasing across two filter
modes (segregating + hemizygous-only).
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.offspring_phasing import (  # noqa: E402
    OffspringPhasingRecord,
    VariantMarker,
    carrier_frequency,
    cluster_offspring_two_class,
    encode_dosage,
    filter_variants,
    hamming_like_distance,
    het_fraction,
    phase_offspring_interval,
    score_variant_class_separation,
    write_class_assignments_tsv,
    write_phasing_tsv,
)


# ----------------------------------------------------------------------
# encode_dosage
# ----------------------------------------------------------------------


class TestEncodeDosage(unittest.TestCase):
    def test_gt_strings(self):
        self.assertEqual(encode_dosage("0/0"), 0.0)
        self.assertEqual(encode_dosage("0/1"), 0.5)
        self.assertEqual(encode_dosage("1/0"), 0.5)
        self.assertEqual(encode_dosage("1/1"), 1.0)
        self.assertIsNone(encode_dosage("./."))
        self.assertIsNone(encode_dosage(""))
        self.assertIsNone(encode_dosage(None))

    def test_numeric(self):
        self.assertEqual(encode_dosage(0), 0.0)
        self.assertEqual(encode_dosage(1), 1.0)
        self.assertEqual(encode_dosage(2), 1.0)        # /2 normalisation
        self.assertEqual(encode_dosage(0.3), 0.3)
        self.assertIsNone(encode_dosage(-0.1))
        self.assertIsNone(encode_dosage(2.5))
        self.assertIsNone(encode_dosage("garbage"))


# ----------------------------------------------------------------------
# carrier_frequency / het_fraction
# ----------------------------------------------------------------------


class TestFrequencies(unittest.TestCase):
    def test_carrier_freq(self):
        self.assertAlmostEqual(carrier_frequency([0.0, 0.0, 0.5, 1.0]), 0.5)
        self.assertIsNone(carrier_frequency([]))

    def test_het_fraction(self):
        self.assertAlmostEqual(het_fraction([0.5, 0.5, 0.5, 0.0]), 0.75)
        self.assertEqual(het_fraction([0.0, 1.0, 0.0, 1.0]), 0.0)


# ----------------------------------------------------------------------
# filter_variants — both modes.
# ----------------------------------------------------------------------


class TestFilterVariants(unittest.TestCase):
    def _matrix(self):
        # V_seg: classic 0.5 segregator (half present, half absent)
        # V_fix: present in everyone (fixed, drop)
        # V_rare: present in 1 of 4 (below 0.20 threshold, drop)
        # V_hemi: every offspring heterozygous (situation 1)
        # V_high_miss: 80% missing
        samples = ["O1", "O2", "O3", "O4"]
        m = {
            "V_seg":   {"O1": "1/1", "O2": "1/1", "O3": "0/0", "O4": "0/0"},
            "V_fix":   {"O1": "1/1", "O2": "1/1", "O3": "1/1", "O4": "1/1"},
            "V_rare":  {"O1": "0/0", "O2": "0/0", "O3": "0/0", "O4": "0/0"},
            "V_hemi":  {"O1": "0/1", "O2": "0/1", "O3": "0/1", "O4": "0/1"},
            "V_high_miss": {"O1": "0/1"},
        }
        return samples, m

    def test_segregating_mode(self):
        samples, m = self._matrix()
        kept = filter_variants(
            ["V_seg", "V_fix", "V_rare", "V_hemi", "V_high_miss"],
            m, samples, mode="segregating",
        )
        self.assertIn("V_seg", kept)
        # V_hemi has dosage 0.5 for every offspring → carrier freq 1.0,
        # above max_freq 0.80 → excluded from segregating band.
        self.assertNotIn("V_hemi", kept)
        self.assertNotIn("V_fix", kept)
        self.assertNotIn("V_rare", kept)
        self.assertNotIn("V_high_miss", kept)

    def test_hemizygous_only_mode(self):
        samples, m = self._matrix()
        kept = filter_variants(
            ["V_seg", "V_fix", "V_rare", "V_hemi", "V_high_miss"],
            m, samples, mode="hemizygous_only",
        )
        self.assertIn("V_hemi", kept)
        self.assertNotIn("V_seg", kept)
        self.assertNotIn("V_fix", kept)
        self.assertNotIn("V_rare", kept)

    def test_unknown_mode_raises(self):
        with self.assertRaises(ValueError):
            filter_variants(["V"], {}, ["S"], mode="nope")


# Patch the earlier test that incorrectly expected V_hemi in segregating.
class TestSegregatingFilterStrict(unittest.TestCase):
    def test_hemi_excluded_from_segregating(self):
        # carrier freq for V_hemi is 1.0 (everyone has dosage > 0) →
        # excluded from segregating band [0.20, 0.80].
        samples = ["O1", "O2", "O3", "O4"]
        m = {"V_hemi": {s: "0/1" for s in samples}}
        kept = filter_variants(["V_hemi"], m, samples, mode="segregating")
        self.assertNotIn("V_hemi", kept)


# ----------------------------------------------------------------------
# Hamming-like distance.
# ----------------------------------------------------------------------


class TestDistance(unittest.TestCase):
    def test_basic(self):
        d, n = hamming_like_distance([0.0, 1.0, 0.5], [1.0, 1.0, 0.5])
        self.assertAlmostEqual(d, 1.0)
        self.assertEqual(n, 3)

    def test_missing_skipped(self):
        d, n = hamming_like_distance([0.0, None, 1.0], [1.0, 1.0, None])
        self.assertAlmostEqual(d, 1.0)
        self.assertEqual(n, 1)


# ----------------------------------------------------------------------
# Two-class clustering.
# ----------------------------------------------------------------------


class TestTwoClassClustering(unittest.TestCase):
    def test_clean_two_class_split(self):
        # 4 variants perfectly separate 4 offspring into 2 + 2.
        samples = ["O1", "O2", "O3", "O4"]
        variants = ["V1", "V2", "V3", "V4"]
        # O1, O2 carry the variants; O3, O4 do not.
        matrix = {
            v: {"O1": "1/1", "O2": "1/1", "O3": "0/0", "O4": "0/0"}
            for v in variants
        }
        classes = cluster_offspring_two_class(samples, variants, matrix)
        a_members = {s for s, c in classes.items() if c == "A"}
        b_members = {s for s, c in classes.items() if c == "B"}
        self.assertEqual(len(a_members), 2)
        self.assertEqual(len(b_members), 2)
        # O1/O2 should be together; O3/O4 should be together.
        self.assertTrue({"O1", "O2"} in (a_members, b_members))
        self.assertTrue({"O3", "O4"} in (a_members, b_members))

    def test_no_variants_unassigned(self):
        classes = cluster_offspring_two_class(
            ["O1", "O2"], [], {},
        )
        self.assertEqual(classes, {"O1": "U", "O2": "U"})


# ----------------------------------------------------------------------
# Variant scoring.
# ----------------------------------------------------------------------


class TestScoring(unittest.TestCase):
    def test_perfect_separator(self):
        samples = ["O1", "O2", "O3", "O4"]
        matrix = {"V": {"O1": "1/1", "O2": "1/1",
                         "O3": "0/0", "O4": "0/0"}}
        classes = {"O1": "A", "O2": "A", "O3": "B", "O4": "B"}
        n_a, n_b, ma, mb, score, n_obs, miss = (
            score_variant_class_separation("V", samples, matrix, classes)
        )
        self.assertEqual((n_a, n_b), (2, 2))
        self.assertAlmostEqual(score, 1.0)
        self.assertAlmostEqual(miss, 0.0)

    def test_useless_variant(self):
        samples = ["O1", "O2", "O3", "O4"]
        matrix = {"V": {s: "0/1" for s in samples}}
        classes = {"O1": "A", "O2": "A", "O3": "B", "O4": "B"}
        _, _, ma, mb, score, _, _ = score_variant_class_separation(
            "V", samples, matrix, classes,
        )
        self.assertAlmostEqual(score, 0.0)


# ----------------------------------------------------------------------
# Orchestrator.
# ----------------------------------------------------------------------


class TestPhaseOffspringInterval(unittest.TestCase):
    def _make_markers(self):
        return [
            VariantMarker("V1", "Chr1", 1_000, 1_001, "DEL"),
            VariantMarker("V2", "Chr1", 2_000, 2_001, "DEL"),
            VariantMarker("V3", "Chr1", 3_000, 3_001, "DUP"),
            VariantMarker("V_out", "Chr2", 999_000, 999_001, "DEL"),
        ]

    def test_segregating_returns_ranked_markers(self):
        samples = ["O1", "O2", "O3", "O4", "O5", "O6"]
        # All three Chr1 markers cleanly split O1-O3 from O4-O6.
        matrix = {
            "V1": {"O1": "1/1", "O2": "1/1", "O3": "1/1",
                   "O4": "0/0", "O5": "0/0", "O6": "0/0"},
            "V2": {"O1": "1/1", "O2": "1/1", "O3": "1/1",
                   "O4": "0/0", "O5": "0/0", "O6": "0/0"},
            "V3": {"O1": "1/1", "O2": "1/1", "O3": "1/1",
                   "O4": "0/0", "O5": "0/0", "O6": "0/0"},
            "V_out": {s: "0/0" for s in samples},
        }
        classes, records = phase_offspring_interval(
            family_id="FAM01", chrom="Chr1",
            interval_start=0, interval_end=10_000,
            markers=self._make_markers(),
            dosage_matrix=matrix,
            offspring=samples,
            mode="segregating",
        )
        # 3 in-interval markers retained.
        self.assertEqual(len(records), 3)
        # Top score should be perfect (1.0).
        self.assertAlmostEqual(records[0].marker_score, 1.0)
        # Sorted descending.
        scores = [r.marker_score for r in records]
        self.assertEqual(scores, sorted(scores, reverse=True))
        # Classes split into 3+3.
        class_a = [s for s, c in classes.items() if c == "A"]
        class_b = [s for s, c in classes.items() if c == "B"]
        self.assertEqual(len(class_a), 3)
        self.assertEqual(len(class_b), 3)

    def test_hemizygous_only_picks_all_het_variant(self):
        samples = ["O1", "O2", "O3", "O4"]
        markers = [
            VariantMarker("V_hemi_1", "Chr1", 100, 200, "DEL"),
            VariantMarker("V_hemi_2", "Chr1", 300, 400, "DEL"),
            VariantMarker("V_hemi_3", "Chr1", 500, 600, "DEL"),
            VariantMarker("V_seg", "Chr1", 700, 800, "DEL"),
        ]
        matrix = {
            "V_hemi_1": {s: "0/1" for s in samples},
            "V_hemi_2": {s: "0/1" for s in samples},
            "V_hemi_3": {s: "0/1" for s in samples},
            "V_seg":    {"O1": "1/1", "O2": "1/1",
                          "O3": "0/0", "O4": "0/0"},
        }
        _, records = phase_offspring_interval(
            family_id="FAM01", chrom="Chr1",
            interval_start=0, interval_end=10_000,
            markers=markers,
            dosage_matrix=matrix,
            offspring=samples,
            mode="hemizygous_only",
        )
        vids = {r.variant_id for r in records}
        self.assertEqual(vids, {"V_hemi_1", "V_hemi_2", "V_hemi_3"})
        # parental state hint should mention hemizygous depleted pattern.
        for r in records:
            self.assertEqual(r.inferred_parental_state,
                             "het_x_absent_or_hom_depleted")

    def test_interval_filter_drops_off_chrom(self):
        samples = ["O1", "O2", "O3", "O4"]
        markers = [
            VariantMarker("V_in", "Chr1", 1_500, 1_600, "DEL"),
            VariantMarker("V_out", "Chr2", 1_500, 1_600, "DEL"),
        ]
        matrix = {
            "V_in":  {"O1": "1/1", "O2": "0/0", "O3": "1/1", "O4": "0/0"},
            "V_out": {"O1": "1/1", "O2": "0/0", "O3": "1/1", "O4": "0/0"},
        }
        _, recs = phase_offspring_interval(
            family_id="FAM01", chrom="Chr1",
            interval_start=0, interval_end=10_000,
            markers=markers, dosage_matrix=matrix,
            offspring=samples, mode="segregating",
        )
        self.assertEqual({r.variant_id for r in recs}, {"V_in"})


# ----------------------------------------------------------------------
# TSV writers round-trip.
# ----------------------------------------------------------------------


class TestWriters(unittest.TestCase):
    def test_phasing_tsv_round_trip(self):
        rec = OffspringPhasingRecord(
            family_id="FAM01", chrom="Chr1",
            interval_start=0, interval_end=10_000,
            variant_id="V1", variant_type="DEL", pos=1_000,
            n_class_a=3, n_class_b=3,
            mean_dosage_class_a=1.0, mean_dosage_class_b=0.0,
            marker_score=1.0, n_informative=6, missingness=0.0,
            filter_mode="segregating",
            segregation_pattern="1:1",
            inferred_parental_state="AA_x_AB",
        )
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            write_phasing_tsv(p, [rec])
            lines = p.read_text().splitlines()
            self.assertIn("variant_id", lines[0])
            self.assertIn("FAM01", lines[1])
            self.assertIn("DEL", lines[1])
            self.assertIn("segregating", lines[1])
        finally:
            p.unlink()

    def test_class_assignments_tsv(self):
        classes = {"O1": "A", "O2": "B", "O3": "U"}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            p = Path(fh.name)
        try:
            write_class_assignments_tsv(
                p, "FAM01", "Chr1", 0, 10_000, classes, "segregating",
            )
            lines = p.read_text().splitlines()
            self.assertEqual(len(lines), 4)            # header + 3
            self.assertIn("sample_id", lines[0])
            self.assertIn("O1", lines[1])
        finally:
            p.unlink()


if __name__ == "__main__":
    unittest.main()
