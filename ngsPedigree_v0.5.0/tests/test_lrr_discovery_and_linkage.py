"""
test_lrr_discovery_and_linkage — de novo LRR discovery from DEL
correlations (bloc 16) + per-DEL arrangement-linkage classifier (bloc 17).
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.del_arrangement_linkage import (  # noqa: E402
    classify_del_linkage,
)
from hpp.del_inheritance import DelMarkerLocus  # noqa: E402
from hpp.lrr_discovery import (  # noqa: E402
    discover_candidate_lrrs,
    marker_pair_correlation,
    pearson_corr,
    window_mean_correlation,
)


# ----------------------------------------------------------------------
# Pearson primitives.
# ----------------------------------------------------------------------


class TestPearson(unittest.TestCase):
    def test_perfect_positive(self):
        self.assertAlmostEqual(pearson_corr([0, 1, 2, 0, 1, 2],
                                              [0, 1, 2, 0, 1, 2]), 1.0)

    def test_perfect_negative(self):
        self.assertAlmostEqual(pearson_corr([0, 1, 2], [2, 1, 0]), -1.0)

    def test_zero_correlation(self):
        # constant in y → variance 0 → None
        self.assertIsNone(pearson_corr([0, 1, 2], [1, 1, 1]))

    def test_too_few_samples(self):
        self.assertIsNone(pearson_corr([0, 1], [0, 1]))


# ----------------------------------------------------------------------
# Window correlation.
# ----------------------------------------------------------------------


class TestWindowCorrelation(unittest.TestCase):
    def test_high_correlation_recognised(self):
        # 4 markers fully co-segregating across 5 samples
        gmatrix = {
            f"M{i}": {f"S{j}": gt for j, gt in enumerate(
                ["0/0", "0/1", "1/1", "0/0", "0/1"])}
            for i in range(4)
        }
        m, n_pairs = window_mean_correlation(
            ["M0", "M1", "M2", "M3"], gmatrix, ["S0", "S1", "S2", "S3", "S4"],
        )
        self.assertEqual(n_pairs, 6)
        self.assertAlmostEqual(m, 1.0)

    def test_uncorrelated_markers(self):
        # 4 random markers — average correlation should be low
        gmatrix = {
            "M0": {f"S{j}": gt for j, gt in enumerate(["0/0", "0/1", "1/1", "0/0", "0/1"])},
            "M1": {f"S{j}": gt for j, gt in enumerate(["0/1", "1/1", "0/0", "0/1", "0/0"])},
            "M2": {f"S{j}": gt for j, gt in enumerate(["1/1", "0/0", "0/1", "1/1", "0/0"])},
            "M3": {f"S{j}": gt for j, gt in enumerate(["0/0", "1/1", "0/1", "1/1", "0/1"])},
        }
        m, n_pairs = window_mean_correlation(
            ["M0", "M1", "M2", "M3"], gmatrix, ["S0", "S1", "S2", "S3", "S4"],
        )
        self.assertLess(m, 0.99)


# ----------------------------------------------------------------------
# De novo LRR discovery.
# ----------------------------------------------------------------------


class TestDiscovery(unittest.TestCase):
    def test_discovers_high_correlation_cluster(self):
        # Build 6 markers on Chr1 spanning a 200kb window, all
        # perfectly co-segregating across 6 samples.
        loci = [DelMarkerLocus(f"DEL_{i:02d}", "Chr1",
                                 i * 30_000, i * 30_000 + 500)
                for i in range(6)]
        samples = [f"S{j}" for j in range(6)]
        # Perfect co-segregation: same dosage at all markers per sample
        per_sample = ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"]
        gmatrix = {l.marker_id: {samples[j]: per_sample[j]
                                  for j in range(len(samples))}
                   for l in loci}
        cands = discover_candidate_lrrs(
            loci=loci, genotype_matrix=gmatrix, samples=samples,
            window_size=200_000, step=100_000,
            min_markers_per_window=4,
            correlation_threshold=0.8,
        )
        self.assertGreaterEqual(len(cands), 1)
        self.assertEqual(cands[0].chrom, "Chr1")
        self.assertGreater(cands[0].mean_pairwise_correlation, 0.9)

    def test_skips_when_no_correlation_block(self):
        loci = [DelMarkerLocus(f"DEL_{i:02d}", "Chr1",
                                 i * 30_000, i * 30_000 + 500)
                for i in range(6)]
        samples = [f"S{j}" for j in range(6)]
        # Random uncorrelated patterns — different per marker
        gmatrix = {}
        patterns = [
            ["0/0", "0/1", "1/1", "0/0", "0/1", "1/1"],
            ["0/1", "1/1", "0/0", "0/1", "0/0", "1/1"],
            ["1/1", "0/0", "0/1", "1/1", "0/0", "0/1"],
            ["0/0", "1/1", "0/1", "1/1", "0/1", "0/0"],
            ["0/1", "0/0", "1/1", "0/0", "1/1", "0/1"],
            ["1/1", "0/1", "0/0", "0/1", "1/1", "0/0"],
        ]
        for i, l in enumerate(loci):
            gmatrix[l.marker_id] = {samples[j]: patterns[i][j] for j in range(6)}
        cands = discover_candidate_lrrs(
            loci=loci, genotype_matrix=gmatrix, samples=samples,
            window_size=200_000, step=100_000,
            min_markers_per_window=4,
            correlation_threshold=0.8,
        )
        self.assertEqual(len(cands), 0)


# ----------------------------------------------------------------------
# Per-DEL × per-LRR arrangement linkage (situation 1).
# ----------------------------------------------------------------------


class TestLinkage(unittest.TestCase):
    def test_arrangement_1_marker_clean(self):
        # 5 hom_ref samples DEL-absent, 5 het samples DEL hemizygous,
        # 5 hom_inv samples DEL homozygous → arrangement-1-linked.
        sample_arr = {}
        del_gts = {}
        for i in range(5):
            s = f"HR{i}"
            sample_arr[s] = "HOM_REF"
            del_gts[s] = "0/0"
        for i in range(5):
            s = f"H{i}"
            sample_arr[s] = "HET"
            del_gts[s] = "0/1"
        for i in range(5):
            s = f"HI{i}"
            sample_arr[s] = "HOM_INV"
            del_gts[s] = "1/1"
        r = classify_del_linkage(
            del_id="DEL_X", lrr_id="LRR_1", chrom="Chr1",
            del_start=1000, del_end=2000,
            sample_arrangement=sample_arr,
            del_genotypes=del_gts,
        )
        self.assertEqual(r.interpretation, "arrangement_1_marker")
        self.assertEqual(r.del_freq_hom_ref, 0.0)
        self.assertEqual(r.del_freq_hom_inv, 1.0)

    def test_situation_1_hom_inv_depleted(self):
        # Arrangement-1-linked DEL but no HOM_INV samples observed.
        sample_arr = {}
        del_gts = {}
        for i in range(5):
            s = f"HR{i}"
            sample_arr[s] = "HOM_REF"
            del_gts[s] = "0/0"
        for i in range(8):
            s = f"H{i}"
            sample_arr[s] = "HET"
            del_gts[s] = "0/1"
        # no HOM_INV observed at all
        r = classify_del_linkage(
            del_id="DEL_X", lrr_id="LRR_1", chrom="Chr1",
            del_start=1000, del_end=2000,
            sample_arrangement=sample_arr,
            del_genotypes=del_gts,
        )
        self.assertEqual(r.interpretation, "arrangement_1_marker_hom_depleted")
        self.assertIn("depletion", r.notes)
        # Honesty: the disclaimer must explicitly say "not lethality".
        self.assertIn("not lethality", r.notes.lower())

    def test_unlinked_del(self):
        # DEL frequency similar (~0.3) across all three classes.
        sample_arr, del_gts = {}, {}
        for i, (cls, gt) in enumerate([
            ("HOM_REF", "0/0"), ("HOM_REF", "0/1"), ("HOM_REF", "0/1"),
            ("HOM_REF", "0/0"), ("HOM_REF", "0/1"),
            ("HET", "0/1"), ("HET", "0/0"), ("HET", "0/1"),
            ("HET", "0/0"), ("HET", "0/1"),
            ("HOM_INV", "0/1"), ("HOM_INV", "0/0"), ("HOM_INV", "0/1"),
            ("HOM_INV", "0/0"), ("HOM_INV", "0/1"),
        ]):
            s = f"S{i}"
            sample_arr[s] = cls
            del_gts[s] = gt
        r = classify_del_linkage(
            del_id="DEL_X", lrr_id="LRR_1", chrom="Chr1",
            del_start=1000, del_end=2000,
            sample_arrangement=sample_arr,
            del_genotypes=del_gts,
        )
        self.assertEqual(r.interpretation, "unlinked")

    def test_ambiguous_with_one_class_only(self):
        # All samples HET → only one class observed → ambiguous.
        sample_arr = {f"S{i}": "HET" for i in range(5)}
        del_gts = {f"S{i}": "0/1" for i in range(5)}
        r = classify_del_linkage(
            del_id="DEL_X", lrr_id="LRR_1", chrom="Chr1",
            del_start=1000, del_end=2000,
            sample_arrangement=sample_arr,
            del_genotypes=del_gts,
        )
        self.assertEqual(r.interpretation, "ambiguous")


# ----------------------------------------------------------------------
# Master CLI smoke test.
# ----------------------------------------------------------------------


FIX_VCF = THIS_DIR / "fixtures" / "synthetic_svvcfs"


class TestFullPipelineCLI(unittest.TestCase):
    def test_runs_with_discovery(self):
        script = THIS_DIR.parent / "scripts" / "10_run_full_pipeline.py"
        with tempfile.TemporaryDirectory() as td:
            outdir = Path(td)
            res = subprocess.run(
                [sys.executable, str(script),
                 "--delly", str(FIX_VCF / "delly.vcf"),
                 "--manta", str(FIX_VCF / "manta.vcf"),
                 "--outdir", str(outdir),
                 "--discover-lrrs",
                 "--window-size", "50000",
                 "--min-markers-per-window", "3",
                 "--correlation-threshold", "0.4",
                 "--min-informative-pair", "5"],
                capture_output=True, text=True, check=True,
            )
            for f in ("merged_del_catalogue.json",
                      "pairwise_relationship_classification.tsv",
                      "candidate_PO_pairs.tsv",
                      "detected_triads.tsv",
                      "inheritance_segments.tsv",
                      "discovered_lrrs.tsv"):
                self.assertTrue((outdir / f).exists(), f"missing {f}")

    def test_runs_without_lrr_inputs(self):
        # No --list-of-lrr, no --discover-lrrs. Should still produce
        # the pedigree+inheritance side, skipping LRR analyses.
        script = THIS_DIR.parent / "scripts" / "10_run_full_pipeline.py"
        with tempfile.TemporaryDirectory() as td:
            outdir = Path(td)
            subprocess.run(
                [sys.executable, str(script),
                 "--delly", str(FIX_VCF / "delly.vcf"),
                 "--manta", str(FIX_VCF / "manta.vcf"),
                 "--outdir", str(outdir),
                 "--min-informative-pair", "5"],
                capture_output=True, text=True, check=True,
            )
            self.assertTrue((outdir / "candidate_PO_pairs.tsv").exists())
            self.assertFalse((outdir / "discovered_lrrs.tsv").exists())
            self.assertFalse((outdir / "lrr_enrichment.tsv").exists())
