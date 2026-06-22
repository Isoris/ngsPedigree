"""
test_lrr_enrichment — family-based OR enrichment for candidate LRRs.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.catalogue_merge import merge_two_callers, to_genotype_matrix  # noqa: E402
from hpp.del_inheritance import DelMarkerLocus  # noqa: E402
from hpp.lrr_enrichment import (  # noqa: E402
    LRRInterval,
    Relationship,
    classify_region_dyad,
    classify_region_triad,
    compute_enrichment,
    compute_odds_ratio,
    load_lrr_list,
    sample_background_regions,
)
from hpp.vcf_sv import read_del_calls  # noqa: E402
import random


FIX_VCF = THIS_DIR / "fixtures" / "synthetic_svvcfs"
FIX_LRR = THIS_DIR / "fixtures" / "synthetic_lrrs"


# ----------------------------------------------------------------------
# LRR TSV loader.
# ----------------------------------------------------------------------


class TestLoader(unittest.TestCase):
    def test_loads_two_lrrs(self):
        lrrs = load_lrr_list(FIX_LRR / "lrr_list.tsv")
        self.assertEqual(len(lrrs), 2)
        self.assertEqual(lrrs[0].lrr_id, "LRR_001")
        self.assertEqual(lrrs[0].chrom, "Chr1")
        self.assertEqual(lrrs[0].length, 90000)


# ----------------------------------------------------------------------
# Odds-ratio + CI + Haldane correction.
# ----------------------------------------------------------------------


class TestOddsRatio(unittest.TestCase):
    def test_perfect_enrichment(self):
        # 10 inside compatible, 0 incompatible vs 0 outside compatible, 10 incompat
        # → Haldane-corrected huge OR
        o = compute_odds_ratio(10, 0, 0, 10)
        self.assertTrue(o.haldane_corrected)
        self.assertGreater(o.odds_ratio, 1.0)

    def test_no_enrichment(self):
        o = compute_odds_ratio(5, 5, 5, 5)
        self.assertAlmostEqual(o.odds_ratio, 1.0, places=4)
        self.assertFalse(o.haldane_corrected)
        # CI for OR=1 spans 1
        self.assertLess(o.ci_low, 1.0)
        self.assertGreater(o.ci_high, 1.0)


# ----------------------------------------------------------------------
# Background sampler.
# ----------------------------------------------------------------------


class TestBackgroundSampler(unittest.TestCase):
    def test_avoids_excluded_intervals(self):
        rng = random.Random(0)
        bg = sample_background_regions(
            chrom="Chr1", chrom_length=100000, region_length=10000,
            n_regions=5, exclude=[(20000, 80000)],
            rng=rng,
        )
        # all background windows must fall outside (20000, 80000)
        for lo, hi in bg:
            self.assertTrue(hi <= 20000 or lo >= 80000)

    def test_zero_excl_returns_n_regions(self):
        rng = random.Random(0)
        bg = sample_background_regions(
            chrom="Chr1", chrom_length=1000000, region_length=10000,
            n_regions=10, exclude=[],
            rng=rng,
        )
        self.assertEqual(len(bg), 10)


# ----------------------------------------------------------------------
# Region classification (dyad + triad).
# ----------------------------------------------------------------------


class TestRegionClassification(unittest.TestCase):
    def setUp(self):
        merged = merge_two_callers(
            read_del_calls(FIX_VCF / "delly.vcf"),
            read_del_calls(FIX_VCF / "manta.vcf"),
        )
        self.loci = [
            DelMarkerLocus(m.marker_id, m.chrom, m.start, m.end)
            for m in merged
        ]
        self.gmatrix = to_genotype_matrix(merged)

    def test_triad_classification_on_lrr(self):
        # LRR Chr1:0-90000 — the triad P_F + P_M → C is Mendelian-
        # consistent here. With dominance threshold low enough we should
        # see block_compatible classifications on at least one side.
        n_inf, n_dom, dom, err, ok = classify_region_triad(
            paternal_sample_id="P_F",
            maternal_sample_id="P_M",
            offspring_sample_id="C",
            chrom="Chr1", start=0, end=90000,
            loci=self.loci, genotype_matrix=self.gmatrix,
            min_markers=3, dominance_threshold=0.5,
        )
        self.assertGreater(n_inf, 0)
        self.assertEqual(err, 0)


# ----------------------------------------------------------------------
# End-to-end: compute_enrichment over the fixture's LRRs + triad.
# ----------------------------------------------------------------------


class TestComputeEnrichment(unittest.TestCase):
    def setUp(self):
        merged = merge_two_callers(
            read_del_calls(FIX_VCF / "delly.vcf"),
            read_del_calls(FIX_VCF / "manta.vcf"),
        )
        self.loci = [
            DelMarkerLocus(m.marker_id, m.chrom, m.start, m.end)
            for m in merged
        ]
        self.gmatrix = to_genotype_matrix(merged)
        self.lrrs = load_lrr_list(FIX_LRR / "lrr_list.tsv")
        self.rels = [
            Relationship(
                relationship_id="triad_T1",
                relationship_type="triad",
                paternal_sample_id="P_F",
                maternal_sample_id="P_M",
                parent_sample_id=None,
                offspring_sample_id="C",
            ),
        ]

    def test_emits_one_row_per_lrr(self):
        enrichments, classifications = compute_enrichment(
            lrrs=self.lrrs,
            relationships=self.rels,
            loci=self.loci,
            genotype_matrix=self.gmatrix,
            n_background_per_lrr=5,
            seed=1,
            min_markers=2,
            dominance_threshold=0.5,
        )
        self.assertEqual(len(enrichments), len(self.lrrs))
        for e in enrichments:
            self.assertEqual(e.n_triads, 1)
            self.assertEqual(e.n_dyads, 0)
            # combined OR may be None when the background classification
            # produces zero informative cells; just check the field exists.
            self.assertIn("a", e.combined.to_dict())
