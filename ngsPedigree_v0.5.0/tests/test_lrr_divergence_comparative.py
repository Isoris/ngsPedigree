"""
test_lrr_divergence_comparative — cross-LRR test of the pseudo-
overdominance hypothesis (situation 1 generalised genome-wide).
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.lrr_divergence_comparative import (  # noqa: E402
    LRRComparativeMetrics,
    compute_dxy_between_arrangements,
    compute_fis,
    compute_het_excess,
    compute_hom_het_ratios,
    compute_lrr_metrics,
    fit_comparative_regression,
    linear_regression,
)


# ----------------------------------------------------------------------
# FIS, het excess, hom/het ratios.
# ----------------------------------------------------------------------


class TestFis(unittest.TestCase):
    def test_hwe_population_fis_zero(self):
        # HWE expected: n_hom0=25, n_het=50, n_hom1=25 → p=0.5, H_obs=0.5,
        # H_exp=0.5 → FIS=0.
        self.assertAlmostEqual(compute_fis(25, 50, 25), 0.0, places=6)

    def test_pure_heterozygotes_fis_minus_one(self):
        # n_het=100, no homs → H_obs=1.0; p=0.5, H_exp=0.5 → FIS=1-2=-1
        self.assertAlmostEqual(compute_fis(0, 100, 0), -1.0, places=6)

    def test_homozygote_excess_positive_fis(self):
        # n_hom0=50, n_hom1=50, no het → H_obs=0; p=0.5, H_exp=0.5 → FIS=1
        self.assertAlmostEqual(compute_fis(50, 0, 50), 1.0, places=6)

    def test_zero_total_returns_none(self):
        self.assertIsNone(compute_fis(0, 0, 0))


class TestHetExcess(unittest.TestCase):
    def test_het_excess_above_one(self):
        self.assertAlmostEqual(compute_het_excess(0, 100, 0), 2.0, places=6)

    def test_het_consistent(self):
        self.assertAlmostEqual(compute_het_excess(25, 50, 25), 1.0, places=6)

    def test_het_deficit_below_one(self):
        # All homozygous → H_obs=0 → het excess = 0
        self.assertAlmostEqual(compute_het_excess(50, 0, 50), 0.0, places=6)


class TestHomHetRatios(unittest.TestCase):
    def test_haldane_corrected_zero_does_not_explode(self):
        r0, r1, r_min, miss = compute_hom_het_ratios(0, 20, 0)
        # both missing → flagged
        self.assertEqual(miss, "BOTH")
        # ratios are 0.5/20.5 ≈ 0.024, not infinity
        self.assertLess(r0, 0.1)
        self.assertLess(r1, 0.1)

    def test_situation_1_pattern(self):
        # HOM_INV missing
        r0, r1, r_min, miss = compute_hom_het_ratios(10, 30, 0)
        self.assertEqual(miss, "HOM_INV")
        self.assertGreater(r0, r1)
        self.assertEqual(r_min, r1)


# ----------------------------------------------------------------------
# dXY between arrangements.
# ----------------------------------------------------------------------


class TestDxy(unittest.TestCase):
    def test_identical_groups_dxy_zero(self):
        # Both arrangement groups have the same allele frequency.
        gmatrix = {
            "M1": {"A1": "0/0", "A2": "0/0", "B1": "0/0", "B2": "0/0"},
            "M2": {"A1": "1/1", "A2": "1/1", "B1": "1/1", "B2": "1/1"},
        }
        dxy, n = compute_dxy_between_arrangements(
            arr0_samples=["A1", "A2"], arr1_samples=["B1", "B2"],
            marker_ids=["M1", "M2"], genotype_matrix=gmatrix,
        )
        self.assertEqual(n, 2)
        self.assertAlmostEqual(dxy, 0.0, places=6)

    def test_perfect_divergence_dxy_one(self):
        # arr0 all 0/0, arr1 all 1/1 → maximum divergence
        gmatrix = {
            "M1": {"A1": "0/0", "A2": "0/0", "B1": "1/1", "B2": "1/1"},
            "M2": {"A1": "0/0", "A2": "0/0", "B1": "1/1", "B2": "1/1"},
        }
        dxy, n = compute_dxy_between_arrangements(
            arr0_samples=["A1", "A2"], arr1_samples=["B1", "B2"],
            marker_ids=["M1", "M2"], genotype_matrix=gmatrix,
        )
        self.assertAlmostEqual(dxy, 1.0, places=6)

    def test_empty_group_returns_none(self):
        gmatrix = {"M1": {"A1": "0/0"}}
        dxy, n = compute_dxy_between_arrangements(
            arr0_samples=["A1"], arr1_samples=[],
            marker_ids=["M1"], genotype_matrix=gmatrix,
        )
        self.assertIsNone(dxy)


# ----------------------------------------------------------------------
# Per-LRR metric driver.
# ----------------------------------------------------------------------


class _MiniLrr:
    def __init__(self, lrr_id, chrom, start, end):
        self.lrr_id = lrr_id
        self.chrom = chrom
        self.start = start
        self.end = end


class TestLrrMetrics(unittest.TestCase):
    def test_situation_1_missing_hom_inv(self):
        lrr = _MiniLrr("LRR_X", "Chr1", 0, 100000)
        arrangement = {f"S{i}": "HOM_REF" for i in range(5)}
        arrangement.update({f"H{i}": "HET" for i in range(20)})
        # no HOM_INV samples
        gmatrix = {f"M{j}": {sid: ("0/0" if sid.startswith("S") else "0/1")
                              for sid in arrangement}
                   for j in range(5)}
        m = compute_lrr_metrics(
            lrr=lrr, arrangement_by_sample=arrangement,
            marker_ids_in_lrr=[f"M{j}" for j in range(5)],
            genotype_matrix=gmatrix,
        )
        self.assertEqual(m.n_hom1, 0)
        self.assertEqual(m.missing_hom_class, "HOM_INV")
        # Het excess: H_obs = 20/25 = 0.8; p = 20/(2*25) = 0.4
        # H_exp = 2*0.4*0.6 = 0.48; H_obs/H_exp ≈ 1.67
        self.assertAlmostEqual(m.heterokaryotype_enrichment, 0.8 / 0.48, places=4)
        # FIS = 1 - 0.8/0.48 ≈ -0.667
        self.assertAlmostEqual(m.FIS, 1 - 0.8 / 0.48, places=4)


# ----------------------------------------------------------------------
# Linear regression + cross-LRR fit.
# ----------------------------------------------------------------------


class TestLinearRegression(unittest.TestCase):
    def test_perfect_positive_slope(self):
        xs = [1, 2, 3, 4, 5]
        ys = [2, 4, 6, 8, 10]
        slope, intercept, r, r2, t, p = linear_regression(xs, ys)
        self.assertAlmostEqual(slope, 2.0, places=6)
        self.assertAlmostEqual(intercept, 0.0, places=6)
        self.assertAlmostEqual(r, 1.0, places=6)

    def test_perfect_negative_slope(self):
        xs = [1, 2, 3, 4, 5]
        ys = [10, 8, 6, 4, 2]
        slope, intercept, r, r2, t, p = linear_regression(xs, ys)
        self.assertAlmostEqual(slope, -2.0, places=6)
        self.assertAlmostEqual(r, -1.0, places=6)

    def test_too_few_points(self):
        slope, *_ = linear_regression([1, 2], [3, 4])
        self.assertIsNone(slope)


class TestCrossLrrRegression(unittest.TestCase):
    def test_hypothesis_signal_recovered(self):
        """Synthesise 10 LRRs with monotonically increasing dXY and
        monotonically decreasing log(min_hom_het_ratio). The regression
        should recover a negative slope with r close to -1."""
        ms = []
        for i in range(10):
            dxy = 0.05 * (i + 1)
            min_ratio = math.exp(-dxy * 10)  # decreasing
            ms.append(LRRComparativeMetrics(
                lrr_id=f"LRR_{i:02d}", chrom="Chr1",
                start=i * 1_000_000, end=(i + 1) * 1_000_000,
                length=1_000_000, n_samples_called=20,
                n_hom0=5, n_het=10, n_hom1=5,
                freq_hom0=0.25, freq_het=0.5, freq_hom1=0.25,
                p_allele_1=0.5, H_obs=0.5, H_exp=0.5, FIS=0.0,
                heterokaryotype_enrichment=1.0,
                hom0_het_ratio=min_ratio, hom1_het_ratio=min_ratio,
                min_hom_het_ratio=min_ratio,
                missing_hom_class=None,
                n_markers_for_dxy=10,
                dxy_between_arrangements=dxy,
            ))
        reg = fit_comparative_regression(ms, y_kind="log_min_hom_het_ratio")
        self.assertEqual(reg.n_lrrs, 10)
        self.assertLess(reg.slope, 0)
        self.assertLess(reg.r, -0.99)
        # asymptotic caveat must be on, n < 30
        self.assertIn("asymptotic", reg.asymptotic_caveat)

    def test_no_signal_yields_flat_slope(self):
        # Random monotonic dXY but constant log(min_hom_het_ratio)
        ms = []
        for i in range(10):
            dxy = 0.05 * (i + 1)
            ms.append(LRRComparativeMetrics(
                lrr_id=f"LRR_{i:02d}", chrom="Chr1",
                start=0, end=1, length=1, n_samples_called=20,
                n_hom0=5, n_het=10, n_hom1=5,
                freq_hom0=0.25, freq_het=0.5, freq_hom1=0.25,
                p_allele_1=0.5, H_obs=0.5, H_exp=0.5, FIS=0.0,
                heterokaryotype_enrichment=1.0,
                hom0_het_ratio=1.0, hom1_het_ratio=1.0,
                min_hom_het_ratio=1.0, missing_hom_class=None,
                n_markers_for_dxy=10,
                dxy_between_arrangements=dxy,
            ))
        reg = fit_comparative_regression(ms, y_kind="log_min_hom_het_ratio")
        # log(1.0) = 0 for every LRR → zero variance in y → r undefined
        self.assertIsNone(reg.r)

    def test_fis_response_variable(self):
        ms = []
        for i in range(10):
            dxy = 0.05 * (i + 1)
            fis = -0.1 * (i + 1)  # increasingly negative FIS with dXY
            ms.append(LRRComparativeMetrics(
                lrr_id=f"LRR_{i:02d}", chrom="Chr1",
                start=0, end=1, length=1, n_samples_called=20,
                n_hom0=5, n_het=10, n_hom1=5,
                freq_hom0=0.25, freq_het=0.5, freq_hom1=0.25,
                p_allele_1=0.5, H_obs=0.5, H_exp=0.5, FIS=fis,
                heterokaryotype_enrichment=1.0,
                hom0_het_ratio=1.0, hom1_het_ratio=1.0,
                min_hom_het_ratio=1.0, missing_hom_class=None,
                n_markers_for_dxy=10,
                dxy_between_arrangements=dxy,
            ))
        reg = fit_comparative_regression(ms, y_kind="fis")
        self.assertLess(reg.slope, 0)   # FIS decreases (becomes more negative) with dXY
        self.assertLess(reg.r, -0.99)
