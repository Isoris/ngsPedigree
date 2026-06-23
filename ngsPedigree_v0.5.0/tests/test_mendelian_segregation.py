"""
test_mendelian_segregation — per-family × per-LRR cross-type
classification + goodness-of-fit (framework steps 7–8).
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.mendelian_segregation import (  # noqa: E402
    HET, HOM_INV, HOM_REF,
    analyze_family_lrr,
    binomial_two_sided_pvalue,
    chi2_gof_df2,
    cross_type,
    expected_distribution,
    is_informative_cross,
    summarise_lrr,
)


# ----------------------------------------------------------------------
# Cross-type table sanity.
# ----------------------------------------------------------------------


class TestCrossTable(unittest.TestCase):
    def test_canonical_cross_independent_of_order(self):
        self.assertEqual(cross_type(HOM_REF, HET), cross_type(HET, HOM_REF))
        self.assertEqual(cross_type(HOM_INV, HET), cross_type(HET, HOM_INV))

    def test_hom_x_hom_opposite_yields_het(self):
        self.assertEqual(expected_distribution(HOM_REF, HOM_INV),
                         {HET: 1.0})

    def test_het_x_het_yields_25_50_25(self):
        d = expected_distribution(HET, HET)
        self.assertAlmostEqual(d[HOM_REF], 0.25)
        self.assertAlmostEqual(d[HET], 0.5)
        self.assertAlmostEqual(d[HOM_INV], 0.25)

    def test_informativeness(self):
        self.assertFalse(is_informative_cross(HOM_REF, HOM_REF))
        self.assertFalse(is_informative_cross(HOM_INV, HOM_INV))
        self.assertFalse(is_informative_cross(HOM_REF, HOM_INV))  # fixed-het
        self.assertTrue(is_informative_cross(HOM_REF, HET))
        self.assertTrue(is_informative_cross(HET, HOM_INV))
        self.assertTrue(is_informative_cross(HET, HET))


# ----------------------------------------------------------------------
# Chi-square + binomial primitives.
# ----------------------------------------------------------------------


class TestStats(unittest.TestCase):
    def test_chi2_zero_when_obs_eq_exp(self):
        chi2, p = chi2_gof_df2([10, 20, 10], [10, 20, 10])
        self.assertAlmostEqual(chi2, 0.0)
        self.assertAlmostEqual(p, 1.0)

    def test_chi2_p_value_decreases_with_chi2(self):
        chi2a, pa = chi2_gof_df2([10, 20, 10], [10, 20, 10])
        chi2b, pb = chi2_gof_df2([0, 40, 0], [10, 20, 10])
        self.assertGreater(chi2b, chi2a)
        self.assertLess(pb, pa)

    def test_binomial_5050(self):
        self.assertAlmostEqual(binomial_two_sided_pvalue(5, 10, 0.5), 1.0)


# ----------------------------------------------------------------------
# Per-family per-LRR analysis paths.
# ----------------------------------------------------------------------


class TestAnalyzeFamilyLrr(unittest.TestCase):
    def test_fixed_inheritance_validation(self):
        # HOM_REF × HOM_INV → 100% HET expected; 5 offspring all HET.
        r = analyze_family_lrr(
            family_id="famA", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HOM_REF, maternal_arrangement=HOM_INV,
            offspring_arrangements=[HET] * 5,
        )
        self.assertEqual(r.interpretation, "fixed_inheritance_validation")
        self.assertEqual(r.test_used, "fixed_check")
        self.assertEqual(r.p_value, 1.0)

    def test_fixed_violation_pedigree_error(self):
        # HOM_REF × HOM_INV → should be 100% HET. A HOM_REF offspring is
        # a Mendelian violation; flag as fixed_violation.
        r = analyze_family_lrr(
            family_id="famB", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HOM_REF, maternal_arrangement=HOM_INV,
            offspring_arrangements=[HET, HET, HET, HOM_REF],
        )
        self.assertEqual(r.interpretation, "fixed_violation")
        self.assertEqual(r.p_value, 0.0)

    def test_segregation_consistent_het_x_hom(self):
        # HET × HOM_REF → expected 50/50. Observe 5 / 5.
        r = analyze_family_lrr(
            family_id="famC", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HET, maternal_arrangement=HOM_REF,
            offspring_arrangements=[HOM_REF] * 5 + [HET] * 5,
        )
        self.assertEqual(r.interpretation, "segregation_consistent")
        self.assertEqual(r.test_used, "binomial")

    def test_segregation_distorted_het_x_hom(self):
        # HET × HOM_REF → 50/50 expected. Observe 10/0 — extreme skew.
        r = analyze_family_lrr(
            family_id="famD", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HET, maternal_arrangement=HOM_REF,
            offspring_arrangements=[HOM_REF] * 10,
        )
        self.assertEqual(r.interpretation, "segregation_distorted")
        self.assertLess(r.p_value, 0.01)

    def test_het_x_het_homozygote_depletion(self):
        # HET × HET → expected 25/50/25 = 5/10/5 (for n=20). Observe
        # 0 / 20 / 0 — total homozygote depletion (pseudo-overdom signature).
        r = analyze_family_lrr(
            family_id="famE", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HET, maternal_arrangement=HET,
            offspring_arrangements=[HET] * 20,
        )
        self.assertEqual(r.interpretation, "homozygote_depletion")
        self.assertLess(r.p_value, 0.05)

    def test_het_x_het_consistent(self):
        # 5/10/5 observed under 5/10/5 expected → consistent.
        r = analyze_family_lrr(
            family_id="famF", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HET, maternal_arrangement=HET,
            offspring_arrangements=([HOM_REF] * 5 + [HET] * 10 + [HOM_INV] * 5),
        )
        self.assertEqual(r.interpretation, "segregation_consistent")

    def test_small_n_carve_out(self):
        r = analyze_family_lrr(
            family_id="famG", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement=HET, maternal_arrangement=HOM_REF,
            offspring_arrangements=[HET, HOM_REF],
            min_offspring_for_test=4,
        )
        self.assertEqual(r.interpretation, "small_n")

    def test_ambiguous_parent_ambiguous_record(self):
        r = analyze_family_lrr(
            family_id="famH", lrr_id="LRR_001",
            paternal_sample_id="P_F", maternal_sample_id="P_M",
            paternal_arrangement="uncalled",
            maternal_arrangement=HOM_REF,
            offspring_arrangements=[HET, HET],
        )
        self.assertEqual(r.interpretation, "ambiguous")


# ----------------------------------------------------------------------
# Cohort summary aggregation.
# ----------------------------------------------------------------------


class TestCohortSummary(unittest.TestCase):
    def test_aggregates_three_families(self):
        records = [
            analyze_family_lrr(
                family_id="A", lrr_id="L",
                paternal_sample_id="P_A_F", maternal_sample_id="P_A_M",
                paternal_arrangement=HET, maternal_arrangement=HET,
                offspring_arrangements=[HET] * 20,
            ),  # homozygote_depletion
            analyze_family_lrr(
                family_id="B", lrr_id="L",
                paternal_sample_id="P_B_F", maternal_sample_id="P_B_M",
                paternal_arrangement=HET, maternal_arrangement=HET,
                offspring_arrangements=([HOM_REF] * 5 + [HET] * 10
                                          + [HOM_INV] * 5),
            ),  # segregation_consistent
            analyze_family_lrr(
                family_id="C", lrr_id="L",
                paternal_sample_id="P_C_F", maternal_sample_id="P_C_M",
                paternal_arrangement=HOM_REF, maternal_arrangement=HOM_INV,
                offspring_arrangements=[HET] * 5,
            ),  # fixed_inheritance_validation
        ]
        s = summarise_lrr("L", records)
        self.assertEqual(s.n_families, 3)
        self.assertEqual(s.n_homozygote_depletion, 1)
        self.assertEqual(s.n_segregation_consistent, 1)
        self.assertEqual(s.n_fixed_inheritance_validation, 1)
        self.assertEqual(s.n_informative_families, 2)
        # HET × HET aggregate (across families A + B):
        # observed: 5 / 30 / 5 (out of 40)
        # expected under 25/50/25: 10/20/10
        # so this is distorted at the cohort level.
        self.assertEqual(s.het_het_n_total, 40)
        self.assertEqual(s.het_het_obs_HET, 30)
        self.assertIsNotNone(s.het_het_p_value)
