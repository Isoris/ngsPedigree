"""
test_parental_haps — parental haplotype variant builder assertions.

Synthetic dyad and triad fixtures exercise every branch of
``build_parental_hap_variants``:

  - hom alt (1/1) → both haps
  - het with phase hap 1 → hap1
  - het with phase hap 2 → hap2
  - het with no phase → unphased
  - hom ref (0/0) → skipped
  - missing (./.) → skipped
"""

from __future__ import annotations

import io
import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.parental_haps import (  # noqa: E402
    build_for_parents,
    build_parental_hap_variants,
)
from hpp.stage3_placeholder import load_parent_phase  # noqa: E402
from hpp.vcf_lite import VariantRecord, read_vcf  # noqa: E402

DYAD_FIX = THIS_DIR / "fixtures" / "synthetic_dyad"
TRIAD_FIX = THIS_DIR / "fixtures" / "synthetic_triad"


def _ids(variants):
    return sorted(v.variant_id for v in variants)


class TestDyadBuild(unittest.TestCase):
    def setUp(self):
        self.variants = list(read_vcf(DYAD_FIX / "joint.vcf"))
        self.phase = load_parent_phase(DYAD_FIX / "parent_phase.tsv")

    def test_vcf_parses_six_variants(self):
        self.assertEqual(len(self.variants), 6)

    def test_normalised_gts(self):
        gts = {v.variant_id: v.genotypes["P_A"] for v in self.variants}
        self.assertEqual(gts["LG01:100:A:T"], "0/1")
        self.assertEqual(gts["LG01:300:C:A"], "1/1")
        self.assertEqual(gts["LG01:500:A:C"], "0/0")

    def test_hap1_has_phased_het_and_hom(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        self.assertEqual(
            _ids(rec.hap1),
            ["LG01:100:A:T", "LG01:300:C:A", "LG01:400:T:G"],
        )

    def test_hap2_has_phased_het_and_hom(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        self.assertEqual(_ids(rec.hap2), ["LG01:200:G:C", "LG01:300:C:A"])

    def test_unphased_has_one_het(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        self.assertEqual(_ids(rec.unphased), ["LG01:600:G:T"])

    def test_hom_ref_skipped(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        all_ids = _ids(rec.hap1) + _ids(rec.hap2) + _ids(rec.unphased)
        self.assertNotIn("LG01:500:A:C", all_ids)

    def test_n_total_counts_hom_twice(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        # hap1=3 (including hom), hap2=2 (including hom), unphased=1 → 6
        self.assertEqual(rec.n_total(), 6)

    def test_carries_on_hap(self):
        rec = build_parental_hap_variants("P_A", self.variants, self.phase)
        self.assertTrue(rec.carries_on_hap("LG01:100:A:T", hap_no=1))
        self.assertFalse(rec.carries_on_hap("LG01:100:A:T", hap_no=2))
        self.assertTrue(rec.carries_on_hap("LG01:300:C:A", hap_no=1))
        self.assertTrue(rec.carries_on_hap("LG01:300:C:A", hap_no=2))


class TestTriadBuild(unittest.TestCase):
    def setUp(self):
        self.variants = list(read_vcf(TRIAD_FIX / "joint.vcf"))
        self.phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        self.per_parent = build_for_parents(
            ["P_pat", "P_mat"], self.variants, self.phase
        )

    def test_paternal_haps(self):
        rec = self.per_parent["P_pat"]
        self.assertEqual(_ids(rec.hap1), ["LG01:100:A:T", "LG01:300:G:C"])
        self.assertEqual(_ids(rec.hap2), ["LG01:300:G:C", "LG01:500:C:G"])
        self.assertEqual(_ids(rec.unphased), [])

    def test_maternal_haps(self):
        rec = self.per_parent["P_mat"]
        self.assertEqual(_ids(rec.hap1), ["LG01:500:C:G", "LG01:800:T:A"])
        self.assertEqual(_ids(rec.hap2), ["LG01:100:A:T", "LG01:500:C:G"])
        self.assertEqual(_ids(rec.unphased), [])

    def test_both_parents_have_300_on_paternal_only(self):
        # P_pat is 1/1 at LG01:300 (both haps), P_mat is 0/0 (skipped)
        self.assertTrue(
            self.per_parent["P_pat"].carries_on_hap("LG01:300:G:C", 1)
        )
        self.assertTrue(
            self.per_parent["P_pat"].carries_on_hap("LG01:300:G:C", 2)
        )
        self.assertFalse(
            self.per_parent["P_mat"].carries_on_hap("LG01:300:G:C", 1)
        )
        self.assertFalse(
            self.per_parent["P_mat"].carries_on_hap("LG01:300:G:C", 2)
        )


class TestBuildEdgeCases(unittest.TestCase):
    def test_missing_genotype_skipped(self):
        v = VariantRecord(
            chrom="LG01",
            pos=1,
            ref="A",
            alt="T",
            genotypes={"P": "./."},
        )
        rec = build_parental_hap_variants("P", [v], parent_phase={})
        self.assertEqual(rec.n_total(), 0)

    def test_phased_vcf_separator_normalised(self):
        # vcf_lite should normalise 0|1 to 0/1
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP\n"
            "LG01\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0|1\n"
        )
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
            fh.write(vcf_text)
            tmp = Path(fh.name)
        try:
            v = list(read_vcf(tmp))[0]
            self.assertEqual(v.genotypes["P"], "0/1")
        finally:
            tmp.unlink()

    def test_unknown_genotype_raises(self):
        v = VariantRecord(
            chrom="LG01",
            pos=1,
            ref="A",
            alt="T",
            genotypes={"P": "2/2"},
        )
        with self.assertRaises(ValueError):
            build_parental_hap_variants("P", [v], parent_phase={})


if __name__ == "__main__":
    unittest.main()
