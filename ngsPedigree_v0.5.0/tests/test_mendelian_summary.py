"""
test_mendelian_summary — MVP 4 Mendelian checks + transmission summaries.

Exercises:
  - mendelian_expected_set covers Punnett-square enumeration
  - triad_consistent on consistent and inconsistent triples
  - offspring_carries_novel_allele de novo flag
  - dyad_partial_consistent on parent-hom and parent-het cases
  - status_from_counts pass/warn/fail/untestable boundaries
  - summarise_dyad and summarise_triad column alignment with schemas C, D
  - end-to-end run on synthetic_dyad and synthetic_triad fixtures
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.gene_status import classify_gene_status  # noqa: E402
from hpp.mendelian import (  # noqa: E402
    dyad_partial_consistent,
    is_called,
    mendelian_expected_set,
    offspring_carries_novel_allele,
    status_from_counts,
    triad_consistent,
)
from hpp.parental_haps import build_for_parents, build_parental_hap_variants  # noqa: E402
from hpp.project import project_dyad_to_offspring, project_triad_to_offspring  # noqa: E402
from hpp.stage3_placeholder import load_dyad_map, load_parent_phase, load_triad_map  # noqa: E402
from hpp.summary import (  # noqa: E402
    TABLE_C_COLUMNS,
    TABLE_D_COLUMNS,
    summarise_dyad,
    summarise_triad,
)
from hpp.variant_master import PlaceholderVariantMaster, VariantAnnotation  # noqa: E402
from hpp.vcf_lite import read_vcf  # noqa: E402

DYAD_FIX = THIS_DIR / "fixtures" / "synthetic_dyad"
TRIAD_FIX = THIS_DIR / "fixtures" / "synthetic_triad"
SCHEMAS = THIS_DIR.parent / "schemas"


def _load_variant_master_tsv(path: Path) -> PlaceholderVariantMaster:
    vm = PlaceholderVariantMaster()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
            vm.add(VariantAnnotation(
                variant_id=row["variant_id"],
                gene_id=row["gene_id"] or None,
                transcript_id=row["transcript_id"] or None,
                consequence=row["consequence"] or None,
                impact=row["impact"] or None,
                sift_class=row["sift_class"] or None,
                vesm_llr=float(row["vesm_llr"]) if row["vesm_llr"] else None,
                splice_subclass=row["splice_subclass"] or None,
            ))
    return vm


# ----------------------------------------------------------------------
# Mendelian primitives.
# ----------------------------------------------------------------------


class TestMendelianPrimitives(unittest.TestCase):
    def test_expected_set_het_x_het(self):
        self.assertEqual(
            mendelian_expected_set("0/1", "0/1"),
            {"0/0", "0/1", "1/1"},
        )

    def test_expected_set_hom_x_hom(self):
        self.assertEqual(mendelian_expected_set("0/0", "1/1"), {"0/1"})

    def test_expected_set_hom_x_het(self):
        self.assertEqual(mendelian_expected_set("1/1", "0/1"), {"0/1", "1/1"})

    def test_triad_consistent_true(self):
        self.assertTrue(triad_consistent("0/1", "0/1", "1/1"))

    def test_triad_consistent_false(self):
        self.assertFalse(triad_consistent("0/0", "0/0", "0/1"))

    def test_triad_consistent_uncalled_returns_true(self):
        # Missing GTs are not inconsistencies.
        self.assertTrue(triad_consistent("./.", "0/1", "0/0"))

    def test_de_novo_detected(self):
        # Both parents 0/0 but offspring carries 1 → novel allele.
        self.assertTrue(offspring_carries_novel_allele("0/0", "0/0", "0/1"))

    def test_de_novo_not_detected_when_parent_has_it(self):
        self.assertFalse(offspring_carries_novel_allele("0/1", "0/0", "0/1"))

    def test_dyad_partial_parent_hom(self):
        self.assertTrue(dyad_partial_consistent("1/1", "0/1"))
        self.assertFalse(dyad_partial_consistent("1/1", "0/0"))
        self.assertTrue(dyad_partial_consistent("0/0", "0/1"))
        self.assertFalse(dyad_partial_consistent("0/0", "1/1"))

    def test_dyad_partial_parent_het_is_untestable(self):
        self.assertIsNone(dyad_partial_consistent("0/1", "0/0"))

    def test_dyad_partial_missing_is_untestable(self):
        self.assertIsNone(dyad_partial_consistent("./.", "0/1"))


# ----------------------------------------------------------------------
# Status threshold rule.
# ----------------------------------------------------------------------


class TestStatusFromCounts(unittest.TestCase):
    def test_untestable_when_zero(self):
        self.assertEqual(
            status_from_counts(inconsistent=0, testable=0), "untestable"
        )

    def test_pass(self):
        self.assertEqual(
            status_from_counts(inconsistent=0, testable=10), "pass"
        )

    def test_warn(self):
        self.assertEqual(
            status_from_counts(inconsistent=1, testable=10), "warn"
        )
        self.assertEqual(
            status_from_counts(inconsistent=2, testable=10), "warn"
        )

    def test_fail(self):
        self.assertEqual(
            status_from_counts(inconsistent=3, testable=10), "fail"
        )
        self.assertEqual(
            status_from_counts(inconsistent=99, testable=100), "fail"
        )


# ----------------------------------------------------------------------
# Schema alignment.
# ----------------------------------------------------------------------


class TestSchemaAlignment(unittest.TestCase):
    def test_table_C_columns(self):
        schema = json.loads(
            (SCHEMAS / "C_hpp_dyad_transmission_summary.schema.json").read_text()
        )
        self.assertEqual(
            TABLE_C_COLUMNS, schema["properties"]["columns"]["items"]["enum"]
        )

    def test_table_D_columns(self):
        schema = json.loads(
            (SCHEMAS / "D_hpp_triad_transmission_summary.schema.json").read_text()
        )
        self.assertEqual(
            TABLE_D_COLUMNS, schema["properties"]["columns"]["items"]["enum"]
        )


# ----------------------------------------------------------------------
# End-to-end on fixtures.
# ----------------------------------------------------------------------


class TestSummariseDyad(unittest.TestCase):
    def setUp(self):
        self.segments = load_dyad_map(DYAD_FIX / "inheritance_map_dyad.tsv")
        phase = load_parent_phase(DYAD_FIX / "parent_phase.tsv")
        self.variants = list(read_vcf(DYAD_FIX / "joint.vcf"))
        haps = build_parental_hap_variants("P_A", self.variants, phase)
        self.vm = _load_variant_master_tsv(DYAD_FIX / "variant_master.tsv")
        self.table_a = project_dyad_to_offspring(
            segments=self.segments, parent_haps=haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        self.table_b = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )

    def test_returns_one_row(self):
        row = summarise_dyad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        self.assertEqual(row.dyad_id, "dyad_PA_OA")
        self.assertEqual(row.n_segments_total, 2)
        self.assertEqual(row.n_segments_Gold, 1)
        self.assertEqual(row.n_segments_Silver, 1)

    def test_n_damaging_in_parent(self):
        row = summarise_dyad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        # T1 damaging in P_A: LG01:100 (stop_gained, 0/1) and
        # LG01:400 (frameshift_variant, 0/1). LG01:300 is synonymous
        # (1/1 but not a LoF). So 2 T1 damaging in parent.
        self.assertEqual(row.n_damaging_variants_in_parent, 2)

    def test_partial_mendelian_pass(self):
        row = summarise_dyad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        # The fixture is internally consistent; partial Mendelian should pass.
        self.assertEqual(row.mendelian_consistency_status, "pass")
        self.assertEqual(row.mendelian_inconsistent_sites, 0)


class TestSummariseTriad(unittest.TestCase):
    def setUp(self):
        self.segments = load_triad_map(TRIAD_FIX / "inheritance_map_triad.tsv")
        phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        self.variants = list(read_vcf(TRIAD_FIX / "joint.vcf"))
        per_parent = build_for_parents(
            ["P_pat", "P_mat"], self.variants, phase
        )
        self.vm = _load_variant_master_tsv(TRIAD_FIX / "variant_master.tsv")
        self.table_a = project_triad_to_offspring(
            segments=self.segments,
            paternal_haps=per_parent["P_pat"],
            maternal_haps=per_parent["P_mat"],
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        self.table_b = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )

    def test_returns_one_row(self):
        row = summarise_triad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        self.assertEqual(row.triad_id, "triad_T1")

    def test_segment_confidence_counts(self):
        row = summarise_triad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        self.assertEqual(row.n_segments_Gold, 1)
        self.assertEqual(row.n_segments_Bronze, 1)
        self.assertEqual(row.n_segments_Silver, 0)

    def test_full_mendelian_status(self):
        row = summarise_triad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        # Synthetic fixture is internally consistent → pass.
        self.assertEqual(row.mendelian_consistency_status, "pass")
        self.assertEqual(row.mendelian_inconsistent_sites, 0)
        self.assertEqual(row.mendelian_inconsistent_damaging_sites, 0)
        self.assertEqual(row.n_de_novo_candidates, 0)

    def test_damaging_in_both_parents_counted(self):
        row = summarise_triad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )
        # P_pat T1 damaging: LG01:100 (0/1), LG01:500 (0/1) → 2.
        # P_mat T1 damaging: LG01:100 (0/1), LG01:500 (1/1) → 2.
        self.assertEqual(row.n_damaging_variants_in_paternal, 2)
        self.assertEqual(row.n_damaging_variants_in_maternal, 2)


if __name__ == "__main__":
    unittest.main()
