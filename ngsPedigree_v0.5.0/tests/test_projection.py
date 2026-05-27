"""
test_projection — MVP 2 dyad and triad projection into Table A.

Exercises:
  - composite_confidence enumeration
  - dyad projection on synthetic_dyad fixture at tier T1
  - dyad projection at tier T3 picks up the extra missense
  - triad projection on synthetic_triad fixture at tier T1
  - unphased-het row → unassigned / unresolved
  - Bronze segment downgrades confidence to 'low'
  - TableARow column list matches the Table A schema
  - TSV round-trip via io.write_tsv
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp import io  # noqa: E402
from hpp.parental_haps import build_for_parents, build_parental_hap_variants  # noqa: E402
from hpp.project import (  # noqa: E402
    TABLE_A_COLUMNS,
    TableARow,
    composite_confidence,
    project_dyad_to_offspring,
    project_triad_to_offspring,
)
from hpp.stage3_placeholder import load_dyad_map, load_parent_phase, load_triad_map  # noqa: E402
from hpp.variant_master import PlaceholderVariantMaster, VariantAnnotation  # noqa: E402
from hpp.vcf_lite import read_vcf  # noqa: E402

DYAD_FIX = THIS_DIR / "fixtures" / "synthetic_dyad"
TRIAD_FIX = THIS_DIR / "fixtures" / "synthetic_triad"
SCHEMAS = THIS_DIR.parent / "schemas"


# ----------------------------------------------------------------------
# Fixture helpers.
# ----------------------------------------------------------------------


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


def _row_by_variant(rows, variant_id, hap_copy=None):
    for r in rows:
        if r.variant_id == variant_id and (hap_copy is None or r.hap_copy == hap_copy):
            return r
    return None


# ----------------------------------------------------------------------
# Confidence rules — SPEC §5.
# ----------------------------------------------------------------------


class TestCompositeConfidence(unittest.TestCase):
    def test_bronze_always_low(self):
        for src in ("parent_homozygous", "parent_heterozygous_phased",
                    "parent_heterozygous_unphased"):
            self.assertEqual(composite_confidence(src, "Bronze"), "low")

    def test_homozygous_is_high_outside_bronze(self):
        self.assertEqual(composite_confidence("parent_homozygous", "Gold"), "high")
        self.assertEqual(composite_confidence("parent_homozygous", "Silver"), "high")

    def test_phased_het_tracks_segment(self):
        self.assertEqual(
            composite_confidence("parent_heterozygous_phased", "Gold"), "high")
        self.assertEqual(
            composite_confidence("parent_heterozygous_phased", "Silver"), "medium")

    def test_unphased_het_is_unresolved(self):
        self.assertEqual(
            composite_confidence("parent_heterozygous_unphased", "Gold"),
            "unresolved")
        self.assertEqual(
            composite_confidence("parent_heterozygous_unphased", "Silver"),
            "unresolved")


# ----------------------------------------------------------------------
# Dyad projection.
# ----------------------------------------------------------------------


class TestDyadProjection(unittest.TestCase):
    def setUp(self):
        self.segments = load_dyad_map(DYAD_FIX / "inheritance_map_dyad.tsv")
        self.phase = load_parent_phase(DYAD_FIX / "parent_phase.tsv")
        self.variants = list(read_vcf(DYAD_FIX / "joint.vcf"))
        self.parent_haps = build_parental_hap_variants(
            "P_A", self.variants, self.phase
        )
        self.vm = _load_variant_master_tsv(DYAD_FIX / "variant_master.tsv")

    def test_T1_produces_four_rows(self):
        rows = project_dyad_to_offspring(
            segments=self.segments,
            parent_haps=self.parent_haps,
            all_variants=self.variants,
            variant_master=self.vm,
            damaging_tier="T1",
        )
        self.assertEqual(len(rows), 4)

    def test_T1_missense_excluded(self):
        rows = project_dyad_to_offspring(
            segments=self.segments,
            parent_haps=self.parent_haps,
            all_variants=self.variants,
            variant_master=self.vm,
            damaging_tier="T1",
        )
        ids = {r.variant_id for r in rows}
        self.assertNotIn("LG01:200:G:C", ids)
        self.assertNotIn("LG01:500:A:C", ids)  # unannotated + parent 0/0

    def test_T3_includes_missense(self):
        rows = project_dyad_to_offspring(
            segments=self.segments,
            parent_haps=self.parent_haps,
            all_variants=self.variants,
            variant_master=self.vm,
            damaging_tier="T3",
        )
        # LG01:200 is missense + SIFT del + VESM<-7 → T3 damaging.
        # parent_phase puts it on hap 2; inherited hap is 1 in segment 1 →
        # row should be hap_from_P1, ref.
        row = _row_by_variant(rows, "LG01:200:G:C")
        self.assertIsNotNone(row)
        self.assertEqual(row.hap_copy, "hap_from_P1")
        self.assertEqual(row.allele_state, "ref")
        self.assertEqual(row.projection_source, "parent_heterozygous_phased")
        self.assertEqual(row.confidence, "high")

    def test_homozygous_row(self):
        rows = project_dyad_to_offspring(
            segments=self.segments, parent_haps=self.parent_haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        row = _row_by_variant(rows, "LG01:300:C:A")
        self.assertEqual(row.projection_source, "parent_homozygous")
        self.assertEqual(row.allele_state, "alt")
        self.assertEqual(row.hap_copy, "hap_from_P1")
        self.assertEqual(row.confidence, "high")

    def test_phased_het_on_inherited(self):
        rows = project_dyad_to_offspring(
            segments=self.segments, parent_haps=self.parent_haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        # LG01:100 — phase=hap 1, inherited=hap 1 → alt
        row = _row_by_variant(rows, "LG01:100:A:T")
        self.assertEqual(row.projection_source, "parent_heterozygous_phased")
        self.assertEqual(row.allele_state, "alt")
        self.assertEqual(row.confidence, "high")

    def test_unphased_het_is_unassigned(self):
        rows = project_dyad_to_offspring(
            segments=self.segments, parent_haps=self.parent_haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        # LG01:600 — parent is het, no parent_phase entry → unassigned
        row = _row_by_variant(rows, "LG01:600:G:T")
        self.assertEqual(row.hap_copy, "unassigned")
        self.assertEqual(row.allele_state, "unknown")
        self.assertEqual(row.projection_source, "parent_heterozygous_unphased")
        self.assertEqual(row.confidence, "unresolved")

    def test_relationship_fields(self):
        rows = project_dyad_to_offspring(
            segments=self.segments, parent_haps=self.parent_haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        for r in rows:
            self.assertEqual(r.relationship_id, "dyad_PA_OA")
            self.assertEqual(r.relationship_type, "dyad")
            self.assertEqual(r.offspring_sample_id, "O_A")
            self.assertEqual(r.chrom, "LG01")

    def test_only_hap_from_P1_for_dyad(self):
        rows = project_dyad_to_offspring(
            segments=self.segments, parent_haps=self.parent_haps,
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        haps = {r.hap_copy for r in rows}
        # dyad MVP 2: no hap_from_P2 rows (unknown parent side)
        self.assertNotIn("hap_from_P2", haps)


# ----------------------------------------------------------------------
# Triad projection.
# ----------------------------------------------------------------------


class TestTriadProjection(unittest.TestCase):
    def setUp(self):
        self.segments = load_triad_map(TRIAD_FIX / "inheritance_map_triad.tsv")
        self.phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        self.variants = list(read_vcf(TRIAD_FIX / "joint.vcf"))
        self.per_parent = build_for_parents(
            ["P_pat", "P_mat"], self.variants, self.phase
        )
        self.vm = _load_variant_master_tsv(TRIAD_FIX / "variant_master.tsv")
        self.rows = project_triad_to_offspring(
            segments=self.segments,
            paternal_haps=self.per_parent["P_pat"],
            maternal_haps=self.per_parent["P_mat"],
            all_variants=self.variants,
            variant_master=self.vm,
            damaging_tier="T1",
        )

    def test_T1_produces_five_rows(self):
        self.assertEqual(len(self.rows), 5)

    def test_paternal_split(self):
        paternal = [r for r in self.rows if r.hap_copy == "hap_from_P1"]
        # P_pat contributes LG01:100 (alt, phased→inherited) and
        # LG01:500 (ref, phased→other-hap). LG01:800 P_pat is 0/0 → skipped.
        self.assertEqual(len(paternal), 2)
        ids = {r.variant_id: r for r in paternal}
        self.assertEqual(ids["LG01:100:A:T"].allele_state, "alt")
        self.assertEqual(ids["LG01:500:C:G"].allele_state, "ref")

    def test_maternal_split(self):
        maternal = [r for r in self.rows if r.hap_copy == "hap_from_P2"]
        # P_mat contributes LG01:100 (alt, phased→inherited),
        # LG01:500 (alt, hom), LG01:800 (ref, phased→other-hap).
        self.assertEqual(len(maternal), 3)
        ids = {r.variant_id: r for r in maternal}
        self.assertEqual(ids["LG01:500:C:G"].projection_source, "parent_homozygous")
        self.assertEqual(ids["LG01:800:T:A"].allele_state, "ref")

    def test_T1_missense_excluded(self):
        self.assertFalse(any(r.variant_id == "LG01:300:G:C" for r in self.rows))

    def test_relationship_fields(self):
        for r in self.rows:
            self.assertEqual(r.relationship_id, "triad_T1")
            self.assertEqual(r.relationship_type, "triad")
            self.assertEqual(r.offspring_sample_id, "O_T1")


# ----------------------------------------------------------------------
# Bronze-segment downgrade.
# ----------------------------------------------------------------------


class TestBronzeDowngrade(unittest.TestCase):
    def test_bronze_segment_emits_low_confidence(self):
        # Build a synthetic Bronze-only triad segment that overlaps a
        # variant the maternal parent is hom-alt on. The composite
        # confidence rule says Bronze always → low.
        from hpp.stage3_placeholder import TriadSegment

        seg = TriadSegment(
            triad_id="t_bronze",
            paternal_sample_id="P_pat",
            maternal_sample_id="P_mat",
            offspring_sample_id="O_T1",
            chrom="LG01",
            seg_start=0, seg_end=1000,
            paternal_hap_inherited="1",
            maternal_hap_inherited="2",
            segment_confidence="Bronze",
            recomb_event_left=False, recomb_event_right=False,
            n_informative_markers=1,
        )
        variants = list(read_vcf(TRIAD_FIX / "joint.vcf"))
        phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        per_parent = build_for_parents(["P_pat", "P_mat"], variants, phase)
        vm = _load_variant_master_tsv(TRIAD_FIX / "variant_master.tsv")
        rows = project_triad_to_offspring(
            segments=[seg],
            paternal_haps=per_parent["P_pat"],
            maternal_haps=per_parent["P_mat"],
            all_variants=variants,
            variant_master=vm,
            damaging_tier="T1",
        )
        self.assertTrue(rows)
        # Every row must have confidence='low' because Bronze.
        self.assertTrue(all(r.confidence == "low" for r in rows))


# ----------------------------------------------------------------------
# Table A column contract.
# ----------------------------------------------------------------------


class TestTableASchemaAlignment(unittest.TestCase):
    def test_columns_match_schema(self):
        schema = json.loads((SCHEMAS / "A_hpp_offspring_haplotype_variants.schema.json").read_text())
        schema_cols = schema["properties"]["columns"]["items"]["enum"]
        self.assertEqual(TABLE_A_COLUMNS, schema_cols)

    def test_tsv_round_trip(self):
        # smoke: emit a row through io.write_tsv and check the header line.
        row = TableARow(
            offspring_sample_id="O", relationship_id="d", relationship_type="dyad",
            variant_id="X:1:A:T", chrom="X", pos=1, ref="A", alt="T",
            gene_id="G", transcript_id="T", consequence="stop_gained",
            impact="HIGH", sift_class=None, vesm_llr=None, splice_subclass=None,
            hap_copy="hap_from_P1", allele_state="alt",
            projection_source="parent_homozygous",
            segment_confidence="Gold", confidence="high",
        )
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            path = Path(fh.name)
        try:
            io.write_tsv(path, columns=TABLE_A_COLUMNS, rows=[row.to_dict()])
            lines = path.read_text().splitlines()
            self.assertEqual(lines[0], "\t".join(TABLE_A_COLUMNS))
            self.assertIn("X:1:A:T", lines[1])
            self.assertIn("hap_from_P1", lines[1])
        finally:
            path.unlink()


if __name__ == "__main__":
    unittest.main()
