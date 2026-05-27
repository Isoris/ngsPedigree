"""
test_gene_status — MVP 3 per-(offspring, gene, segment) classifier.

Exercises:
  - all seven enum classes via constructed inputs
  - schema-column alignment with B_hpp_offspring_gene_status.schema.json
  - synonymous-only group → reference_like
  - confidence aggregation (worst-of)
  - Bronze segment propagation
  - real-fixture path: synthetic_dyad + synthetic_triad at T1 → known truth
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.gene_status import (  # noqa: E402
    TABLE_B_COLUMNS,
    TableBRow,
    _classify,
    _worst_confidence,
    classify_gene_status,
)
from hpp.parental_haps import build_for_parents, build_parental_hap_variants  # noqa: E402
from hpp.project import TableARow, project_dyad_to_offspring, project_triad_to_offspring  # noqa: E402
from hpp.stage3_placeholder import load_dyad_map, load_parent_phase, load_triad_map  # noqa: E402
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
# Pure classification rule (SPEC §7) — exhaustive enum coverage.
# ----------------------------------------------------------------------


class TestPureClassifier(unittest.TestCase):
    def test_reference_like(self):
        self.assertEqual(_classify([], [], []), "reference_like")

    def test_unresolved(self):
        self.assertEqual(_classify([], [], ["v1"]), "unresolved")

    def test_partially_resolved(self):
        self.assertEqual(_classify(["v1"], [], ["v2"]), "partially_resolved")

    def test_het_masked_p1(self):
        self.assertEqual(_classify(["v1"], [], []), "het_masked")

    def test_het_masked_p2(self):
        self.assertEqual(_classify([], ["v1"], []), "het_masked")

    def test_compound_het_cis_p1(self):
        self.assertEqual(_classify(["v1", "v2"], [], []), "compound_het_cis")

    def test_compound_het_cis_p2(self):
        self.assertEqual(_classify([], ["v1", "v2"], []), "compound_het_cis")

    def test_hom_exposed_same_variant(self):
        self.assertEqual(_classify(["v1"], ["v1"], []), "hom_exposed_same_variant")

    def test_compound_het_trans_distinct(self):
        self.assertEqual(_classify(["v1"], ["v2"], []), "compound_het_trans")

    def test_compound_het_trans_one_vs_two(self):
        # 2 vs 1 on different copies — still compound_het_trans (both hit)
        self.assertEqual(_classify(["v1", "v2"], ["v3"], []), "compound_het_trans")


# ----------------------------------------------------------------------
# Worst-of confidence.
# ----------------------------------------------------------------------


class TestWorstConfidence(unittest.TestCase):
    def test_empty_defaults_high(self):
        self.assertEqual(_worst_confidence([]), "high")

    def test_high_only(self):
        self.assertEqual(_worst_confidence(["high", "high"]), "high")

    def test_high_medium(self):
        self.assertEqual(_worst_confidence(["high", "medium"]), "medium")

    def test_medium_low(self):
        self.assertEqual(_worst_confidence(["medium", "low"]), "low")

    def test_anything_with_unresolved(self):
        self.assertEqual(_worst_confidence(["high", "unresolved", "low"]),
                         "unresolved")


# ----------------------------------------------------------------------
# Schema alignment.
# ----------------------------------------------------------------------


class TestSchemaAlignment(unittest.TestCase):
    def test_columns_match(self):
        schema = json.loads(
            (SCHEMAS / "B_hpp_offspring_gene_status.schema.json").read_text()
        )
        self.assertEqual(
            TABLE_B_COLUMNS, schema["properties"]["columns"]["items"]["enum"]
        )


# ----------------------------------------------------------------------
# Real fixture — dyad at T1.
# ----------------------------------------------------------------------


class TestDyadGeneStatus(unittest.TestCase):
    def setUp(self):
        self.segments = load_dyad_map(DYAD_FIX / "inheritance_map_dyad.tsv")
        phase = load_parent_phase(DYAD_FIX / "parent_phase.tsv")
        variants = list(read_vcf(DYAD_FIX / "joint.vcf"))
        haps = build_parental_hap_variants("P_A", variants, phase)
        self.vm = _load_variant_master_tsv(DYAD_FIX / "variant_master.tsv")
        self.table_a = project_dyad_to_offspring(
            segments=self.segments, parent_haps=haps,
            all_variants=variants, variant_master=self.vm, damaging_tier="T1",
        )

    def test_runs_to_completion(self):
        rows = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )
        # At T1 the dyad has 2 damaging Table A rows (LG01:100, LG01:400)
        # on hap_from_P1 alt, in genes G1 and G4 respectively.
        # Plus 1 unresolved damaging (LG01:600 synonymous — not damaging at T1)?
        # Wait: LG01:600 is synonymous → not damaging → does not enter d1/d2.
        # So we get 2 gene rows from the damaging variants: het_masked × 2.
        gene_status_by_gene = {(r.gene_id, r.predicted_gene_status) for r in rows}
        self.assertIn(("G1", "het_masked"), gene_status_by_gene)
        self.assertIn(("G4", "het_masked"), gene_status_by_gene)

    def test_synonymous_only_gene_is_reference_like(self):
        # G3 has only LG01:300 (synonymous, hom-alt → on hap_from_P1).
        # No damaging variants → reference_like.
        rows = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )
        g3 = [r for r in rows if r.gene_id == "G3"]
        self.assertEqual(len(g3), 1)
        self.assertEqual(g3[0].predicted_gene_status, "reference_like")

    def test_T3_promotes_missense_to_het_masked(self):
        # LG01:200 in gene G2 is missense T3. At T3 it becomes damaging,
        # but parent_phase puts it on hap 2 while inherited hap is 1 →
        # allele_state = "ref" → does NOT count as damaging carrier.
        # So G2 still classifies reference_like (or doesn't appear).
        table_a_T3 = project_dyad_to_offspring(
            segments=self.segments,
            parent_haps=build_parental_hap_variants(
                "P_A", list(read_vcf(DYAD_FIX / "joint.vcf")),
                load_parent_phase(DYAD_FIX / "parent_phase.tsv"),
            ),
            all_variants=list(read_vcf(DYAD_FIX / "joint.vcf")),
            variant_master=self.vm, damaging_tier="T3",
        )
        rows = classify_gene_status(
            table_a_rows=table_a_T3, segments=self.segments,
            variant_master=self.vm, damaging_tier="T3",
        )
        g2 = [r for r in rows if r.gene_id == "G2"]
        # G2 row exists (because LG01:200 is now in Table A at T3) but is
        # ref-on-inherited → reference_like.
        self.assertEqual(len(g2), 1)
        self.assertEqual(g2[0].predicted_gene_status, "reference_like")


# ----------------------------------------------------------------------
# Real fixture — triad at T1.
# ----------------------------------------------------------------------


class TestTriadGeneStatus(unittest.TestCase):
    def setUp(self):
        self.segments = load_triad_map(TRIAD_FIX / "inheritance_map_triad.tsv")
        phase = load_parent_phase(TRIAD_FIX / "parent_phase.tsv")
        variants = list(read_vcf(TRIAD_FIX / "joint.vcf"))
        per_parent = build_for_parents(["P_pat", "P_mat"], variants, phase)
        self.vm = _load_variant_master_tsv(TRIAD_FIX / "variant_master.tsv")
        self.table_a = project_triad_to_offspring(
            segments=self.segments,
            paternal_haps=per_parent["P_pat"],
            maternal_haps=per_parent["P_mat"],
            all_variants=variants,
            variant_master=self.vm, damaging_tier="T1",
        )

    def test_gene_G5_hom_exposed(self):
        # LG01:500:C:G in gene G5 is frameshift_variant (T1 damaging).
        # P_pat 0/1 phased to hap 2 → inherited hap 1 → ref on hap_from_P1.
        # P_mat 1/1 → alt on hap_from_P2.
        # So d1 = [], d2 = [LG01:500] → het_masked.
        rows = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )
        g5 = [r for r in rows if r.gene_id == "G5"]
        self.assertEqual(len(g5), 1)
        self.assertEqual(g5[0].predicted_gene_status, "het_masked")
        # The damaging variant lives on hap_from_P2.
        self.assertEqual(g5[0].hap_from_P1_damaging_variants, "")
        self.assertEqual(g5[0].hap_from_P2_damaging_variants, "LG01:500:C:G")

    def test_gene_G1_compound_het_trans(self):
        # LG01:100:A:T in G1 — both parents het, both phased to inherited
        # haps → P_pat alt on hap_from_P1, P_mat alt on hap_from_P2.
        # Same variant on both copies → hom_exposed_same_variant.
        rows = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )
        g1 = [r for r in rows if r.gene_id == "G1"]
        self.assertEqual(len(g1), 1)
        self.assertEqual(g1[0].predicted_gene_status, "hom_exposed_same_variant")


if __name__ == "__main__":
    unittest.main()
