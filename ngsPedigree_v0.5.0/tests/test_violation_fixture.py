"""
test_violation_fixture — adversarial Mendelian-violation triad.

The synthetic_triad fixture is internally consistent and only exercises
the "0 inconsistent / 0 de novo / status=pass" path. This fixture
exercises the warn/fail thresholds and the de-novo flag for real.

See fixtures/synthetic_triad_violation/EXPECTED.md for the per-variant
truth table.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.gene_status import classify_gene_status  # noqa: E402
from hpp.parental_haps import build_for_parents  # noqa: E402
from hpp.project import project_triad_to_offspring  # noqa: E402
from hpp.stage3_placeholder import load_parent_phase, load_triad_map  # noqa: E402
from hpp.summary import summarise_triad  # noqa: E402
from hpp.variant_master import PlaceholderVariantMaster, VariantAnnotation  # noqa: E402
from hpp.vcf_lite import read_vcf  # noqa: E402

FIX = THIS_DIR / "fixtures" / "synthetic_triad_violation"


def _load_vm(path):
    vm = PlaceholderVariantMaster()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
            vm.add(VariantAnnotation(
                variant_id=row["variant_id"],
                gene_id=row.get("gene_id") or None,
                transcript_id=row.get("transcript_id") or None,
                consequence=row.get("consequence") or None,
                impact=row.get("impact") or None,
                sift_class=row.get("sift_class") or None,
                vesm_llr=float(row["vesm_llr"]) if row.get("vesm_llr") else None,
                splice_subclass=row.get("splice_subclass") or None,
            ))
    return vm


class TestViolationFixtureT1(unittest.TestCase):
    def setUp(self):
        self.segments = load_triad_map(FIX / "inheritance_map_triad.tsv")
        phase = load_parent_phase(FIX / "parent_phase.tsv")
        self.variants = list(read_vcf(FIX / "joint.vcf"))
        per_parent = build_for_parents(
            ["P_v_pat", "P_v_mat"], self.variants, phase
        )
        self.vm = _load_vm(FIX / "variant_master.tsv")
        self.table_a = project_triad_to_offspring(
            segments=self.segments,
            paternal_haps=per_parent["P_v_pat"],
            maternal_haps=per_parent["P_v_mat"],
            all_variants=self.variants, variant_master=self.vm,
            damaging_tier="T1",
        )
        self.table_b = classify_gene_status(
            table_a_rows=self.table_a, segments=self.segments,
            variant_master=self.vm, damaging_tier="T1",
        )
        self.row = summarise_triad(
            segments=self.segments, table_a_rows=self.table_a,
            table_b_rows=self.table_b, all_variants=self.variants,
            variant_master=self.vm, damaging_tier="T1",
        )

    def test_inconsistent_sites_is_three(self):
        # LG02:100 (de novo "1"), LG02:200 (1/1 from 0/0+0/1),
        # LG02:400 (0/0 from 1/1+0/0) → three inconsistent.
        self.assertEqual(self.row.mendelian_inconsistent_sites, 3)

    def test_inconsistent_damaging_at_T1_is_two(self):
        # LG02:100 and LG02:200 are T1 damaging (stop_gained / frameshift).
        # LG02:400 is T3-only missense → not counted at T1.
        self.assertEqual(self.row.mendelian_inconsistent_damaging_sites, 2)

    def test_de_novo_candidates_is_one(self):
        # Only LG02:100 has a novel allele (both parents 0/0, offspring 0/1).
        # LG02:200's "1" exists in P_v_mat. LG02:400's "0" exists in P_v_mat.
        self.assertEqual(self.row.n_de_novo_candidates, 1)

    def test_status_is_fail(self):
        # 3 inconsistencies ≥ 3 threshold → fail.
        self.assertEqual(self.row.mendelian_consistency_status, "fail")

    def test_n_damaging_in_parents(self):
        # T1 damaging in P_v_pat: none (LG02:100 is 0/0, LG02:200 is 0/0,
        # LG02:400 is missense T3 only). → 0.
        self.assertEqual(self.row.n_damaging_variants_in_paternal, 0)
        # T1 damaging in P_v_mat: LG02:200 (frameshift 0/1). → 1.
        self.assertEqual(self.row.n_damaging_variants_in_maternal, 1)


class TestViolationFixtureT3(unittest.TestCase):
    def setUp(self):
        segments = load_triad_map(FIX / "inheritance_map_triad.tsv")
        phase = load_parent_phase(FIX / "parent_phase.tsv")
        variants = list(read_vcf(FIX / "joint.vcf"))
        per_parent = build_for_parents(["P_v_pat", "P_v_mat"], variants, phase)
        vm = _load_vm(FIX / "variant_master.tsv")
        table_a = project_triad_to_offspring(
            segments=segments,
            paternal_haps=per_parent["P_v_pat"],
            maternal_haps=per_parent["P_v_mat"],
            all_variants=variants, variant_master=vm, damaging_tier="T3",
        )
        table_b = classify_gene_status(
            table_a_rows=table_a, segments=segments,
            variant_master=vm, damaging_tier="T3",
        )
        self.row = summarise_triad(
            segments=segments, table_a_rows=table_a, table_b_rows=table_b,
            all_variants=variants, variant_master=vm, damaging_tier="T3",
        )

    def test_T3_picks_up_missense_violation(self):
        # At T3, LG02:400 (missense, SIFT del, VESM<-7) is now damaging.
        # → inconsistent_damaging_sites rises to 3.
        self.assertEqual(self.row.mendelian_inconsistent_damaging_sites, 3)
        # Inconsistent_sites and de_novo unchanged.
        self.assertEqual(self.row.mendelian_inconsistent_sites, 3)
        self.assertEqual(self.row.n_de_novo_candidates, 1)
