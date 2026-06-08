"""
test_pedigree_extract — verify the Stage 1/2 → polarization-IN converter
closes the backbone from ngsRelate to ngsTracts.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.pedigree_extract import (  # noqa: E402
    PedigreeExtractError,
    build_polarization_input,
    extract_dyads,
    extract_triads,
    load_karyotype_tsv,
    load_stage1_pairwise,
    load_stage1_roster,
    load_stage2_review_set,
)

FIX = THIS_DIR / "fixtures" / "synthetic_pedigree_pipeline"
INV_FIX = THIS_DIR / "fixtures" / "synthetic_inversion"


# ----------------------------------------------------------------------
# Loaders.
# ----------------------------------------------------------------------


class TestLoaders(unittest.TestCase):
    def test_pairwise_loads(self):
        rows = load_stage1_pairwise(FIX / "stage1_pairwise.tsv")
        self.assertEqual(len(rows), 20)
        po = [r for r in rows if r["edge_class"] == "parent_offspring"]
        self.assertEqual(len(po), 8)

    def test_roster_loads(self):
        rost = load_stage1_roster(FIX / "stage1_roster.tsv")
        self.assertEqual(len(rost), 8)
        self.assertEqual(rost["S001"].possible_role, "father")
        self.assertEqual(rost["S002"].possible_role, "mother")
        self.assertEqual(rost["S004"].possible_role, "parent_a")
        self.assertEqual(rost["S003"].hub_type, "two_parents_with_sibship")

    def test_review_set_loads(self):
        rev = load_stage2_review_set(FIX / "stage2_pairwise.tsv")
        # only S004↔S007 is REVIEW
        self.assertEqual(rev, {("S004", "S007")})

    def test_karyotype_filters_by_inversion(self):
        calls = load_karyotype_tsv(FIX / "karyotype.tsv", "inv_LG01_pod")
        self.assertEqual(len(calls), 8)
        self.assertNotIn("S099", [c.sample_id for c in calls])

    def test_karyotype_missing_required_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("sample_id\nX\n")
            p = Path(fh.name)
        try:
            with self.assertRaises(PedigreeExtractError):
                load_karyotype_tsv(p, "inv_X")
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Dyad extraction.
# ----------------------------------------------------------------------


class TestExtractDyads(unittest.TestCase):
    def setUp(self):
        self.pairwise = load_stage1_pairwise(FIX / "stage1_pairwise.tsv")
        self.roster = load_stage1_roster(FIX / "stage1_roster.tsv")

    def test_eight_dyads_without_review(self):
        dyads = extract_dyads(pairwise=self.pairwise, roster=self.roster)
        self.assertEqual(len(dyads), 8)

    def test_sex_propagated_from_mother_father_roles(self):
        dyads = extract_dyads(pairwise=self.pairwise, roster=self.roster)
        # H001 hub: S001 (father) → S003 / S002 (mother) → S003
        d1 = [d for d in dyads if d.parent_sample_id == "S001"][0]
        d2 = [d for d in dyads if d.parent_sample_id == "S002"][0]
        self.assertEqual(d1.parent_sex, "male")
        self.assertEqual(d2.parent_sex, "female")

    def test_blind_parents_get_no_sex(self):
        dyads = extract_dyads(pairwise=self.pairwise, roster=self.roster)
        for d in dyads:
            if d.parent_sample_id in ("S004", "S005"):
                self.assertIsNone(d.parent_sex)

    def test_offspring_correctly_oriented(self):
        dyads = extract_dyads(pairwise=self.pairwise, roster=self.roster)
        for d in dyads:
            self.assertIn(d.parent_sample_id, {"S001", "S002", "S004", "S005"})
            self.assertIn(d.offspring_sample_id, {"S003", "S006", "S007", "S008"})

    def test_review_set_drops_flagged_edge(self):
        rev = load_stage2_review_set(FIX / "stage2_pairwise.tsv")
        dyads = extract_dyads(
            pairwise=self.pairwise, roster=self.roster, review_set=rev,
        )
        self.assertEqual(len(dyads), 7)
        # S004→S007 dyad is dropped; S005→S007 remains.
        pairs = {(d.parent_sample_id, d.offspring_sample_id) for d in dyads}
        self.assertNotIn(("S004", "S007"), pairs)
        self.assertIn(("S005", "S007"), pairs)


# ----------------------------------------------------------------------
# Triad extraction.
# ----------------------------------------------------------------------


class TestExtractTriads(unittest.TestCase):
    def setUp(self):
        self.pairwise = load_stage1_pairwise(FIX / "stage1_pairwise.tsv")
        self.roster = load_stage1_roster(FIX / "stage1_roster.tsv")

    def test_four_triads_without_review(self):
        triads, warnings = extract_triads(
            pairwise=self.pairwise, roster=self.roster,
        )
        self.assertEqual(len(triads), 4)
        # H001 with mother/father assignment.
        t_h001 = [t for t in triads if t.offspring_sample_id == "S003"][0]
        self.assertEqual(t_h001.paternal_sample_id, "S001")
        self.assertEqual(t_h001.maternal_sample_id, "S002")
        # H002 with blind parent_a/parent_b convention.
        h002_pat = {t.paternal_sample_id for t in triads
                    if t.offspring_sample_id in ("S006", "S007", "S008")}
        h002_mat = {t.maternal_sample_id for t in triads
                    if t.offspring_sample_id in ("S006", "S007", "S008")}
        self.assertEqual(h002_pat, {"S004"})
        self.assertEqual(h002_mat, {"S005"})

    def test_blind_mode_warning_emitted(self):
        _, warnings = extract_triads(pairwise=self.pairwise, roster=self.roster)
        joined = " | ".join(warnings)
        self.assertIn("sex unknown", joined)
        self.assertIn("H002", joined)

    def test_review_set_drops_triad_with_missing_PO_support(self):
        rev = load_stage2_review_set(FIX / "stage2_pairwise.tsv")
        triads, _ = extract_triads(
            pairwise=self.pairwise, roster=self.roster, review_set=rev,
        )
        # S007 triad needs PO from BOTH S004 and S005; S004↔S007 is
        # flagged → triad drops.
        self.assertEqual(len(triads), 3)
        off = {t.offspring_sample_id for t in triads}
        self.assertNotIn("S007", off)
        self.assertEqual(off, {"S003", "S006", "S008"})


# ----------------------------------------------------------------------
# Bundle builder.
# ----------------------------------------------------------------------


class TestBuildBundle(unittest.TestCase):
    def test_bundle_in_json_shape(self):
        b = build_polarization_input(
            stage1_edges_path=FIX / "stage1_pairwise.tsv",
            stage1_roster_path=FIX / "stage1_roster.tsv",
            karyotype_path=FIX / "karyotype.tsv",
            inversion_id="inv_LG01_pod",
            polarity_hint="band_0_is_REF",
        )
        doc = b.to_in_json()
        self.assertEqual(doc["schema"], "ngspedigree_karyotype_calls_in_v1")
        self.assertEqual(doc["inversion_id"], "inv_LG01_pod")
        self.assertEqual(doc["polarity_hint"], "band_0_is_REF")
        self.assertEqual(len(doc["karyotype_calls"]), 8)
        self.assertEqual(len(doc["dyads"]), 8)
        self.assertEqual(len(doc["triads"]), 4)

    def test_bad_polarity_hint_raises(self):
        with self.assertRaises(PedigreeExtractError):
            build_polarization_input(
                stage1_edges_path=FIX / "stage1_pairwise.tsv",
                stage1_roster_path=FIX / "stage1_roster.tsv",
                karyotype_path=FIX / "karyotype.tsv",
                inversion_id="inv_LG01_pod",
                polarity_hint="band_0_is_NEITHER",
            )


# ----------------------------------------------------------------------
# End-to-end: 06 → 05 CLI chain (backbone verification).
# ----------------------------------------------------------------------


class TestEndToEndPipeline(unittest.TestCase):
    def _run_chain(self, *, with_stage2=False, with_mtdna=False):
        scripts = THIS_DIR.parent / "scripts"
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            in_json = Path(fh.name)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            out_json = Path(fh.name)
        try:
            cmd_06 = [
                sys.executable, str(scripts / "06_build_polarization_input.py"),
                "--stage1-edges",  str(FIX / "stage1_pairwise.tsv"),
                "--stage1-roster", str(FIX / "stage1_roster.tsv"),
                "--karyotype",     str(FIX / "karyotype.tsv"),
                "--inversion-id",  "inv_LG01_pod",
                "--polarity-hint", "band_0_is_REF",
                "--out",           str(in_json),
            ]
            if with_stage2:
                cmd_06 += ["--stage2-edges", str(FIX / "stage2_pairwise.tsv")]
            subprocess.run(cmd_06, capture_output=True, text=True, check=True)
            cmd_05 = [
                sys.executable, str(scripts / "05_polarize_inversion.py"),
                "--in",  str(in_json),
                "--out", str(out_json),
            ]
            if with_mtdna:
                # Reuse the synthetic_inversion mtDNA fixture — it covers
                # S001..S008 with the same haplotype convention.
                cmd_05 += ["--mtdna", str(INV_FIX / "mtdna_haplotypes.tsv")]
            subprocess.run(cmd_05, capture_output=True, text=True, check=True)
            return json.loads(out_json.read_text())
        finally:
            in_json.unlink(missing_ok=True)
            out_json.unlink(missing_ok=True)

    def test_backbone_no_filters(self):
        doc = self._run_chain()
        self.assertEqual(doc["schema"],
                         "ngspedigree_polarized_transmissions_v1")
        self.assertEqual(doc["intended_consumer"], "ngsTracts")
        # 4 triads × 2 parents = 8 transmission rows (no bare dyads
        # survive triad-override).
        self.assertEqual(len(doc["transmissions"]), 8)
        self.assertEqual(doc["mtdna_validation"], {"supplied": False})

    def test_backbone_with_stage2_review(self):
        doc = self._run_chain(with_stage2=True)
        # Stage 2 drops S007 triad → 3 triads × 2 parents = 6 triad
        # transmissions. Bare dyad S005→S007 still survives (its edge is
        # not flagged) and S004→S007 is removed. So 6 triad rows + 1
        # dyad row = 7.
        n_triad = sum(1 for t in doc["transmissions"]
                      if t["relationship_type"] == "triad")
        n_dyad = sum(1 for t in doc["transmissions"]
                     if t["relationship_type"] == "dyad")
        self.assertEqual(n_triad, 6)
        self.assertEqual(n_dyad, 1)

    def test_backbone_with_mtdna_preflight(self):
        doc = self._run_chain(with_mtdna=True)
        mv = doc["mtdna_validation"]
        self.assertTrue(mv["supplied"])
        # S005→S007 mtDNA-incompatible (M_B vs M_C) → triad 3 (S007)
        # dropped pre-polarization, maternal dyad S005→S007 dropped.
        self.assertEqual(mv["n_incompatible"], 1)
        self.assertEqual(mv["n_triads_excluded"], 1)
