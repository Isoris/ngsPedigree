"""
test_mtdna_check — maternal-lineage validation layer.

Exercises:
  - TSV loader (required columns, optional columns, error paths)
  - pair check: compatible / incompatible / ambiguous
  - Hamming-distance soft-compatibility when sequences are present
  - pedigree-wide check + dedupe across triad/maternal-dyad overlap
  - filter step removes only incompatible relationships
  - integration with the polarization CLI on the synthetic_inversion fixture
    (with --mtdna): triad 3 excluded by mtDNA, triad 4 still excluded by
    nuclear; OUT JSON carries mtdna_validation block; transmission counts
    shift to 7 rows with S004→S007 promoted to a bare paternal dyad row.
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

from hpp.inversion_polarization import DyadPair, TriadTrio  # noqa: E402
from hpp.mtdna_check import (  # noqa: E402
    HAMMING_THRESHOLD_DEFAULT,
    MtdnaContractError,
    MtdnaRecord,
    build_not_supplied_block,
    build_validation_block,
    check_pair,
    check_pedigree,
    filter_by_mtdna,
    load_mtdna_haplotypes,
)

FIX = THIS_DIR / "fixtures" / "synthetic_inversion"


# ----------------------------------------------------------------------
# Loader.
# ----------------------------------------------------------------------


class TestLoader(unittest.TestCase):
    def test_loads_fixture(self):
        recs = load_mtdna_haplotypes(FIX / "mtdna_haplotypes.tsv")
        self.assertEqual(len(recs), 8)
        self.assertEqual(recs["S002"].haplotype, "M_A")
        self.assertEqual(recs["S007"].haplotype, "M_C")
        self.assertIsNone(recs["S002"].sequence)

    def test_loads_with_optional_sequence(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("sample_id\tmtdna_haplotype\tmtdna_sequence\tmtdna_n_sites\n")
            fh.write("X\tHA\tACGT\t4\n")
            fh.write("Y\tHB\tACCT\t4\n")
            p = Path(fh.name)
        try:
            recs = load_mtdna_haplotypes(p)
            self.assertEqual(recs["X"].sequence, "ACGT")
            self.assertEqual(recs["X"].n_sites, 4)
            self.assertEqual(recs["Y"].sequence, "ACCT")
        finally:
            p.unlink()

    def test_missing_required_column_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("sample_id\n")
            fh.write("X\n")
            p = Path(fh.name)
        try:
            with self.assertRaises(MtdnaContractError):
                load_mtdna_haplotypes(p)
        finally:
            p.unlink()

    def test_contradictory_record_raises(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("sample_id\tmtdna_haplotype\n")
            fh.write("X\tHA\n")
            fh.write("X\tHB\n")
            p = Path(fh.name)
        try:
            with self.assertRaises(MtdnaContractError):
                load_mtdna_haplotypes(p)
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# Pair check.
# ----------------------------------------------------------------------


class TestPairCheck(unittest.TestCase):
    def setUp(self):
        self.recs = {
            "M": MtdnaRecord("M", "HA"),
            "O_same": MtdnaRecord("O_same", "HA"),
            "O_diff": MtdnaRecord("O_diff", "HB"),
            "M_seq": MtdnaRecord("M_seq", "HA", sequence="ACGT"),
            "O_close": MtdnaRecord("O_close", "HB", sequence="ACCT"),
            "O_far": MtdnaRecord("O_far", "HB", sequence="TGCA"),
        }

    def test_compatible_on_label_match(self):
        c = check_pair("M", "O_same", self.recs)
        self.assertEqual(c.status, "compatible")
        self.assertEqual(c.distance, 0)

    def test_incompatible_on_label_diff_no_seq(self):
        c = check_pair("M", "O_diff", self.recs)
        self.assertEqual(c.status, "incompatible")
        self.assertIsNone(c.distance)

    def test_hamming_within_threshold_promotes_to_compatible(self):
        # HA "ACGT" vs HB "ACCT" → Hamming=1 ≤ default threshold 2.
        c = check_pair("M_seq", "O_close", self.recs)
        self.assertEqual(c.status, "compatible")
        self.assertEqual(c.distance, 1)
        self.assertIn("Hamming", c.note)

    def test_hamming_above_threshold_stays_incompatible(self):
        # ACGT vs TGCA → all 4 differ.
        c = check_pair("M_seq", "O_far", self.recs)
        self.assertEqual(c.status, "incompatible")
        self.assertEqual(c.distance, 4)

    def test_ambiguous_when_record_missing(self):
        c = check_pair("M", "S_missing", self.recs)
        self.assertEqual(c.status, "ambiguous")
        self.assertIn("S_missing", c.note)


# ----------------------------------------------------------------------
# Pedigree-wide check + filter.
# ----------------------------------------------------------------------


class TestPedigreeCheckAndFilter(unittest.TestCase):
    def setUp(self):
        self.recs = load_mtdna_haplotypes(FIX / "mtdna_haplotypes.tsv")
        self.dyads = [
            DyadPair("S001", "S003", "male"),
            DyadPair("S002", "S003", "female"),
            DyadPair("S004", "S006", "male"),
            DyadPair("S005", "S006", "female"),
            DyadPair("S004", "S007", "male"),
            DyadPair("S005", "S007", "female"),
            DyadPair("S004", "S008", "male"),
            DyadPair("S005", "S008", "female"),
        ]
        self.triads = [
            TriadTrio("S001", "S002", "S003"),
            TriadTrio("S004", "S005", "S006"),
            TriadTrio("S004", "S005", "S007"),
            TriadTrio("S004", "S005", "S008"),
        ]

    def test_check_counts(self):
        checks = check_pedigree(self.recs, self.dyads, self.triads)
        # 4 unique (mother, offspring) pairs after triad/maternal-dyad dedup.
        self.assertEqual(len(checks), 4)
        by_status = {c.status for c in checks}
        self.assertEqual(by_status, {"compatible", "incompatible"})
        # 1 incompatible: (S005, S007)
        incomp = [c for c in checks if c.status == "incompatible"]
        self.assertEqual(len(incomp), 1)
        self.assertEqual(incomp[0].mother_sample_id, "S005")
        self.assertEqual(incomp[0].offspring_sample_id, "S007")

    def test_filter_removes_incompatible_triad_and_maternal_dyad(self):
        checks = check_pedigree(self.recs, self.dyads, self.triads)
        f_dyads, f_triads = filter_by_mtdna(self.dyads, self.triads, checks)
        # Triad 3 (S007) removed.
        self.assertEqual(len(f_triads), 3)
        self.assertNotIn(
            ("S004", "S005", "S007"),
            [(t.paternal_sample_id, t.maternal_sample_id, t.offspring_sample_id)
             for t in f_triads],
        )
        # Maternal dyad (S005→S007) removed; paternal dyad (S004→S007) kept.
        d_pairs = [(d.parent_sample_id, d.offspring_sample_id, d.parent_sex)
                   for d in f_dyads]
        self.assertNotIn(("S005", "S007", "female"), d_pairs)
        self.assertIn(("S004", "S007", "male"), d_pairs)
        self.assertEqual(len(f_dyads), 7)

    def test_ambiguous_does_not_filter(self):
        # Missing record means "no evidence" — should not exclude.
        triads = [TriadTrio("S001", "S_unknown", "S003")]
        dyads = []
        checks = check_pedigree(self.recs, dyads, triads)
        self.assertEqual(len(checks), 1)
        self.assertEqual(checks[0].status, "ambiguous")
        f_dyads, f_triads = filter_by_mtdna(dyads, triads, checks)
        self.assertEqual(len(f_triads), 1)


# ----------------------------------------------------------------------
# Block builders.
# ----------------------------------------------------------------------


class TestValidationBlock(unittest.TestCase):
    def test_not_supplied(self):
        self.assertEqual(build_not_supplied_block(), {"supplied": False})

    def test_validation_block_shape(self):
        recs = load_mtdna_haplotypes(FIX / "mtdna_haplotypes.tsv")
        triads = [
            TriadTrio("S001", "S002", "S003"),
            TriadTrio("S004", "S005", "S007"),
        ]
        checks = check_pedigree(recs, [], triads)
        block = build_validation_block(
            checks, n_triads_excluded=1, n_dyads_excluded=0,
        )
        self.assertTrue(block["supplied"])
        self.assertEqual(block["n_relationships_checked"], 2)
        self.assertEqual(block["n_compatible"], 1)
        self.assertEqual(block["n_incompatible"], 1)
        self.assertEqual(block["n_triads_excluded"], 1)
        self.assertEqual(block["hamming_threshold"], HAMMING_THRESHOLD_DEFAULT)
        self.assertEqual(len(block["checks"]), 2)


# ----------------------------------------------------------------------
# End-to-end via the polarization CLI with --mtdna.
# ----------------------------------------------------------------------


class TestPolarizeCliWithMtdna(unittest.TestCase):
    def test_cli_emits_mtdna_validation_block(self):
        script = THIS_DIR.parent / "scripts" / "05_polarize_inversion.py"
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            out = Path(fh.name)
        try:
            res = subprocess.run(
                [sys.executable, str(script),
                 "--in",     str(FIX / "karyotype_calls.json"),
                 "--mtdna",  str(FIX / "mtdna_haplotypes.tsv"),
                 "--out",    str(out)],
                capture_output=True, text=True, check=True,
            )
            doc = json.loads(out.read_text())
            mv = doc["mtdna_validation"]
            self.assertTrue(mv["supplied"])
            self.assertEqual(mv["n_relationships_checked"], 4)
            self.assertEqual(mv["n_compatible"], 3)
            self.assertEqual(mv["n_incompatible"], 1)
            self.assertEqual(mv["n_triads_excluded"], 1)
            self.assertEqual(mv["n_dyads_excluded"], 1)
            # 7 transmissions after mtDNA-driven triad-3 removal: 6
            # triad-resolved + 1 paternal-dyad fallback (S004→S007).
            self.assertEqual(len(doc["transmissions"]), 7)
            # Paternal-side dyad fallback row exists with HET parent
            # transmitting INV (S007 is HOM_INV).
            fallback = [t for t in doc["transmissions"]
                        if t["relationship_type"] == "dyad"
                        and t["parent_sample_id"] == "S004"
                        and t["offspring_sample_id"] == "S007"]
            self.assertEqual(len(fallback), 1)
            self.assertEqual(fallback[0]["transmitted_arrangement"], "INV")
            self.assertTrue(fallback[0]["informative_for_drive"])
            # Drive: 1 REF + 1 INV informative.
            self.assertEqual(doc["drive_stats"]["n_informative_transmissions"], 2)
            self.assertEqual(doc["drive_stats"]["n_REF_transmitted"], 1)
            self.assertEqual(doc["drive_stats"]["n_INV_transmitted"], 1)
        finally:
            out.unlink()

    def test_cli_without_mtdna_records_not_supplied(self):
        script = THIS_DIR.parent / "scripts" / "05_polarize_inversion.py"
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            out = Path(fh.name)
        try:
            subprocess.run(
                [sys.executable, str(script),
                 "--in",  str(FIX / "karyotype_calls.json"),
                 "--out", str(out)],
                capture_output=True, text=True, check=True,
            )
            doc = json.loads(out.read_text())
            self.assertEqual(doc["mtdna_validation"], {"supplied": False})
            self.assertEqual(len(doc["transmissions"]), 8)  # original count
        finally:
            out.unlink()
