"""
test_karyotype_catalogue — registry-shape JSON adapter for
whole-genome karyotype calls.

Exercises:
  - HOM1/HET/HOM2 → band 0/1/2 mapping
  - per-LRR filtering with optional sample whitelist
  - inversion_id fallback to lrr_id
  - duplicate-sample detection
  - schema-violation paths
  - round-trip writer
  - hand-off to the polarization IN-JSON array shape
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.karyotype_catalogue import (  # noqa: E402
    CATALOGUE_SCHEMA_VERSION,
    BAND_TO_LABEL,
    KARYOTYPE_LABEL_TO_BAND,
    KaryotypeCatalogue,
    KaryotypeCatalogueError,
    catalogue_calls_to_in_json_array,
    load_catalogue,
    parse_catalogue,
    write_catalogue,
)

FIX = THIS_DIR / "fixtures" / "synthetic_catalogue"


class TestLabelMap(unittest.TestCase):
    def test_label_band_round_trip(self):
        for label, band in KARYOTYPE_LABEL_TO_BAND.items():
            self.assertEqual(BAND_TO_LABEL[band], label)
        self.assertEqual(KARYOTYPE_LABEL_TO_BAND, {"HOM1": 0, "HET": 1, "HOM2": 2})


class TestLoader(unittest.TestCase):
    def setUp(self):
        self.cat = load_catalogue(FIX / "karyotype_catalogue.json")

    def test_rows_loaded(self):
        self.assertEqual(len(self.cat.rows), 10)

    def test_lrr_inventory(self):
        self.assertEqual(self.cat.lrrs(), ["LRR_001", "LRR_002", "LRR_003"])

    def test_chrom_inventory(self):
        self.assertEqual(self.cat.chroms(), ["Chr1", "Chr2"])

    def test_sample_inventory(self):
        self.assertEqual(self.cat.samples(), ["S001", "S002", "S003", "S004"])

    def test_coverage_counts(self):
        cov = self.cat.coverage()
        self.assertEqual(cov, {"LRR_001": 4, "LRR_002": 3, "LRR_003": 3})


class TestFilterToInversion(unittest.TestCase):
    def setUp(self):
        self.cat = load_catalogue(FIX / "karyotype_catalogue.json")

    def test_lrr_001_yields_four_calls(self):
        calls = self.cat.filter_to_inversion("LRR_001")
        self.assertEqual(len(calls), 4)
        bands = {c.sample_id: c.band for c in calls}
        self.assertEqual(bands, {"S001": 0, "S002": 2, "S003": 1, "S004": 1})
        self.assertTrue(all(c.inversion_id == "LRR_001" for c in calls))

    def test_lrr_002_three_calls_with_medium_confidence(self):
        calls = self.cat.filter_to_inversion("LRR_002")
        self.assertEqual(len(calls), 3)
        s001 = next(c for c in calls if c.sample_id == "S001")
        self.assertEqual(s001.confidence, "medium")

    def test_inversion_id_override_picks_row_up(self):
        # LRR_003 has S002 tagged with inversion_id="inv_LG02_b"
        calls = self.cat.filter_to_inversion("inv_LG02_b")
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].sample_id, "S002")
        self.assertEqual(calls[0].band, 2)

    def test_sample_whitelist(self):
        calls = self.cat.filter_to_inversion(
            "LRR_001", sample_whitelist=["S001", "S003"],
        )
        self.assertEqual({c.sample_id for c in calls}, {"S001", "S003"})

    def test_duplicate_sample_raises(self):
        bad = {
            "schema": CATALOGUE_SCHEMA_VERSION,
            "rows": [
                {"chrom": "Chr1", "lrr_id": "LRR_X", "sample_id": "S", "karyotype": "HOM1"},
                {"chrom": "Chr1", "lrr_id": "LRR_X", "sample_id": "S", "karyotype": "HOM2"},
            ],
        }
        c = parse_catalogue(bad)
        with self.assertRaises(KaryotypeCatalogueError):
            c.filter_to_inversion("LRR_X")


class TestSchemaValidation(unittest.TestCase):
    def test_bad_schema_raises(self):
        with self.assertRaises(KaryotypeCatalogueError):
            parse_catalogue({"schema": "wrong_version", "rows": []})

    def test_missing_rows_raises(self):
        with self.assertRaises(KaryotypeCatalogueError):
            parse_catalogue({"schema": CATALOGUE_SCHEMA_VERSION})

    def test_bad_karyotype_label_raises(self):
        bad = {
            "schema": CATALOGUE_SCHEMA_VERSION,
            "rows": [
                {"chrom": "Chr1", "lrr_id": "L", "sample_id": "S", "karyotype": "HOMOZYGOUS"},
            ],
        }
        with self.assertRaises(KaryotypeCatalogueError):
            parse_catalogue(bad)

    def test_missing_required_field_raises(self):
        bad = {
            "schema": CATALOGUE_SCHEMA_VERSION,
            "rows": [{"chrom": "Chr1", "lrr_id": "L", "karyotype": "HET"}],
        }
        with self.assertRaises(KaryotypeCatalogueError):
            parse_catalogue(bad)


class TestRoundTrip(unittest.TestCase):
    def test_write_then_load(self):
        cat = load_catalogue(FIX / "karyotype_catalogue.json")
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            p = Path(fh.name)
        try:
            write_catalogue(p, cat)
            cat2 = load_catalogue(p)
            self.assertEqual(len(cat2.rows), len(cat.rows))
            for r1, r2 in zip(cat.rows, cat2.rows):
                self.assertEqual(r1.sample_id, r2.sample_id)
                self.assertEqual(r1.band, r2.band)
                self.assertEqual(r1.inversion_id, r2.inversion_id)
        finally:
            p.unlink()


class TestInJsonHandoff(unittest.TestCase):
    def test_array_shape_matches_polarization_in_json(self):
        cat = load_catalogue(FIX / "karyotype_catalogue.json")
        arr = catalogue_calls_to_in_json_array(cat, "LRR_001")
        self.assertEqual(len(arr), 4)
        for entry in arr:
            self.assertIn("sample_id", entry)
            self.assertIn("band", entry)
            self.assertIn("confidence", entry)
            self.assertIn(entry["band"], {0, 1, 2})
