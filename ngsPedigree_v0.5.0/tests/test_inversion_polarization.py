"""
test_inversion_polarization — layers 2-3 (compatibility + transmission calling)
and the ngsTracts IN/OUT JSON adapters.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.inversion_polarization import (  # noqa: E402
    ARRANGEMENT_HET,
    ARRANGEMENT_HOM_INV,
    ARRANGEMENT_HOM_REF,
    DyadPair,
    KaryotypeCall,
    TriadTrio,
    _arrangement_for_band,
    _dyad_compatible,
    _infer_dyad_transmission,
    _infer_triad_transmission,
    _triad_compatible,
    binomial_two_sided_pvalue,
    call_transmissions,
    drive_test,
    polarize,
)
from hpp.ngstracts_io import (  # noqa: E402
    NgsTractsAdapterError,
    load_karyotype_calls,
    write_polarized_transmissions,
)

FIX = THIS_DIR / "fixtures" / "synthetic_inversion"
SCHEMAS = THIS_DIR.parent / "schemas"


# ----------------------------------------------------------------------
# Compatibility primitives.
# ----------------------------------------------------------------------


class TestArrangementForBand(unittest.TestCase):
    def test_band_1_always_het(self):
        self.assertEqual(_arrangement_for_band(1, "band_0_is_REF"), ARRANGEMENT_HET)
        self.assertEqual(_arrangement_for_band(1, "band_0_is_INV"), ARRANGEMENT_HET)

    def test_band_0_under_each_polarity(self):
        self.assertEqual(_arrangement_for_band(0, "band_0_is_REF"), ARRANGEMENT_HOM_REF)
        self.assertEqual(_arrangement_for_band(0, "band_0_is_INV"), ARRANGEMENT_HOM_INV)

    def test_band_2_under_each_polarity(self):
        self.assertEqual(_arrangement_for_band(2, "band_0_is_REF"), ARRANGEMENT_HOM_INV)
        self.assertEqual(_arrangement_for_band(2, "band_0_is_INV"), ARRANGEMENT_HOM_REF)

    def test_bad_polarity_raises(self):
        with self.assertRaises(ValueError):
            _arrangement_for_band(0, "band_0_is_NEITHER")


class TestDyadCompatible(unittest.TestCase):
    def test_hom_ref_parent(self):
        self.assertTrue(_dyad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_REF))
        self.assertTrue(_dyad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HET))
        self.assertFalse(_dyad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV))

    def test_hom_inv_parent(self):
        self.assertFalse(_dyad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_REF))
        self.assertTrue(_dyad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HET))
        self.assertTrue(_dyad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_INV))

    def test_het_parent_always_ok(self):
        for o in (ARRANGEMENT_HOM_REF, ARRANGEMENT_HET, ARRANGEMENT_HOM_INV):
            self.assertTrue(_dyad_compatible(ARRANGEMENT_HET, o))


class TestTriadCompatible(unittest.TestCase):
    def test_hom_x_hom_opposite(self):
        # HOM_REF × HOM_INV → only HET allowed.
        self.assertTrue(_triad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV, ARRANGEMENT_HET))
        self.assertFalse(_triad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_REF))
        self.assertFalse(_triad_compatible(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_INV))

    def test_het_x_het_allows_all(self):
        for o in (ARRANGEMENT_HOM_REF, ARRANGEMENT_HET, ARRANGEMENT_HOM_INV):
            self.assertTrue(_triad_compatible(ARRANGEMENT_HET, ARRANGEMENT_HET, o))

    def test_hom_inv_x_het_disallows_hom_ref(self):
        self.assertFalse(_triad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HET, ARRANGEMENT_HOM_REF))
        self.assertTrue(_triad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HET, ARRANGEMENT_HET))
        self.assertTrue(_triad_compatible(ARRANGEMENT_HOM_INV, ARRANGEMENT_HET, ARRANGEMENT_HOM_INV))


# ----------------------------------------------------------------------
# Transmission inference.
# ----------------------------------------------------------------------


class TestDyadTransmission(unittest.TestCase):
    def test_hom_ref_transmits_REF(self):
        for o in (ARRANGEMENT_HOM_REF, ARRANGEMENT_HET):
            self.assertEqual(_infer_dyad_transmission(ARRANGEMENT_HOM_REF, o), "REF")

    def test_hom_inv_transmits_INV(self):
        for o in (ARRANGEMENT_HOM_INV, ARRANGEMENT_HET):
            self.assertEqual(_infer_dyad_transmission(ARRANGEMENT_HOM_INV, o), "INV")

    def test_het_to_homozygous_offspring_resolvable(self):
        self.assertEqual(_infer_dyad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HOM_REF), "REF")
        self.assertEqual(_infer_dyad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HOM_INV), "INV")

    def test_het_to_het_ambiguous(self):
        self.assertEqual(_infer_dyad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HET), "ambiguous")


class TestTriadTransmission(unittest.TestCase):
    def test_hom_x_hom_opposite_het_offspring(self):
        # HOM_REF × HOM_INV → HET. Father → REF, Mother → INV.
        self.assertEqual(
            _infer_triad_transmission(ARRANGEMENT_HOM_REF, ARRANGEMENT_HOM_INV, ARRANGEMENT_HET),
            ("REF", "INV"),
        )

    def test_het_x_hom_inv_het_off(self):
        # HET × HOM_INV → HET. Father (het) transmitted REF; mother transmitted INV.
        self.assertEqual(
            _infer_triad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HOM_INV, ARRANGEMENT_HET),
            ("REF", "INV"),
        )

    def test_het_x_hom_inv_hom_inv_off(self):
        # Father transmitted INV.
        self.assertEqual(
            _infer_triad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HOM_INV, ARRANGEMENT_HOM_INV),
            ("INV", "INV"),
        )

    def test_het_x_het_het_off_ambiguous(self):
        self.assertEqual(
            _infer_triad_transmission(ARRANGEMENT_HET, ARRANGEMENT_HET, ARRANGEMENT_HET),
            ("ambiguous", "ambiguous"),
        )

    def test_contradiction_propagated(self):
        # HOM_INV × HET → HOM_REF is impossible.
        self.assertEqual(
            _infer_triad_transmission(ARRANGEMENT_HOM_INV, ARRANGEMENT_HET, ARRANGEMENT_HOM_REF),
            ("contradiction", "contradiction"),
        )


# ----------------------------------------------------------------------
# Binomial test.
# ----------------------------------------------------------------------


class TestBinomial(unittest.TestCase):
    def test_5050_passes(self):
        self.assertAlmostEqual(binomial_two_sided_pvalue(5, 10, 0.5), 1.0, places=6)

    def test_extreme_low_p(self):
        # 0 of 10 → p = 2 * (0.5^10) = ~0.001953
        self.assertAlmostEqual(binomial_two_sided_pvalue(0, 10, 0.5), 2 * (0.5 ** 10), places=6)

    def test_extreme_high_p(self):
        # 10 of 10 → same by symmetry
        self.assertAlmostEqual(binomial_two_sided_pvalue(10, 10, 0.5), 2 * (0.5 ** 10), places=6)

    def test_n_zero_returns_one(self):
        self.assertEqual(binomial_two_sided_pvalue(0, 0, 0.5), 1.0)

    def test_bad_k_raises(self):
        with self.assertRaises(ValueError):
            binomial_two_sided_pvalue(11, 10, 0.5)


# ----------------------------------------------------------------------
# End-to-end on synthetic_inversion fixture.
# ----------------------------------------------------------------------


class TestFixture(unittest.TestCase):
    def setUp(self):
        (
            self.inversion_id, self.polarity_hint,
            self.calls, self.dyads, self.triads,
        ) = load_karyotype_calls(FIX / "karyotype_calls.json")

    def test_load_counts(self):
        self.assertEqual(self.inversion_id, "inv_LG01_pod")
        self.assertEqual(self.polarity_hint, "band_0_is_REF")
        self.assertEqual(len(self.calls), 8)
        self.assertEqual(len(self.dyads), 8)
        self.assertEqual(len(self.triads), 4)

    def test_polarization_contradictions(self):
        res = polarize(
            inversion_id=self.inversion_id,
            karyotype_calls=self.calls,
            dyads=self.dyads,
            triads=self.triads,
            polarity_hint=self.polarity_hint,
        )
        # 1 dyad contradiction (S005→S008), 1 triad contradiction (S008 trio).
        self.assertEqual(res.dyad_compat["band_0_is_REF"].n_incompatible, 1)
        self.assertEqual(res.triad_compat["band_0_is_REF"].n_incompatible, 1)
        # The two orientations agree on contradiction counts (pure-data symmetry).
        self.assertTrue(res.polarities_symmetric)
        self.assertEqual(res.chosen_polarity, "band_0_is_REF")
        self.assertEqual(res.band_counts, {0: 2, 1: 3, 2: 3})

    def test_incompatible_dyad_is_S005_S008(self):
        res = polarize(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            polarity_hint=self.polarity_hint,
        )
        self.assertIn(("S005", "S008"), res.incompatible_dyads)
        # S004→S008 alone is compatible (HET parent could transmit REF).
        self.assertNotIn(("S004", "S008"), res.incompatible_dyads)

    def test_incompatible_triad_is_S008(self):
        res = polarize(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            polarity_hint=self.polarity_hint,
        )
        self.assertEqual(res.incompatible_triads, [("S004", "S005", "S008")])

    def test_transmissions_triad_overrides_dyad(self):
        rows = call_transmissions(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            chosen_polarity=self.polarity_hint,
        )
        # All 4 triads × 2 parents = 8 transmission rows. The dyads that
        # also appear in triads are suppressed.
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(r.relationship_type == "triad" for r in rows))

    def test_informative_drive_rows_are_two(self):
        rows = call_transmissions(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            chosen_polarity=self.polarity_hint,
        )
        informative = [r for r in rows if r.informative_for_drive]
        # S004 (HET) in triads 2 and 3 produces resolvable transmissions:
        #   triad 2 (offspring HET, co-parent HOM_INV) → S004 transmitted REF
        #   triad 3 (offspring HOM_INV, co-parent HOM_INV) → S004 transmitted INV
        # S004 in triad 4 produces a contradiction → NOT informative.
        self.assertEqual(len(informative), 2)
        parents = sorted(r.parent_sample_id for r in informative)
        self.assertEqual(parents, ["S004", "S004"])
        # One REF, one INV
        kinds = sorted(r.transmitted_arrangement for r in informative)
        self.assertEqual(kinds, ["INV", "REF"])

    def test_contradiction_rows_are_two(self):
        rows = call_transmissions(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            chosen_polarity=self.polarity_hint,
        )
        contras = [r for r in rows if r.transmitted_arrangement == "contradiction"]
        # Triad 4: both parents' transmissions tagged contradiction.
        self.assertEqual(len(contras), 2)
        self.assertEqual(sorted(r.offspring_sample_id for r in contras),
                         ["S008", "S008"])

    def test_drive_stats(self):
        rows = call_transmissions(
            inversion_id=self.inversion_id, karyotype_calls=self.calls,
            dyads=self.dyads, triads=self.triads,
            chosen_polarity=self.polarity_hint,
        )
        stats = drive_test(rows)
        self.assertEqual(stats.n_informative_transmissions, 2)
        self.assertEqual(stats.n_REF_transmitted, 1)
        self.assertEqual(stats.n_INV_transmitted, 1)
        self.assertAlmostEqual(stats.binomial_pvalue, 1.0, places=6)
        self.assertEqual(stats.paternal_n, 2)
        self.assertEqual(stats.maternal_n, 0)


# ----------------------------------------------------------------------
# IN adapter error paths.
# ----------------------------------------------------------------------


class TestInAdapterErrors(unittest.TestCase):
    def _write(self, body: dict) -> Path:
        fh = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(body, fh)
        fh.close()
        return Path(fh.name)

    def test_bad_schema_raises(self):
        p = self._write({"schema": "wrong_v0", "inversion_id": "i",
                         "polarity_hint": "band_0_is_REF",
                         "karyotype_calls": []})
        try:
            with self.assertRaises(NgsTractsAdapterError):
                load_karyotype_calls(p)
        finally:
            p.unlink()

    def test_bad_polarity_hint_raises(self):
        p = self._write({"schema": "ngspedigree_karyotype_calls_in_v1",
                         "inversion_id": "i",
                         "polarity_hint": "band_0_is_other",
                         "karyotype_calls": []})
        try:
            with self.assertRaises(NgsTractsAdapterError):
                load_karyotype_calls(p)
        finally:
            p.unlink()

    def test_bad_band_raises(self):
        p = self._write({"schema": "ngspedigree_karyotype_calls_in_v1",
                         "inversion_id": "i",
                         "polarity_hint": "band_0_is_REF",
                         "karyotype_calls": [{"sample_id": "S", "band": 3}]})
        try:
            with self.assertRaises(NgsTractsAdapterError):
                load_karyotype_calls(p)
        finally:
            p.unlink()


# ----------------------------------------------------------------------
# OUT adapter round-trip.
# ----------------------------------------------------------------------


class TestOutAdapterRoundTrip(unittest.TestCase):
    def test_round_trip(self):
        (
            inversion_id, polarity_hint, calls, dyads, triads,
        ) = load_karyotype_calls(FIX / "karyotype_calls.json")
        res = polarize(
            inversion_id=inversion_id, karyotype_calls=calls,
            dyads=dyads, triads=triads, polarity_hint=polarity_hint,
        )
        rows = call_transmissions(
            inversion_id=inversion_id, karyotype_calls=calls,
            dyads=dyads, triads=triads,
            chosen_polarity=res.chosen_polarity,
        )
        stats = drive_test(rows)
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            out = Path(fh.name)
        try:
            write_polarized_transmissions(
                out, result=res, transmissions=rows, drive_stats=stats,
            )
            doc = json.loads(out.read_text())
            self.assertEqual(doc["schema"],
                             "ngspedigree_polarized_transmissions_v1")
            self.assertEqual(doc["intended_consumer"], "ngsTracts")
            self.assertEqual(doc["inversion_id"], "inv_LG01_pod")
            self.assertEqual(doc["polarization"]["chosen_polarity"],
                             "band_0_is_REF")
            self.assertEqual(len(doc["transmissions"]), 8)
            self.assertEqual(doc["drive_stats"]["n_informative_transmissions"], 2)
        finally:
            out.unlink()


# ----------------------------------------------------------------------
# Schema-column alignment with the OUT JSON Schema.
# ----------------------------------------------------------------------


class TestOutSchemaSurface(unittest.TestCase):
    def test_out_schema_lists_required_top_level(self):
        s = json.loads(
            (SCHEMAS / "polarized_transmissions.out.schema.json").read_text()
        )
        self.assertEqual(
            set(s["required"]),
            {"schema", "intended_consumer", "inversion_id",
             "polarization", "transmissions", "drive_stats"},
        )

    def test_in_schema_lists_required_top_level(self):
        s = json.loads(
            (SCHEMAS / "karyotype_calls.in.schema.json").read_text()
        )
        self.assertEqual(
            set(s["required"]),
            {"schema", "inversion_id", "polarity_hint", "karyotype_calls"},
        )
