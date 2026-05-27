"""
test_adapters — adapter surface assertions.

Verifies:
  - stage3_real raises Stage3RealNotReadyError (swap-target stub)
  - variant_master placeholder classifies T1 / T2 / T3 correctly
  - kbc_adapter placeholder round-trips records
  - io.default_*_adapter wiring returns the expected concrete classes
  - io.write_tsv emits the columns in declared order
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp import io  # noqa: E402
from hpp.kbc_adapter import (  # noqa: E402
    KbcAssignment,
    KbcNotReadyError,
    PlaceholderKbcAdapter,
    load_kbc_table_b,
)
from hpp.stage3_real import Stage3RealNotReadyError  # noqa: E402
from hpp.stage3_real import (  # noqa: E402
    load as real_load,
)
from hpp.stage3_real import (  # noqa: E402
    load_dyad_map as real_load_dyad,
)
from hpp.stage3_real import (  # noqa: E402
    load_parent_phase as real_load_phase,
)
from hpp.stage3_real import (  # noqa: E402
    load_triad_map as real_load_triad,
)
from hpp.variant_master import (  # noqa: E402
    PlaceholderVariantMaster,
    VariantAnnotation,
    VariantMasterNotReadyError,
    load_variant_master,
)


class TestStage3RealStub(unittest.TestCase):
    def test_load_raises(self):
        with self.assertRaises(Stage3RealNotReadyError):
            real_load(dyad_map_path="/tmp/nope.tsv")

    def test_load_dyad_raises(self):
        with self.assertRaises(Stage3RealNotReadyError):
            real_load_dyad("/tmp/nope.tsv")

    def test_load_triad_raises(self):
        with self.assertRaises(Stage3RealNotReadyError):
            real_load_triad("/tmp/nope.tsv")

    def test_load_phase_raises(self):
        with self.assertRaises(Stage3RealNotReadyError):
            real_load_phase("/tmp/nope.tsv")

    def test_is_subclass_of_notimplemented(self):
        # downstream callers should be able to except NotImplementedError
        # to handle the swap-not-yet-ready case generically.
        self.assertTrue(issubclass(Stage3RealNotReadyError, NotImplementedError))


class TestVariantMaster(unittest.TestCase):
    def setUp(self):
        self.vm = PlaceholderVariantMaster()
        self.vm.add(VariantAnnotation(
            variant_id="LG01:100:A:T",
            gene_id="G1", transcript_id="T1",
            consequence="stop_gained", impact="HIGH",
            sift_class=None, vesm_llr=None,
            splice_subclass=None,
        ))
        self.vm.add(VariantAnnotation(
            variant_id="LG01:200:G:C",
            gene_id="G2", transcript_id="T2",
            consequence="missense_variant", impact="MODERATE",
            sift_class="deleterious", vesm_llr=-8.5,
            splice_subclass=None,
        ))
        self.vm.add(VariantAnnotation(
            variant_id="LG01:300:T:A",
            gene_id="G3", transcript_id="T3",
            consequence="synonymous_variant", impact="LOW",
            sift_class="tolerated", vesm_llr=-1.2,
            splice_subclass=None,
        ))
        self.vm.add(VariantAnnotation(
            variant_id="LG01:400:G:T",
            gene_id="G4", transcript_id="T4",
            consequence="splice_region_variant", impact="MODERATE",
            sift_class=None, vesm_llr=None,
            splice_subclass="splice_branch_classA",
        ))

    def test_tier1_LoF_only(self):
        self.assertTrue(self.vm.is_damaging("LG01:100:A:T", tier="T1"))   # stop_gained
        self.assertFalse(self.vm.is_damaging("LG01:200:G:C", tier="T1"))  # missense
        self.assertFalse(self.vm.is_damaging("LG01:300:T:A", tier="T1"))  # synonymous
        self.assertFalse(self.vm.is_damaging("LG01:400:G:T", tier="T1"))  # T2-only splice

    def test_tier2_adds_splice_classA(self):
        self.assertTrue(self.vm.is_damaging("LG01:100:A:T", tier="T2"))   # still LoF
        self.assertFalse(self.vm.is_damaging("LG01:200:G:C", tier="T2"))  # still no missense
        self.assertTrue(self.vm.is_damaging("LG01:400:G:T", tier="T2"))   # splice Class A

    def test_tier3_adds_missense(self):
        self.assertTrue(self.vm.is_damaging("LG01:100:A:T", tier="T3"))
        self.assertTrue(self.vm.is_damaging("LG01:200:G:C", tier="T3"))   # SIFT del + VESM<-7
        self.assertFalse(self.vm.is_damaging("LG01:300:T:A", tier="T3"))  # synonymous never
        self.assertTrue(self.vm.is_damaging("LG01:400:G:T", tier="T3"))

    def test_unknown_tier_raises(self):
        with self.assertRaises(ValueError):
            self.vm.is_damaging("LG01:100:A:T", tier="T4")

    def test_unknown_variant_not_damaging(self):
        self.assertFalse(self.vm.is_damaging("LG01:999:A:T", tier="T3"))

    def test_loader_stub_raises(self):
        with self.assertRaises(VariantMasterNotReadyError):
            load_variant_master("/tmp/nope.tsv")


class TestKbcAdapter(unittest.TestCase):
    def test_roundtrip(self):
        adapter = PlaceholderKbcAdapter()
        adapter.add(KbcAssignment(
            variant_id="LG01:100:A:T",
            inversion_id="inv_LG01_A",
            pod_segment="L",
            arrangement_background="A_private",
            assignment_confidence="high",
        ))
        rec = adapter.lookup("LG01:100:A:T")
        self.assertIsNotNone(rec)
        self.assertEqual(rec.arrangement_background, "A_private")
        self.assertIsNone(adapter.lookup("LG01:999:A:T"))

    def test_loader_stub_raises(self):
        with self.assertRaises(KbcNotReadyError):
            load_kbc_table_b("/tmp/nope.tsv")


class TestIoModule(unittest.TestCase):
    def test_default_stage3_returns_placeholder(self):
        adapter = io.default_stage3_adapter()
        # the module itself satisfies the Stage3Adapter protocol
        from hpp import stage3_placeholder
        self.assertIs(adapter, stage3_placeholder)

    def test_default_variant_master(self):
        adapter = io.default_variant_master_adapter()
        self.assertIsInstance(adapter, PlaceholderVariantMaster)

    def test_default_kbc(self):
        adapter = io.default_kbc_adapter()
        self.assertIsInstance(adapter, PlaceholderKbcAdapter)

    def test_write_tsv_preserves_column_order(self):
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            path = Path(fh.name)
        try:
            io.write_tsv(
                path,
                columns=["a", "b", "c"],
                rows=[{"a": "1", "b": "2", "c": "3"},
                      {"a": "4", "b": None, "c": "6"}],
            )
            lines = path.read_text().splitlines()
            self.assertEqual(lines[0], "a\tb\tc")
            self.assertEqual(lines[1], "1\t2\t3")
            self.assertEqual(lines[2], "4\t\t6")  # None → empty
        finally:
            path.unlink()


if __name__ == "__main__":
    unittest.main()
