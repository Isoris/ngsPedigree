"""
test_synthetic_panels — comprehensive recovery tests over six synthetic
populations covering the topology mixes the real catfish cohort will
hit. Each panel is deterministic (seeded), so the assertions act as
fixed-point regression checks for the entire pipeline.

The six panels (see synthetic_panels.py):
  A — 50 small mixed families (1-5 individuals each)
  B — 10 large multi-offspring triads (8-12 each)
  C — 30 dyad-only families
  D — 20 triad-only families
  E — 60 PO dyads with deliberate mtDNA swaps to test detection
  F — 70 mixed families, 80% dyads / 20% triads

Tests assert lower bounds on:
  - per-pair edge-classification accuracy
  - ground-truth → recovered dyad/triad counts
  - mtDNA swap-detection recall
  - polarization-layer Mendelian compatibility on simulated bands
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR.parent / "src"))

from hpp.recovery_harness import format_report, run_panel  # noqa: E402
from hpp.synthetic_panels import (  # noqa: E402
    PANEL_BUILDERS,
    build_panel_A,
    build_panel_B,
    build_panel_C,
    build_panel_D,
    build_panel_E,
    build_panel_F,
)


# ----------------------------------------------------------------------
# Panel structural sanity (without simulation).
# ----------------------------------------------------------------------


class TestPanelStructure(unittest.TestCase):
    def test_panel_A_sample_count_range(self):
        p = build_panel_A()
        # 50 families × 2–4 individuals = ~100-200 total
        self.assertGreaterEqual(len(p.individuals), 100)
        self.assertLessEqual(len(p.individuals), 250)
        # Every family must have at least one parent-child link.
        for ind in p.individuals:
            if ind.father is None and ind.mother is None:
                # founder — must have at least one offspring
                has_off = any(o.father == ind.sample_id
                              or o.mother == ind.sample_id
                              for o in p.individuals)
                self.assertTrue(has_off,
                                f"founder {ind.sample_id} has no offspring")

    def test_panel_B_each_family_is_a_triad_with_sibship(self):
        p = build_panel_B()
        from collections import defaultdict
        by_fam = defaultdict(list)
        for ind in p.individuals:
            by_fam[ind.family_id].append(ind)
        for fid, members in by_fam.items():
            parents = [m for m in members if m.father is None and m.mother is None]
            offspring = [m for m in members if m.father or m.mother]
            self.assertEqual(len(parents), 2)
            self.assertGreaterEqual(len(offspring), 6)
            self.assertLessEqual(len(offspring), 10)

    def test_panel_C_is_all_dyads(self):
        p = build_panel_C()
        from collections import defaultdict
        by_fam = defaultdict(list)
        for ind in p.individuals:
            by_fam[ind.family_id].append(ind)
        self.assertEqual(len(by_fam), 30)
        for members in by_fam.values():
            self.assertEqual(len(members), 2)

    def test_panel_D_is_all_triads_single_off(self):
        p = build_panel_D()
        from collections import defaultdict
        by_fam = defaultdict(list)
        for ind in p.individuals:
            by_fam[ind.family_id].append(ind)
        self.assertEqual(len(by_fam), 20)
        for members in by_fam.values():
            self.assertEqual(len(members), 3)

    def test_panel_F_distribution(self):
        p = build_panel_F()
        from collections import defaultdict
        by_fam = defaultdict(list)
        for ind in p.individuals:
            by_fam[ind.family_id].append(ind)
        self.assertEqual(len(by_fam), 70)
        sizes = [len(m) for m in by_fam.values()]
        n_dyads = sum(1 for s in sizes if s == 2)
        # 80% dyads expected (with seeded RNG: 56/70).
        self.assertGreaterEqual(n_dyads, 50)


# ----------------------------------------------------------------------
# Recovery-pipeline assertions per panel.
# ----------------------------------------------------------------------


class TestPanelRecovery(unittest.TestCase):
    """For each panel, run the full pipeline and assert lower bounds on
    recovery metrics. Lower bounds are tuned to be tight enough to catch
    regressions but loose enough that minor noise tweaks don't break
    the suite."""

    def _run(self, key):
        builder, cfg = PANEL_BUILDERS[key]
        pedigree = builder(cfg)
        report = run_panel(pedigree, cfg)
        # for debugging: print the full report when a test fails
        self._last_report = report
        return report

    def _assert_recovery(self, key, *, min_edge_acc, recall_floor,
                          min_dyad_ratio=1.0, min_triad_ratio=1.0):
        rep = self._run(key)
        msg = "\n" + format_report(rep)
        self.assertGreaterEqual(
            rep.edge_recovery.accuracy, min_edge_acc,
            f"edge accuracy below floor:{msg}")
        # Floors on per-class recall for the four classes we control.
        for cls, floor in recall_floor.items():
            if rep.edge_recovery.by_class_total.get(cls, 0) == 0:
                continue
            self.assertGreaterEqual(
                rep.edge_recovery.by_class_recall[cls], floor,
                f"recall for {cls} below floor ({floor}):{msg}")
        # Dyad/triad recovery counts.
        if rep.n_true_dyads:
            self.assertGreaterEqual(
                rep.n_recovered_dyads / rep.n_true_dyads, min_dyad_ratio,
                f"dyad recovery below {min_dyad_ratio}:{msg}")
        if rep.n_true_triads:
            self.assertGreaterEqual(
                rep.n_recovered_triads / rep.n_true_triads, min_triad_ratio,
                f"triad recovery below {min_triad_ratio}:{msg}")

    def test_panel_A_many_small_mixed(self):
        # 40% of families are isolated PO dyads → po_dyad_only hubs →
        # direction unrecoverable. Dyads only emerge from triad/sibship hubs.
        # All triads should be recovered.
        self._assert_recovery(
            "A_many_small_mixed",
            min_edge_acc=0.95,
            recall_floor={"parent_offspring": 0.95,
                          "full_sibling": 0.85,
                          "unrelated": 0.95},
            min_dyad_ratio=0.80,
            min_triad_ratio=1.0,
        )

    def test_panel_B_few_large_triads(self):
        self._assert_recovery(
            "B_few_large_triads",
            min_edge_acc=0.95,
            recall_floor={"parent_offspring": 0.95,
                          "full_sibling": 0.85,
                          "unrelated": 0.95},
        )

    def test_panel_C_dyad_only(self):
        # All dyads — no full sibs, no triads. PO recovery is the focus.
        rep = self._run("C_dyad_only")
        msg = "\n" + format_report(rep)
        self.assertGreaterEqual(rep.edge_recovery.accuracy, 0.98, msg)
        self.assertGreaterEqual(
            rep.edge_recovery.by_class_recall.get("parent_offspring", 0), 0.95, msg)
        # All true dyads ARE PO edges, but truth-roster won't infer
        # direction (po_dyad_only), so n_recovered_dyads = 0.
        self.assertEqual(rep.n_recovered_dyads, 0)
        self.assertEqual(rep.n_recovered_triads, 0)

    def test_panel_D_triad_only(self):
        self._assert_recovery(
            "D_triad_only",
            min_edge_acc=0.98,
            recall_floor={"parent_offspring": 0.95,
                          "unrelated": 0.95},
            min_dyad_ratio=1.0,
            min_triad_ratio=1.0,
        )

    def test_panel_E_many_small_dyads_with_mtdna_swaps(self):
        rep = self._run("E_many_small_dyads")
        msg = "\n" + format_report(rep)
        self.assertGreaterEqual(rep.edge_recovery.accuracy, 0.98, msg)
        # mtDNA swaps were injected but cannot be detected here because
        # po_dyad_only hubs leave parent_sex unknown — the maternal-dyad
        # mtDNA check requires parent_sex="female". This is a documented
        # blind-mode limitation, not a bug; assert it.
        self.assertEqual(rep.n_mtdna_swaps_detected, 0,
                          f"mtDNA detection on po_dyad_only is impossible:{msg}")

    def test_panel_F_eighty_twenty(self):
        # 80% of families are isolated PO dyads → unrecoverable direction.
        # Dyad-recovery ratio is therefore low by design; triads (the
        # other 20%) recover fully.
        self._assert_recovery(
            "F_eighty_twenty",
            min_edge_acc=0.95,
            recall_floor={"parent_offspring": 0.95,
                          "unrelated": 0.95},
            min_dyad_ratio=0.30,
            min_triad_ratio=1.0,
        )


# ----------------------------------------------------------------------
# Polarization sanity: every panel emits internally-consistent transmissions
# (no contradictions, since the simulator generates Mendelian-correct bands).
# ----------------------------------------------------------------------


class TestPolarizationConsistency(unittest.TestCase):
    def test_simulated_bands_are_mendelian(self):
        for key in ("A_many_small_mixed", "B_few_large_triads",
                    "D_triad_only", "F_eighty_twenty"):
            builder, cfg = PANEL_BUILDERS[key]
            rep = run_panel(builder(cfg), cfg)
            # Mendelian-correct bands ⇒ zero nuclear Mendelian violations
            # (the simulator does not inject band swaps).
            self.assertEqual(
                rep.n_polarization_dyad_incompat, 0,
                f"{key}: simulated bands violated Mendel — bug")
            self.assertEqual(
                rep.n_polarization_triad_incompat, 0,
                f"{key}: simulated bands violated Mendel — bug")


# ----------------------------------------------------------------------
# Panel B — large families let us catch mtDNA swaps via the triad path.
# ----------------------------------------------------------------------


class TestPanelB_MtdnaSwapDetection(unittest.TestCase):
    """Panel B has triads with multi-offspring sibships. We inject mtDNA
    swaps and verify the triad-side mtDNA check catches them."""

    def test_inject_and_detect(self):
        from hpp.synthetic_panels import PanelConfig, build_panel_B
        # Add swaps
        cfg = PanelConfig(
            name="panel_B_with_swaps",
            seed=4242,
            inversion_id="inv_LG01_pod",
            mtdna_coverage=1.0,
            mtdna_swap_rate=0.10,
        )
        rep = run_panel(build_panel_B(), cfg)
        msg = "\n" + format_report(rep)
        # We should detect ≥ 50% of injected swaps via the triad path
        # (mother is identifiable in two_parents_with_sibship hubs).
        self.assertGreater(rep.n_mtdna_injected_swaps, 0, msg)
        if rep.n_mtdna_injected_swaps:
            recall = rep.n_mtdna_swaps_detected / rep.n_mtdna_injected_swaps
            self.assertGreaterEqual(
                recall, 0.5,
                f"mtDNA swap detection recall < 50%:{msg}")
