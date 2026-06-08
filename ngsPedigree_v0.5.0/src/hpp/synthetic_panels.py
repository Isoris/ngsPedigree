"""
Synthetic test panels for ngsPedigree pipeline validation.

Each ``Panel`` is a deterministic (seeded) synthetic cohort with mixed
family topology, varied mtDNA coverage, and known ground truth. Panels
are designed to cover the failure modes the real catfish cohort will
have to handle.

Panel inventory (matches the user-requested coverage):

  Panel A — many_small_mixed
      50 families, 1-5 catfish each. Mix of dyads (40%), single-
      offspring triads (30%), small sibships (30%). ~50% have mtDNA.

  Panel B — few_large_triads
      10 families, 8-12 catfish each. All triads with multi-offspring
      sibships. Full mtDNA coverage.

  Panel C — dyad_only
      30 PO dyads, no triads, no sibships. ~70% have mtDNA.

  Panel D — triad_only
      20 single-offspring triads. Full mtDNA.

  Panel E — many_small_dyads
      60 PO dyads. ~80% have mtDNA, some swapped to inject incompatibility.

  Panel F — eighty_twenty
      70 small families (2-5 catfish), 80% dyads / 20% triads. Mixed mtDNA.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from typing import List, Optional

from .relatedness_sim import (
    Individual,
    Pedigree,
    make_dyad_family,
    make_single_parent_sibship,
    make_triad_family,
)


@dataclass(frozen=True)
class PanelConfig:
    name: str
    seed: int
    inversion_id: str = "inv_LG01_pod"
    mtdna_coverage: float = 1.0       # fraction of samples with mtDNA call
    mtdna_swap_rate: float = 0.0      # fraction of offspring with deliberate maternal mtDNA mismatch
    inv_allele_freq: float = 0.4


def _rng(seed: int) -> random.Random:
    return random.Random(seed)


def build_panel_A(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel A — many_small_mixed: 50 families."""
    cfg = cfg or PanelConfig("panel_A_many_small_mixed", 1001)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(50):
        fid = f"famA{i+1:03d}"
        prefix = f"A{i+1:03d}"
        r = rng.random()
        if r < 0.40:
            sexes = rng.choice([("male", "female"), ("female", "male"),
                                 ("male", "male"), ("female", "female")])
            individuals.extend(make_dyad_family(fid, prefix, sexes))
        elif r < 0.70:
            individuals.extend(make_triad_family(fid, prefix, n_offspring=1, rng=rng))
        else:
            n = rng.randint(2, 3)
            individuals.extend(make_triad_family(fid, prefix, n_offspring=n, rng=rng))
    return Pedigree(name=cfg.name, individuals=individuals)


def build_panel_B(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel B — few_large_triads: 10 families × 8-12 each."""
    cfg = cfg or PanelConfig("panel_B_few_large_triads", 1002)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(10):
        fid = f"famB{i+1:02d}"
        prefix = f"B{i+1:02d}"
        n_off = rng.randint(6, 10)   # 6-10 offspring + 2 parents = 8-12 total
        individuals.extend(make_triad_family(fid, prefix, n_offspring=n_off, rng=rng))
    return Pedigree(name=cfg.name, individuals=individuals)


def build_panel_C(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel C — dyad_only: 30 isolated PO dyads."""
    cfg = cfg or PanelConfig("panel_C_dyad_only", 1003)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(30):
        fid = f"famC{i+1:03d}"
        prefix = f"C{i+1:03d}"
        sexes = rng.choice([("male", "female"), ("female", "male"),
                             ("male", "male"), ("female", "female")])
        individuals.extend(make_dyad_family(fid, prefix, sexes))
    return Pedigree(name=cfg.name, individuals=individuals)


def build_panel_D(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel D — triad_only: 20 single-offspring triads."""
    cfg = cfg or PanelConfig("panel_D_triad_only", 1004)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(20):
        fid = f"famD{i+1:03d}"
        prefix = f"D{i+1:03d}"
        individuals.extend(make_triad_family(fid, prefix, n_offspring=1, rng=rng))
    return Pedigree(name=cfg.name, individuals=individuals)


def build_panel_E(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel E — many_small_dyads: 60 PO dyads."""
    cfg = cfg or PanelConfig("panel_E_many_small_dyads", 1005,
                              mtdna_coverage=0.8, mtdna_swap_rate=0.08)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(60):
        fid = f"famE{i+1:03d}"
        prefix = f"E{i+1:03d}"
        sexes = rng.choice([("male", "female"), ("female", "male"),
                             ("male", "male"), ("female", "female")])
        individuals.extend(make_dyad_family(fid, prefix, sexes))
    return Pedigree(name=cfg.name, individuals=individuals)


def build_panel_F(cfg: Optional[PanelConfig] = None) -> Pedigree:
    """Panel F — eighty_twenty: 70 families, 80% dyads / 20% triads."""
    cfg = cfg or PanelConfig("panel_F_eighty_twenty", 1006)
    rng = _rng(cfg.seed)
    individuals: List[Individual] = []
    for i in range(70):
        fid = f"famF{i+1:03d}"
        prefix = f"F{i+1:03d}"
        if rng.random() < 0.80:
            sexes = rng.choice([("male", "female"), ("female", "male"),
                                 ("male", "male"), ("female", "female")])
            individuals.extend(make_dyad_family(fid, prefix, sexes))
        else:
            n_off = rng.randint(1, 3)
            individuals.extend(make_triad_family(fid, prefix,
                                                  n_offspring=n_off, rng=rng))
    return Pedigree(name=cfg.name, individuals=individuals)


PANEL_BUILDERS = {
    "A_many_small_mixed":   (build_panel_A, PanelConfig("panel_A_many_small_mixed", 1001, mtdna_coverage=0.5)),
    "B_few_large_triads":   (build_panel_B, PanelConfig("panel_B_few_large_triads", 1002, mtdna_coverage=1.0)),
    "C_dyad_only":          (build_panel_C, PanelConfig("panel_C_dyad_only", 1003, mtdna_coverage=0.7)),
    "D_triad_only":         (build_panel_D, PanelConfig("panel_D_triad_only", 1004, mtdna_coverage=1.0)),
    "E_many_small_dyads":   (build_panel_E, PanelConfig("panel_E_many_small_dyads", 1005,
                                                          mtdna_coverage=0.8, mtdna_swap_rate=0.08)),
    "F_eighty_twenty":      (build_panel_F, PanelConfig("panel_F_eighty_twenty", 1006, mtdna_coverage=0.6)),
}
