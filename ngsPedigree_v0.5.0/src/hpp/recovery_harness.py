"""
Recovery-test harness for ngsPedigree synthetic panels.

Runs each Panel through the full pipeline (Stage 1 shadow classifier →
roster construction → pedigree extract → polarization → mtDNA check →
transmission calling) and reports recovery metrics versus the
ground-truth pedigree.

The Stage 1 step here is the stdlib shadow ``classify_edge_stdlib`` to
keep the test suite pandas-free; the pure-function logic is identical
to STEP_PED_01_annotate_relationships.classify_edge.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .inversion_polarization import (
    DyadPair,
    KaryotypeCall,
    TriadTrio,
    call_transmissions,
    drive_test,
    polarize,
)
from .mtdna_check import (
    build_not_supplied_block,
    build_validation_block,
    check_pedigree,
    filter_by_mtdna,
)
from .pedigree_extract import (
    PARENT_ROLES,
    OFFSPRING_ROLES,
)
from .relatedness_sim import (
    DEFAULT_THRESHOLDS,
    Pedigree,
    SimConfig,
    classify_edge_stdlib,
    derive_hubs_and_roster,
    expected_verdict,
    simulate_karyotype,
    simulate_mtdna,
    simulate_pairwise_table,
    true_relationship,
)
from .synthetic_panels import PanelConfig


# ----------------------------------------------------------------------
# Metrics dataclasses.
# ----------------------------------------------------------------------


@dataclass
class EdgeRecoveryReport:
    panel_name: str
    n_pairs: int
    n_correct: int
    n_misclassified: int
    confusion: Dict[Tuple[str, str], int] = field(default_factory=dict)
    accuracy: float = 0.0
    # by-class breakdown
    by_class_total: Dict[str, int] = field(default_factory=dict)
    by_class_correct: Dict[str, int] = field(default_factory=dict)
    by_class_recall: Dict[str, float] = field(default_factory=dict)


@dataclass
class PipelineReport:
    panel_name: str
    n_individuals: int
    n_true_dyads: int          # parent->offspring pairs in ground truth
    n_true_triads: int          # triads in ground truth
    n_recovered_dyads: int      # dyads emerging from extract step
    n_recovered_triads: int
    edge_recovery: EdgeRecoveryReport = field(default_factory=lambda:
        EdgeRecoveryReport("", 0, 0, 0))
    # mtDNA injection / detection
    n_mtdna_injected_swaps: int = 0
    n_mtdna_incompatible: int = 0
    n_mtdna_swaps_detected: int = 0
    # Polarization
    n_polarization_dyad_incompat: int = 0
    n_polarization_triad_incompat: int = 0
    n_transmissions_emitted: int = 0
    n_informative: int = 0


# ----------------------------------------------------------------------
# Edge-classification recovery.
# ----------------------------------------------------------------------


def edge_recovery(panel_name: str, edges: List[Dict]) -> EdgeRecoveryReport:
    rep = EdgeRecoveryReport(panel_name, len(edges), 0, 0)
    for e in edges:
        truth = expected_verdict(e["true_relationship"])
        verdict = classify_edge_stdlib(e, DEFAULT_THRESHOLDS)
        predicted = verdict["edge_class"]
        e["edge_class"] = predicted
        e["confidence"] = verdict["confidence"]
        e["reasons"] = verdict["reasons"]
        if predicted == truth:
            rep.n_correct += 1
            rep.by_class_correct[truth] = rep.by_class_correct.get(truth, 0) + 1
        else:
            rep.n_misclassified += 1
        rep.by_class_total[truth] = rep.by_class_total.get(truth, 0) + 1
        key = (truth, predicted)
        rep.confusion[key] = rep.confusion.get(key, 0) + 1
    rep.accuracy = rep.n_correct / rep.n_pairs if rep.n_pairs else 0.0
    for cls, n in rep.by_class_total.items():
        rep.by_class_recall[cls] = rep.by_class_correct.get(cls, 0) / n if n else 0.0
    return rep


# ----------------------------------------------------------------------
# Pedigree → Stage1 roster (truth-based) → extract dyads/triads via
# the same in-pipeline extract module.
# ----------------------------------------------------------------------


def _roster_to_dict(roster: List[Dict]) -> Dict[str, Dict]:
    return {r["sample_id"]: r for r in roster}


def _extract_dyads_from_roster(
    pedigree: Pedigree, edges: List[Dict], roster: List[Dict],
) -> List[DyadPair]:
    """Mirror pedigree_extract.extract_dyads but on in-memory rows."""
    by_role = _roster_to_dict(roster)
    out: List[DyadPair] = []
    for row in edges:
        if row.get("edge_class") != "parent_offspring":
            continue
        a, b = row["sample_a"], row["sample_b"]
        ra = by_role.get(a)
        rb = by_role.get(b)
        if not ra or not rb:
            continue
        ra_is_p = ra["possible_role"] in PARENT_ROLES
        rb_is_p = rb["possible_role"] in PARENT_ROLES
        ra_is_o = ra["possible_role"] in OFFSPRING_ROLES
        rb_is_o = rb["possible_role"] in OFFSPRING_ROLES
        if ra_is_p and rb_is_o:
            sex = ra["possible_role"] if ra["possible_role"] in ("mother", "father") else None
            sex = {"father": "male", "mother": "female"}.get(sex)
            out.append(DyadPair(parent_sample_id=a, offspring_sample_id=b,
                                 parent_sex=sex))
        elif rb_is_p and ra_is_o:
            sex = rb["possible_role"] if rb["possible_role"] in ("mother", "father") else None
            sex = {"father": "male", "mother": "female"}.get(sex)
            out.append(DyadPair(parent_sample_id=b, offspring_sample_id=a,
                                 parent_sex=sex))
    return out


def _extract_triads_from_roster(
    pedigree: Pedigree, edges: List[Dict], roster: List[Dict],
) -> List[TriadTrio]:
    by_hub: Dict[str, List[Dict]] = {}
    for r in roster:
        by_hub.setdefault(r["hub_id"], []).append(r)
    # PO neighbour set
    po: Dict[str, set] = {}
    for row in edges:
        if row.get("edge_class") != "parent_offspring":
            continue
        a, b = row["sample_a"], row["sample_b"]
        po.setdefault(a, set()).add(b)
        po.setdefault(b, set()).add(a)
    out: List[TriadTrio] = []
    for hub_id, members in by_hub.items():
        if not members or members[0]["hub_type"] != "two_parents_with_sibship":
            continue
        father = [m for m in members if m["possible_role"] == "father"]
        mother = [m for m in members if m["possible_role"] == "mother"]
        offspring = [m for m in members if m["possible_role"] == "possible_offspring"]
        if len(father) != 1 or len(mother) != 1:
            continue
        pat = father[0]["sample_id"]
        mat = mother[0]["sample_id"]
        for o in offspring:
            sid = o["sample_id"]
            if pat in po.get(sid, set()) and mat in po.get(sid, set()):
                out.append(TriadTrio(pat, mat, sid))
    return out


# ----------------------------------------------------------------------
# Ground-truth counts (for assertion baselines).
# ----------------------------------------------------------------------


def ground_truth_counts(pedigree: Pedigree) -> Tuple[int, int]:
    n_dyads = 0
    n_triads = 0
    for ind in pedigree.individuals:
        if ind.father:
            n_dyads += 1
        if ind.mother:
            n_dyads += 1
        if ind.father and ind.mother:
            n_triads += 1
    return n_dyads, n_triads


# ----------------------------------------------------------------------
# Full-pipeline runner.
# ----------------------------------------------------------------------


def run_panel(
    pedigree: Pedigree,
    cfg: PanelConfig,
    *,
    sim_config: Optional[SimConfig] = None,
) -> PipelineReport:
    """Run one panel through the recovery pipeline. Returns the report."""
    rng = random.Random(cfg.seed)
    sim_cfg = sim_config or SimConfig()

    # 1. Simulate pairwise coefficients.
    edges = simulate_pairwise_table(pedigree, rng, sim_cfg)

    # 2. Edge classification recovery (stdlib shadow).
    edge_rep = edge_recovery(cfg.name, edges)

    # 3. Build truth-roster (per-family hubs/roles).
    edges, roster = derive_hubs_and_roster(pedigree, edges)

    # 4. Extract dyads/triads through the pipeline logic.
    dyads = _extract_dyads_from_roster(pedigree, edges, roster)
    triads = _extract_triads_from_roster(pedigree, edges, roster)

    # 5. Simulate karyotype + mtDNA.
    bands = simulate_karyotype(pedigree, cfg.inversion_id, rng,
                                inv_allele_freq=cfg.inv_allele_freq)
    # Pick offspring to swap mtDNA on.
    offspring_with_mother = [
        ind.sample_id for ind in pedigree.individuals if ind.mother
    ]
    n_swap = int(round(cfg.mtdna_swap_rate * len(offspring_with_mother)))
    swap_ids = rng.sample(offspring_with_mother, n_swap) if n_swap else []
    mtdna_full = simulate_mtdna(pedigree, rng, swap_offspring_ids=swap_ids)

    # mtDNA coverage: drop a fraction of samples from the mtDNA table.
    all_samples = list(mtdna_full.keys())
    n_keep = int(round(cfg.mtdna_coverage * len(all_samples)))
    kept = set(rng.sample(all_samples, n_keep))
    mtdna_records = {s: h for s, h in mtdna_full.items() if s in kept}

    # mtDNA check via the in-repo machinery.
    if mtdna_records:
        from .mtdna_check import MtdnaRecord
        mt_recs = {sid: MtdnaRecord(sample_id=sid, haplotype=hap)
                   for sid, hap in mtdna_records.items()}
        mt_checks = check_pedigree(mt_recs, dyads, triads)
        # swap detection: count how many incompatible checks correspond to
        # injected swap offspring whose mother is in the kept mtDNA set.
        swap_set = set(swap_ids)
        n_swap_detected = sum(
            1 for c in mt_checks
            if c.status == "incompatible" and c.offspring_sample_id in swap_set
        )
        n_incompat = sum(1 for c in mt_checks if c.status == "incompatible")
        dyads_filt, triads_filt = filter_by_mtdna(dyads, triads, mt_checks)
    else:
        n_swap_detected = 0
        n_incompat = 0
        dyads_filt, triads_filt = list(dyads), list(triads)

    # 6. Polarization.
    calls = [KaryotypeCall(sid, cfg.inversion_id, b) for sid, b in bands.items()]
    res = polarize(
        inversion_id=cfg.inversion_id, karyotype_calls=calls,
        dyads=dyads_filt, triads=triads_filt,
        polarity_hint="band_0_is_REF",
    )
    transmissions = call_transmissions(
        inversion_id=cfg.inversion_id, karyotype_calls=calls,
        dyads=dyads_filt, triads=triads_filt,
        chosen_polarity=res.chosen_polarity,
    )
    stats = drive_test(transmissions)

    n_true_dyads, n_true_triads = ground_truth_counts(pedigree)
    return PipelineReport(
        panel_name=cfg.name,
        n_individuals=len(pedigree.individuals),
        n_true_dyads=n_true_dyads,
        n_true_triads=n_true_triads,
        n_recovered_dyads=len(dyads),
        n_recovered_triads=len(triads),
        edge_recovery=edge_rep,
        n_mtdna_injected_swaps=len(swap_ids),
        n_mtdna_incompatible=n_incompat,
        n_mtdna_swaps_detected=n_swap_detected,
        n_polarization_dyad_incompat=res.dyad_compat[res.chosen_polarity].n_incompatible,
        n_polarization_triad_incompat=res.triad_compat[res.chosen_polarity].n_incompatible,
        n_transmissions_emitted=len(transmissions),
        n_informative=stats.n_informative_transmissions,
    )


def format_report(rep: PipelineReport) -> str:
    er = rep.edge_recovery
    lines = [
        f"=== {rep.panel_name} ===",
        f"  individuals          : {rep.n_individuals}",
        f"  ground-truth dyads   : {rep.n_true_dyads}",
        f"  ground-truth triads  : {rep.n_true_triads}",
        f"  recovered dyads      : {rep.n_recovered_dyads}",
        f"  recovered triads     : {rep.n_recovered_triads}",
        f"  edge accuracy        : {er.accuracy * 100:.2f}%  "
        f"({er.n_correct}/{er.n_pairs})",
        "  by-class recall:",
    ]
    for cls in sorted(er.by_class_total):
        n = er.by_class_total[cls]
        if n == 0:
            continue
        lines.append(
            f"    {cls:24s}  "
            f"recall={er.by_class_recall[cls] * 100:6.2f}%  ({er.by_class_correct.get(cls, 0)}/{n})"
        )
    lines += [
        f"  mtDNA: injected_swaps={rep.n_mtdna_injected_swaps}  "
        f"incompatible={rep.n_mtdna_incompatible}  "
        f"swaps_detected={rep.n_mtdna_swaps_detected}",
        f"  polarization: dyad_incompat={rep.n_polarization_dyad_incompat}  "
        f"triad_incompat={rep.n_polarization_triad_incompat}",
        f"  transmissions: total={rep.n_transmissions_emitted}  "
        f"informative={rep.n_informative}",
    ]
    return "\n".join(lines)
