"""
Synthetic-population simulator for ngsPedigree pipeline testing.

Generates ground-truth pedigrees with varied topology, simulates
ngsRelate-style (theta, IBS0) coefficients per pair, and synthesizes
PCAngsd karyotype bands + mtDNA haplotypes that respect Mendelian
inheritance. Outputs are Stage-1-compatible TSVs that feed straight
into the existing pipeline.

The simulator is stdlib-only. The shadow ``classify_edge_stdlib``
mirrors STEP_PED_01_annotate_relationships.classify_edge so the
recovery tests can run without pandas installed.

References used for expected (theta, IBS0):
  - Manichaikul et al. 2010 KING-robust thresholds (first/second/third
    degree at 0.177, 0.0884, 0.0442).
  - Standard pedigree-kinship values:
      duplicate / MZ twin:  theta = 0.50
      parent–offspring:      theta = 0.25,  IBS0 ≈ 0   (always share ≥1 allele)
      full sibling:          theta = 0.25,  IBS0 ≈ 0.0125
      half sibling:          theta = 0.125
      first cousin:          theta = 0.0625
      unrelated:             theta = 0.0
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


# ----------------------------------------------------------------------
# Pedigree representation.
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class Individual:
    sample_id: str
    sex: str            # "male" | "female"
    father: Optional[str]
    mother: Optional[str]
    family_id: str
    generation: int     # 0 = founder, 1 = F1, ...


@dataclass
class Pedigree:
    name: str
    individuals: List[Individual] = field(default_factory=list)

    def by_id(self) -> Dict[str, Individual]:
        return {i.sample_id: i for i in self.individuals}

    def offspring_of(self, sample_id: str) -> List[Individual]:
        return [i for i in self.individuals
                if i.father == sample_id or i.mother == sample_id]

    def founders(self) -> List[Individual]:
        return [i for i in self.individuals
                if i.father is None and i.mother is None]


# ----------------------------------------------------------------------
# Pedigree builders.
# ----------------------------------------------------------------------


def make_dyad_family(family_id: str, prefix: str, sexes: Tuple[str, str]) -> List[Individual]:
    """One parent + one offspring (single PO dyad)."""
    p_sex = sexes[0]
    o_sex = sexes[1]
    parent = Individual(
        sample_id=f"{prefix}_P", sex=p_sex,
        father=None, mother=None, family_id=family_id, generation=0,
    )
    if p_sex == "male":
        offspring = Individual(
            sample_id=f"{prefix}_O", sex=o_sex,
            father=parent.sample_id, mother=None,
            family_id=family_id, generation=1,
        )
    else:
        offspring = Individual(
            sample_id=f"{prefix}_O", sex=o_sex,
            father=None, mother=parent.sample_id,
            family_id=family_id, generation=1,
        )
    return [parent, offspring]


def make_triad_family(
    family_id: str, prefix: str, n_offspring: int = 1,
    rng: Optional[random.Random] = None,
) -> List[Individual]:
    """Two unrelated parents (one male, one female) + n_offspring."""
    rng = rng or random.Random()
    father = Individual(
        sample_id=f"{prefix}_F", sex="male",
        father=None, mother=None, family_id=family_id, generation=0,
    )
    mother = Individual(
        sample_id=f"{prefix}_M", sex="female",
        father=None, mother=None, family_id=family_id, generation=0,
    )
    offspring: List[Individual] = []
    for i in range(n_offspring):
        offspring.append(Individual(
            sample_id=f"{prefix}_O{i+1}",
            sex=rng.choice(["male", "female"]),
            father=father.sample_id,
            mother=mother.sample_id,
            family_id=family_id, generation=1,
        ))
    return [father, mother] + offspring


def make_single_parent_sibship(
    family_id: str, prefix: str, n_offspring: int,
    rng: Optional[random.Random] = None,
) -> List[Individual]:
    """One known parent + sibship of n_offspring (other parent absent)."""
    rng = rng or random.Random()
    p_sex = rng.choice(["male", "female"])
    parent = Individual(
        sample_id=f"{prefix}_P", sex=p_sex,
        father=None, mother=None, family_id=family_id, generation=0,
    )
    offspring: List[Individual] = []
    for i in range(n_offspring):
        os = rng.choice(["male", "female"])
        offspring.append(Individual(
            sample_id=f"{prefix}_O{i+1}", sex=os,
            father=parent.sample_id if p_sex == "male" else None,
            mother=parent.sample_id if p_sex == "female" else None,
            family_id=family_id, generation=1,
        ))
    return [parent] + offspring


# ----------------------------------------------------------------------
# Pedigree-relationship lookup.
# ----------------------------------------------------------------------

# Relationship classes used by the classifier
REL_DUPLICATE = "duplicate_or_clone"
REL_PO = "parent_offspring"
REL_FS = "full_sibling"
REL_HS = "half_sibling"          # second-degree
REL_AVUNCULAR = "avuncular"       # second-degree
REL_FIRST_COUSIN = "first_cousin" # third-degree
REL_UNRELATED = "unrelated"


def true_relationship(pedigree: Pedigree, a: str, b: str) -> str:
    """Compute the canonical pedigree relationship between two samples
    by walking the pedigree."""
    by = pedigree.by_id()
    ia, ib = by[a], by[b]
    if ia.sample_id == ib.sample_id:
        return REL_DUPLICATE
    # Parent-offspring (direct).
    if ib.father == ia.sample_id or ib.mother == ia.sample_id:
        return REL_PO
    if ia.father == ib.sample_id or ia.mother == ib.sample_id:
        return REL_PO
    # Full siblings share both parents.
    if (ia.father is not None and ia.father == ib.father
            and ia.mother is not None and ia.mother == ib.mother):
        return REL_FS
    # Half siblings share one parent.
    if (ia.father is not None and ia.father == ib.father) or \
       (ia.mother is not None and ia.mother == ib.mother):
        return REL_HS
    # No deeper relationships are constructed in these panels (no
    # grandparental or cousin links). Anyone else: unrelated.
    return REL_UNRELATED


# Expected (theta, IBS0) per relationship class.
_TRUE_COEFFS = {
    REL_DUPLICATE:     (0.50, 0.0000),
    REL_PO:            (0.25, 0.0002),
    REL_FS:            (0.25, 0.0125),
    REL_HS:            (0.125, 0.04),
    REL_AVUNCULAR:     (0.125, 0.04),
    REL_FIRST_COUSIN:  (0.0625, 0.08),
    REL_UNRELATED:     (0.0,  0.12),
}


def expected_coeffs(rel: str) -> Tuple[float, float]:
    return _TRUE_COEFFS.get(rel, _TRUE_COEFFS[REL_UNRELATED])


# ----------------------------------------------------------------------
# Coefficient simulator.
# ----------------------------------------------------------------------


@dataclass
class SimConfig:
    """Defaults match ngsRelate at ~100k+ sites of genome-wide data.

    sigma_ibs0 is intentionally small (0.0008): for parent-offspring
    pairs the true IBS0 ≈ 0.0002 and the standard ngsRelate estimator's
    SE at this sample size is sub-0.001, so realistic noise leaves
    PO well below the 0.005 ibs0_po_max threshold. The noise model
    is the dominant lever on per-pair classification recall — tighter
    sigma → higher recall, matching reality at higher nSites.
    """
    n_sites_total: int = 100000
    sigma_theta: float = 0.010
    sigma_ibs0: float = 0.0008
    # Per-chromosome simulation
    n_chromosomes: int = 30
    low_data_chrom_fraction: float = 0.0
    low_data_n_sites: int = 500


def _truncated_normal(rng: random.Random, mu: float, sigma: float,
                      lo: float, hi: float) -> float:
    for _ in range(20):
        x = rng.gauss(mu, sigma)
        if lo <= x <= hi:
            return x
    return max(lo, min(hi, mu))


def simulate_pair_coeffs(
    rel: str, rng: random.Random, cfg: SimConfig
) -> Dict[str, float]:
    mu_theta, mu_ibs0 = expected_coeffs(rel)
    theta = _truncated_normal(rng, mu_theta, cfg.sigma_theta, 0.0, 0.6)
    ibs0 = _truncated_normal(rng, mu_ibs0, cfg.sigma_ibs0, 0.0, 0.3)
    # KING-robust ~ kinship for our purposes; add small noise.
    king = _truncated_normal(rng, mu_theta, cfg.sigma_theta, -0.05, 0.6)
    return {
        "nSites": cfg.n_sites_total,
        "theta": theta,
        "IBS0": ibs0,
        "KING": king,
    }


def simulate_pairwise_table(
    pedigree: Pedigree,
    rng: random.Random,
    cfg: Optional[SimConfig] = None,
) -> List[Dict]:
    """Return a list of dict rows ready to be written as a Stage-1 .res-style
    TSV. Each pair gets its true relationship class and a noisy
    coefficient draw consistent with that class.
    """
    cfg = cfg or SimConfig()
    samples = [i.sample_id for i in pedigree.individuals]
    rows: List[Dict] = []
    for i, a in enumerate(samples):
        for j, b in enumerate(samples[i+1:], start=i+1):
            rel = true_relationship(pedigree, a, b)
            row = simulate_pair_coeffs(rel, rng, cfg)
            row["a"] = i
            row["b"] = j
            row["sample_a"] = a
            row["sample_b"] = b
            row["true_relationship"] = rel
            rows.append(row)
    return rows


# ----------------------------------------------------------------------
# Karyotype simulator (Mendelian band inheritance).
# ----------------------------------------------------------------------


def simulate_karyotype(
    pedigree: Pedigree,
    inversion_id: str,
    rng: random.Random,
    *,
    inv_allele_freq: float = 0.4,
) -> Dict[str, int]:
    """Mendelianly inherit a biallelic inversion (REF=0, INV=1).
    Returns {sample_id: band} where band ∈ {0, 1, 2}:
      band 0 = HOM_REF
      band 1 = HET
      band 2 = HOM_INV
    Founders draw alleles from Bernoulli(inv_allele_freq).
    Offspring inherit one allele per parent at random.
    """
    haps: Dict[str, Tuple[int, int]] = {}
    # Topological order: founders first.
    pending = list(pedigree.individuals)
    while pending:
        progress = False
        next_pending = []
        for ind in pending:
            if ind.father is None and ind.mother is None:
                a = 1 if rng.random() < inv_allele_freq else 0
                b = 1 if rng.random() < inv_allele_freq else 0
                haps[ind.sample_id] = (a, b)
                progress = True
                continue
            # Need both parents present to inherit cleanly.
            f_hap = None
            m_hap = None
            if ind.father is not None and ind.father in haps:
                f_hap = haps[ind.father]
            elif ind.father is None:
                # Unknown parent: draw from population freq.
                f_hap = (1 if rng.random() < inv_allele_freq else 0,
                         1 if rng.random() < inv_allele_freq else 0)
            if ind.mother is not None and ind.mother in haps:
                m_hap = haps[ind.mother]
            elif ind.mother is None:
                m_hap = (1 if rng.random() < inv_allele_freq else 0,
                         1 if rng.random() < inv_allele_freq else 0)
            if f_hap is None or m_hap is None:
                next_pending.append(ind)
                continue
            a = rng.choice(f_hap)
            b = rng.choice(m_hap)
            haps[ind.sample_id] = (a, b)
            progress = True
        pending = next_pending
        if not progress:
            # Should not happen on well-formed pedigrees.
            break
    return {sid: sum(hap) for sid, hap in haps.items()}   # 0/1/2 = band


# ----------------------------------------------------------------------
# mtDNA simulator (maternal inheritance + optional swaps).
# ----------------------------------------------------------------------


def simulate_mtdna(
    pedigree: Pedigree,
    rng: random.Random,
    *,
    swap_offspring_ids: Sequence[str] = (),
) -> Dict[str, str]:
    """Assign mtDNA haplotypes. Each founder female gets a unique
    haplotype; her matrilineal descendants share it. Founder males get
    irrelevant labels.

    ``swap_offspring_ids`` is a list of offspring whose mtDNA is
    deliberately swapped to an unrelated haplotype (to inject
    mtDNA-incompatible mother-offspring pairs for recovery tests).
    """
    out: Dict[str, str] = {}
    # founder-female haplotype counter
    hap_counter = 0
    founder_male_counter = 0

    def _new_hap() -> str:
        nonlocal hap_counter
        hap_counter += 1
        return f"M{hap_counter:03d}"

    def _new_male_hap() -> str:
        nonlocal founder_male_counter
        founder_male_counter += 1
        return f"MF{founder_male_counter:03d}"

    # Assign founders.
    by = pedigree.by_id()
    for ind in pedigree.individuals:
        if ind.father is None and ind.mother is None:
            if ind.sex == "female":
                out[ind.sample_id] = _new_hap()
            else:
                out[ind.sample_id] = _new_male_hap()

    # Propagate down (everyone except founders).
    pending = [i for i in pedigree.individuals if i.sample_id not in out]
    while pending:
        progress = False
        next_pending = []
        for ind in pending:
            mom_hap = None
            if ind.mother is not None and ind.mother in out:
                mom_hap = out[ind.mother]
            elif ind.mother is None and ind.father is not None and ind.father in out:
                # Unknown mother → assign a fresh haplotype
                mom_hap = _new_hap()
            if mom_hap is None:
                next_pending.append(ind)
                continue
            out[ind.sample_id] = mom_hap
            progress = True
        pending = next_pending
        if not progress:
            break

    # Inject swaps for incompatibility tests.
    for sid in swap_offspring_ids:
        if sid in out:
            out[sid] = _new_hap()
    return out


# ----------------------------------------------------------------------
# Stdlib shadow classifier (mirror of STEP_PED_01.classify_edge).
# ----------------------------------------------------------------------


DEFAULT_THRESHOLDS = {
    "theta_first": 0.177,
    "theta_second": 0.0884,
    "theta_third": 0.0442,
    "theta_dup_min": 0.45,
    "ibs0_po_max": 0.005,
}


def _isnan(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def classify_edge_stdlib(row: Dict, thresholds: Optional[Dict] = None) -> Dict[str, str]:
    """Stdlib mirror of STEP_PED_01_annotate_relationships.classify_edge
    (same decision order, same thresholds). Returns
    {edge_class, confidence, reasons}.
    """
    t = thresholds or DEFAULT_THRESHOLDS
    theta = row.get("theta")
    ibs0 = row.get("IBS0")
    king = row.get("KING")
    n_sites = row.get("nSites")

    reasons: List[str] = []
    confidence = "high"

    if _isnan(theta):
        return {"edge_class": "undetermined", "confidence": "low",
                "reasons": "theta_missing"}

    if n_sites is not None and not _isnan(n_sites) and n_sites < 1000:
        confidence = "low"
        reasons.append(f"low_n_sites_{int(n_sites)}")

    if theta >= t["theta_dup_min"]:
        if _isnan(ibs0):
            reasons.append("ibs0_missing_dup_inferred_from_theta")
            confidence = "medium" if confidence == "high" else confidence
        elif ibs0 >= 0.001:
            reasons.append("high_theta_unexpected_ibs0")
            confidence = "low"
        return {"edge_class": "duplicate_or_clone",
                "confidence": confidence,
                "reasons": ";".join(reasons)}

    if theta >= t["theta_first"]:
        if _isnan(ibs0):
            reasons.append("ibs0_missing_cannot_separate_po_fs")
            return {"edge_class": "ambiguous_first_degree",
                    "confidence": "low", "reasons": ";".join(reasons)}
        if ibs0 < t["ibs0_po_max"]:
            if king is not None and not _isnan(king) and king < 0.15:
                reasons.append(f"king_low_for_first_degree_{king:.3f}")
                confidence = "medium" if confidence == "high" else confidence
            return {"edge_class": "parent_offspring",
                    "confidence": confidence, "reasons": ";".join(reasons)}
        return {"edge_class": "full_sibling",
                "confidence": confidence, "reasons": ";".join(reasons)}

    if theta >= t["theta_second"]:
        return {"edge_class": "second_degree",
                "confidence": confidence, "reasons": ";".join(reasons)}
    if theta >= t["theta_third"]:
        return {"edge_class": "third_degree",
                "confidence": confidence, "reasons": ";".join(reasons)}
    return {"edge_class": "unrelated",
            "confidence": confidence, "reasons": ";".join(reasons)}


# ----------------------------------------------------------------------
# Map true relationship → expected classifier verdict.
# ----------------------------------------------------------------------


_EXPECTED_VERDICT = {
    REL_DUPLICATE: "duplicate_or_clone",
    REL_PO: "parent_offspring",
    REL_FS: "full_sibling",
    REL_HS: "second_degree",
    REL_AVUNCULAR: "second_degree",
    REL_FIRST_COUSIN: "third_degree",
    REL_UNRELATED: "unrelated",
}


def expected_verdict(rel: str) -> str:
    return _EXPECTED_VERDICT.get(rel, "unrelated")


# ----------------------------------------------------------------------
# Stage-1-output extractor (stdlib).
# ----------------------------------------------------------------------


def derive_hubs_and_roster(
    pedigree: Pedigree,
    edges: List[Dict],
) -> Tuple[List[Dict], List[Dict]]:
    """Given simulated edges with edge_class assigned, build
    (pairwise_rows, roster_rows) ready to be written as Stage-1 TSVs.

    Hub assignment is taken from family_id (perfect recovery — this
    function does not re-solve hub topology, it just translates the
    known pedigree into Stage 1 output schema).

    Role assignment uses the pedigree's father/mother fields:
      - both parents known → father/mother in two_parents hub
      - one parent known → likely_parent/forced_parent in parent_with_sibship
        (forced_parent if >=3 offspring, else likely_parent)
      - dyads (one PO pair only) → ambiguous_first_degree_PO
    """
    by_family: Dict[str, List[Individual]] = {}
    for ind in pedigree.individuals:
        by_family.setdefault(ind.family_id, []).append(ind)

    roster: List[Dict] = []
    hub_counter = 0
    for fid, members in by_family.items():
        hub_counter += 1
        hub_id = f"H{hub_counter:03d}"
        parents = [m for m in members
                   if any(o.father == m.sample_id or o.mother == m.sample_id
                          for o in members)]
        offspring = [m for m in members if m.father or m.mother]
        n_off = len(offspring)
        if len(parents) == 2 and len(offspring) >= 1:
            hub_type = "two_parents_with_sibship"
            for m in members:
                if m.sex == "male" and m in parents:
                    role = "father"
                elif m.sex == "female" and m in parents:
                    role = "mother"
                elif m in offspring:
                    role = "possible_offspring"
                else:
                    role = "ambiguous"
                roster.append({
                    "sample_id": m.sample_id, "hub_id": hub_id,
                    "hub_type": hub_type, "hub_size": len(members),
                    "possible_role": role, "role_confidence": "high",
                    "reason": "simulated_truth",
                })
        elif len(parents) == 1 and n_off >= 2:
            hub_type = "parent_with_sibship"
            parent = parents[0]
            for m in members:
                if m is parent:
                    role = "forced_parent" if n_off >= 3 else "likely_parent"
                elif m in offspring:
                    role = "possible_offspring"
                else:
                    role = "ambiguous"
                roster.append({
                    "sample_id": m.sample_id, "hub_id": hub_id,
                    "hub_type": hub_type, "hub_size": len(members),
                    "possible_role": role, "role_confidence": "high",
                    "reason": "simulated_truth",
                })
        elif len(parents) == 1 and n_off == 1:
            # One PO dyad — direction NOT recoverable from edge alone.
            hub_type = "po_dyad_only"
            for m in members:
                roster.append({
                    "sample_id": m.sample_id, "hub_id": hub_id,
                    "hub_type": hub_type, "hub_size": len(members),
                    "possible_role": "ambiguous_first_degree_PO",
                    "role_confidence": "low",
                    "reason": "dyad_no_triangulation",
                })
        else:
            hub_type = "mixed_or_complex"
            for m in members:
                roster.append({
                    "sample_id": m.sample_id, "hub_id": hub_id,
                    "hub_type": hub_type, "hub_size": len(members),
                    "possible_role": "ambiguous", "role_confidence": "low",
                    "reason": "simulated",
                })

    return edges, roster


# ----------------------------------------------------------------------
# TSV writers.
# ----------------------------------------------------------------------


def write_pairwise_tsv(
    path, edges: List[Dict], include_truth: bool = False
) -> None:
    cols = ["sample_a", "sample_b", "a", "b", "nSites", "theta", "IBS0",
            "KING", "edge_class", "confidence", "reasons"]
    if include_truth:
        cols.append("true_relationship")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in edges:
            row = []
            for c in cols:
                v = e.get(c, "")
                if isinstance(v, float):
                    row.append(f"{v:.6f}")
                else:
                    row.append(str(v))
            fh.write("\t".join(row) + "\n")


def write_roster_tsv(path, roster: List[Dict]) -> None:
    cols = ["sample_id", "hub_id", "hub_type", "hub_size",
            "possible_role", "role_confidence", "reason"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in roster:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def write_karyotype_tsv(
    path, bands: Dict[str, int], inversion_id: str
) -> None:
    with open(path, "w") as fh:
        fh.write("sample_id\tinversion_id\tband\tconfidence\n")
        for sid, b in bands.items():
            fh.write(f"{sid}\t{inversion_id}\t{b}\thigh\n")


def write_mtdna_tsv(path, haps: Dict[str, str]) -> None:
    with open(path, "w") as fh:
        fh.write("sample_id\tmtdna_haplotype\n")
        for sid, h in haps.items():
            fh.write(f"{sid}\t{h}\n")
