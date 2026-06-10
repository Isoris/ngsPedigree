"""
Stage 1/2 → polarization-IN extraction.

Reads ngsPedigree Stage 1 outputs (and optionally Stage 2's review
flags), extracts directional dyads and triads, and combines them with
a per-sample karyotype call table to emit the polarization IN JSON
expected by ``scripts/05_polarize_inversion.py``.

This is the bridge between the pedigree-classification half of
ngsPedigree (which works from ngsRelate coefficients) and the
inversion-polarization half (which feeds ngsTracts).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .inversion_polarization import (
    ALLOWED_BANDS,
    ALLOWED_POLARITY_HINTS,
    DyadPair,
    KaryotypeCall,
    TriadTrio,
)


PARENT_ROLES = frozenset({
    "mother", "father",
    "parent_a", "parent_b",
    "forced_parent", "likely_parent",
})
OFFSPRING_ROLES = frozenset({"possible_offspring"})
MALE_ROLES = frozenset({"father"})
FEMALE_ROLES = frozenset({"mother"})


class PedigreeExtractError(ValueError):
    pass


@dataclass(frozen=True)
class RosterRow:
    sample_id: str
    hub_id: str
    hub_type: str
    possible_role: str
    role_confidence: str


# ----------------------------------------------------------------------
# TSV readers.
# ----------------------------------------------------------------------


def _read_tsv(path: str | Path, required: Iterable[str]) -> List[Dict[str, str]]:
    path = Path(path)
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        if not header_line:
            raise PedigreeExtractError(f"{path}: empty file")
        header = header_line.split("\t")
        for col in required:
            if col not in header:
                raise PedigreeExtractError(
                    f"{path}: missing required column {col!r}"
                )
        rows: List[Dict[str, str]] = []
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                raise PedigreeExtractError(
                    f"{path}:{n}: expected {len(header)} fields, got {len(parts)}"
                )
            rows.append(dict(zip(header, parts)))
    return rows


def load_stage1_pairwise(path: str | Path) -> List[Dict[str, str]]:
    return _read_tsv(path, required=("sample_a", "sample_b", "edge_class"))


def load_stage1_roster(path: str | Path) -> Dict[str, RosterRow]:
    rows = _read_tsv(path, required=("sample_id", "hub_id", "hub_type",
                                     "possible_role"))
    out: Dict[str, RosterRow] = {}
    for r in rows:
        out[r["sample_id"]] = RosterRow(
            sample_id=r["sample_id"],
            hub_id=r["hub_id"],
            hub_type=r["hub_type"],
            possible_role=r["possible_role"],
            role_confidence=r.get("role_confidence", ""),
        )
    return out


def load_stage2_review_set(path: str | Path) -> Set[Tuple[str, str]]:
    """Return {(sample_a, sample_b)} for pairs flagged by Stage 2 QC."""
    rows = _read_tsv(path, required=("sample_a", "sample_b"))
    out: Set[Tuple[str, str]] = set()
    for r in rows:
        flag = r.get("pair_review_flag", "OK")
        if flag and flag != "OK":
            out.add((r["sample_a"], r["sample_b"]))
    return out


def load_karyotype_tsv(
    path: str | Path,
    inversion_id: str,
) -> List[KaryotypeCall]:
    """Read a karyotype-calls TSV.

    Required columns: ``sample_id``, ``band``.
    Optional columns: ``confidence``, ``inversion_id``.

    If the file has an ``inversion_id`` column, rows are filtered to
    those matching ``inversion_id``. Otherwise every row is taken to
    belong to ``inversion_id``.
    """
    rows = _read_tsv(path, required=("sample_id", "band"))
    has_inv_col = "inversion_id" in rows[0] if rows else False
    out: List[KaryotypeCall] = []
    for r in rows:
        if has_inv_col and r.get("inversion_id") != inversion_id:
            continue
        try:
            band = int(r["band"])
        except ValueError as exc:
            raise PedigreeExtractError(
                f"{path}: band must be 0/1/2; got {r['band']!r}"
            ) from exc
        if band not in ALLOWED_BANDS:
            raise PedigreeExtractError(
                f"{path}: band must be in {ALLOWED_BANDS}; got {band!r}"
            )
        out.append(KaryotypeCall(
            sample_id=r["sample_id"],
            inversion_id=inversion_id,
            band=band,
            confidence=r.get("confidence") or "high",
        ))
    return out


# ----------------------------------------------------------------------
# Direction inference and triad assembly.
# ----------------------------------------------------------------------


def _parent_sex_from_role(role: str) -> Optional[str]:
    if role in MALE_ROLES:
        return "male"
    if role in FEMALE_ROLES:
        return "female"
    return None


def extract_dyads(
    *,
    pairwise: List[Dict[str, str]],
    roster: Dict[str, RosterRow],
    review_set: Optional[Set[Tuple[str, str]]] = None,
) -> List[DyadPair]:
    """Emit one DyadPair per PO edge whose direction is decidable from
    the roster. Skip edges with ambiguous direction or with both
    samples missing from the roster."""
    dyads: List[DyadPair] = []
    for row in pairwise:
        if row.get("edge_class") != "parent_offspring":
            continue
        a, b = row["sample_a"], row["sample_b"]
        if review_set is not None and (a, b) in review_set:
            continue
        ra = roster.get(a)
        rb = roster.get(b)
        if ra is None or rb is None:
            continue
        a_is_parent = ra.possible_role in PARENT_ROLES
        b_is_parent = rb.possible_role in PARENT_ROLES
        a_is_off = ra.possible_role in OFFSPRING_ROLES
        b_is_off = rb.possible_role in OFFSPRING_ROLES
        if a_is_parent and b_is_off:
            dyads.append(DyadPair(
                parent_sample_id=a,
                offspring_sample_id=b,
                parent_sex=_parent_sex_from_role(ra.possible_role),
            ))
        elif b_is_parent and a_is_off:
            dyads.append(DyadPair(
                parent_sample_id=b,
                offspring_sample_id=a,
                parent_sex=_parent_sex_from_role(rb.possible_role),
            ))
        # ambiguous_first_degree_PO and other unresolved pairs are skipped.
    return dyads


def extract_triads(
    *,
    pairwise: List[Dict[str, str]],
    roster: Dict[str, RosterRow],
    review_set: Optional[Set[Tuple[str, str]]] = None,
) -> Tuple[List[TriadTrio], List[str]]:
    """Emit one TriadTrio per (paternal, maternal, offspring) triple
    derivable from a `two_parents_with_sibship` hub. Returns (triads,
    warnings).

    Sex disambiguation: if the hub carries `mother` + `father` roles
    use those. If only `parent_a` + `parent_b` (blind mode) use the
    convention paternal = parent_a, maternal = parent_b, with a
    warning so the caller knows the sex assignment is arbitrary.
    """
    # Build PO neighbour sets per sample (for verifying the triad is
    # supported by edges).
    po_neighbors: Dict[str, Set[str]] = {}
    for row in pairwise:
        if row.get("edge_class") != "parent_offspring":
            continue
        a, b = row["sample_a"], row["sample_b"]
        if review_set is not None and (a, b) in review_set:
            continue
        po_neighbors.setdefault(a, set()).add(b)
        po_neighbors.setdefault(b, set()).add(a)

    # Group roster by hub.
    hubs: Dict[str, List[RosterRow]] = {}
    for r in roster.values():
        hubs.setdefault(r.hub_id, []).append(r)

    triads: List[TriadTrio] = []
    warnings: List[str] = []

    for hub_id, members in hubs.items():
        if not members:
            continue
        hub_type = members[0].hub_type
        if hub_type != "two_parents_with_sibship":
            continue
        # Identify the two parents.
        mothers = [m for m in members if m.possible_role == "mother"]
        fathers = [m for m in members if m.possible_role == "father"]
        parent_as = [m for m in members if m.possible_role == "parent_a"]
        parent_bs = [m for m in members if m.possible_role == "parent_b"]
        offspring = [m for m in members if m.possible_role == "possible_offspring"]

        if len(mothers) == 1 and len(fathers) == 1:
            pat = fathers[0].sample_id
            mat = mothers[0].sample_id
        elif len(parent_as) == 1 and len(parent_bs) == 1:
            pat = parent_as[0].sample_id
            mat = parent_bs[0].sample_id
            warnings.append(
                f"hub {hub_id}: sex unknown for the two parents; "
                f"using parent_a→paternal, parent_b→maternal as convention "
                f"(sex-stratified drive test will be meaningless)"
            )
        else:
            warnings.append(
                f"hub {hub_id}: could not identify exactly two parents "
                f"(mothers={len(mothers)}, fathers={len(fathers)}, "
                f"parent_a={len(parent_as)}, parent_b={len(parent_bs)})"
            )
            continue

        # Emit one triad per offspring PO-connected to both parents.
        for o in offspring:
            sid = o.sample_id
            if pat in po_neighbors.get(sid, set()) and \
               mat in po_neighbors.get(sid, set()):
                triads.append(TriadTrio(
                    paternal_sample_id=pat,
                    maternal_sample_id=mat,
                    offspring_sample_id=sid,
                ))
    return triads, warnings


# ----------------------------------------------------------------------
# Top-level assembler.
# ----------------------------------------------------------------------


@dataclass
class PolarizationInputBundle:
    inversion_id: str
    polarity_hint: str
    karyotype_calls: List[KaryotypeCall]
    dyads: List[DyadPair]
    triads: List[TriadTrio]
    warnings: List[str]

    def to_in_json(self) -> dict:
        return {
            "schema": "ngspedigree_karyotype_calls_in_v1",
            "inversion_id": self.inversion_id,
            "polarity_hint": self.polarity_hint,
            "karyotype_calls": [
                {"sample_id": c.sample_id, "band": c.band,
                 "confidence": c.confidence}
                for c in self.karyotype_calls
            ],
            "dyads": [
                {"parent_sample_id": d.parent_sample_id,
                 "offspring_sample_id": d.offspring_sample_id,
                 "parent_sex": d.parent_sex}
                for d in self.dyads
            ],
            "triads": [
                {"paternal_sample_id": t.paternal_sample_id,
                 "maternal_sample_id": t.maternal_sample_id,
                 "offspring_sample_id": t.offspring_sample_id}
                for t in self.triads
            ],
        }


def build_polarization_input(
    *,
    stage1_edges_path: str | Path,
    stage1_roster_path: str | Path,
    karyotype_path: str | Path,
    inversion_id: str,
    polarity_hint: str,
    stage2_edges_path: Optional[str | Path] = None,
) -> PolarizationInputBundle:
    if polarity_hint not in ALLOWED_POLARITY_HINTS:
        raise PedigreeExtractError(
            f"polarity_hint must be in {ALLOWED_POLARITY_HINTS}; "
            f"got {polarity_hint!r}"
        )
    pairwise = load_stage1_pairwise(stage1_edges_path)
    roster = load_stage1_roster(stage1_roster_path)
    review_set = (load_stage2_review_set(stage2_edges_path)
                  if stage2_edges_path else None)
    dyads = extract_dyads(
        pairwise=pairwise, roster=roster, review_set=review_set,
    )
    triads, warnings = extract_triads(
        pairwise=pairwise, roster=roster, review_set=review_set,
    )
    karyotype_calls = load_karyotype_tsv(karyotype_path, inversion_id)
    return PolarizationInputBundle(
        inversion_id=inversion_id,
        polarity_hint=polarity_hint,
        karyotype_calls=karyotype_calls,
        dyads=dyads,
        triads=triads,
        warnings=warnings,
    )
