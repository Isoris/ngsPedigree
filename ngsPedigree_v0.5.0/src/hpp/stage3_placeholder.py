"""
Stage 3 placeholder loader for HPP MVP 1.

Reads the three inputs Stage 3 will eventually produce (or that the
HANDOFF flags as the open question):

  1. inheritance_map_dyad.tsv   — SPEC_HPP.md §3.2 dyad schema
  2. inheritance_map_triad.tsv  — SPEC_HPP.md §3.2 triad schema
  3. parent_phase.tsv           — sidecar for parent-heterozygous-site
                                   phase (HANDOFF.md open question #1)

When ngsPedigree Stage 3 lands, this file is replaced by
``stage3_real.py``. Downstream code consumes ``Stage3Inputs`` and is
agnostic to which loader produced it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

ALLOWED_PARENTAL_HAP = {"1", "2", "ambiguous"}
ALLOWED_PHASED_HAP = {"1", "2"}
ALLOWED_CONFIDENCE = {"Gold", "Silver", "Bronze"}

DYAD_COLUMNS = [
    "dyad_id",
    "parent_sample_id",
    "offspring_sample_id",
    "chrom",
    "seg_start",
    "seg_end",
    "parental_hap_inherited",
    "segment_confidence",
    "recomb_event_left",
    "recomb_event_right",
    "n_informative_markers",
]

TRIAD_COLUMNS = [
    "triad_id",
    "paternal_sample_id",
    "maternal_sample_id",
    "offspring_sample_id",
    "chrom",
    "seg_start",
    "seg_end",
    "paternal_hap_inherited",
    "maternal_hap_inherited",
    "segment_confidence",
    "recomb_event_left",
    "recomb_event_right",
    "n_informative_markers",
]

PHASE_COLUMNS = ["parent_sample_id", "chrom", "pos", "ref", "alt", "parental_hap"]


@dataclass(frozen=True)
class DyadSegment:
    dyad_id: str
    parent_sample_id: str
    offspring_sample_id: str
    chrom: str
    seg_start: int
    seg_end: int
    parental_hap_inherited: str
    segment_confidence: str
    recomb_event_left: bool
    recomb_event_right: bool
    n_informative_markers: int
    notes: str = ""


@dataclass(frozen=True)
class TriadSegment:
    triad_id: str
    paternal_sample_id: str
    maternal_sample_id: str
    offspring_sample_id: str
    chrom: str
    seg_start: int
    seg_end: int
    paternal_hap_inherited: str
    maternal_hap_inherited: str
    segment_confidence: str
    recomb_event_left: bool
    recomb_event_right: bool
    n_informative_markers: int
    notes: str = ""


PhaseKey = Tuple[str, str]  # (parent_sample_id, variant_id)


@dataclass
class Stage3Inputs:
    dyad_segments: List[DyadSegment] = field(default_factory=list)
    triad_segments: List[TriadSegment] = field(default_factory=list)
    parent_phase: Dict[PhaseKey, str] = field(default_factory=dict)

    def segments_for_dyad(self, dyad_id: str) -> List[DyadSegment]:
        return [s for s in self.dyad_segments if s.dyad_id == dyad_id]

    def segments_for_triad(self, triad_id: str) -> List[TriadSegment]:
        return [s for s in self.triad_segments if s.triad_id == triad_id]

    def dyad_ids(self) -> List[str]:
        seen: List[str] = []
        for s in self.dyad_segments:
            if s.dyad_id not in seen:
                seen.append(s.dyad_id)
        return seen

    def triad_ids(self) -> List[str]:
        seen: List[str] = []
        for s in self.triad_segments:
            if s.triad_id not in seen:
                seen.append(s.triad_id)
        return seen


class Stage3SchemaError(ValueError):
    pass


def _read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path) as fh:
        header_line = fh.readline().rstrip("\n")
        if not header_line:
            raise Stage3SchemaError(f"{path}: empty file")
        header = header_line.split("\t")
        rows: List[Dict[str, str]] = []
        for n, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                raise Stage3SchemaError(
                    f"{path}:{n}: expected {len(header)} fields, got {len(parts)}"
                )
            rows.append(dict(zip(header, parts)))
    return header, rows


def _require_columns(path: Path, header: Iterable[str], required: Iterable[str]) -> None:
    missing = [c for c in required if c not in header]
    if missing:
        raise Stage3SchemaError(f"{path}: missing required columns: {missing}")


def _parse_bool(token: str, path: Path, line: int, col: str) -> bool:
    t = token.strip().lower()
    if t in ("true", "1", "t", "yes"):
        return True
    if t in ("false", "0", "f", "no"):
        return False
    raise Stage3SchemaError(f"{path}:{line}: column {col!r} not boolean: {token!r}")


def _parse_int(token: str, path: Path, line: int, col: str) -> int:
    try:
        return int(token)
    except ValueError as exc:
        raise Stage3SchemaError(
            f"{path}:{line}: column {col!r} not integer: {token!r}"
        ) from exc


def _validate_enum(value: str, allowed: set, path: Path, line: int, col: str) -> str:
    if value not in allowed:
        raise Stage3SchemaError(
            f"{path}:{line}: column {col!r} value {value!r} not in {sorted(allowed)}"
        )
    return value


def load_dyad_map(path: str | Path) -> List[DyadSegment]:
    path = Path(path)
    header, rows = _read_tsv(path)
    _require_columns(path, header, DYAD_COLUMNS)
    segments: List[DyadSegment] = []
    for n, row in enumerate(rows, start=2):
        seg_start = _parse_int(row["seg_start"], path, n, "seg_start")
        seg_end = _parse_int(row["seg_end"], path, n, "seg_end")
        if seg_end <= seg_start:
            raise Stage3SchemaError(
                f"{path}:{n}: seg_end ({seg_end}) must be > seg_start ({seg_start})"
            )
        segments.append(
            DyadSegment(
                dyad_id=row["dyad_id"],
                parent_sample_id=row["parent_sample_id"],
                offspring_sample_id=row["offspring_sample_id"],
                chrom=row["chrom"],
                seg_start=seg_start,
                seg_end=seg_end,
                parental_hap_inherited=_validate_enum(
                    row["parental_hap_inherited"],
                    ALLOWED_PARENTAL_HAP,
                    path,
                    n,
                    "parental_hap_inherited",
                ),
                segment_confidence=_validate_enum(
                    row["segment_confidence"],
                    ALLOWED_CONFIDENCE,
                    path,
                    n,
                    "segment_confidence",
                ),
                recomb_event_left=_parse_bool(
                    row["recomb_event_left"], path, n, "recomb_event_left"
                ),
                recomb_event_right=_parse_bool(
                    row["recomb_event_right"], path, n, "recomb_event_right"
                ),
                n_informative_markers=_parse_int(
                    row["n_informative_markers"], path, n, "n_informative_markers"
                ),
                notes=row.get("notes", ""),
            )
        )
    return segments


def load_triad_map(path: str | Path) -> List[TriadSegment]:
    path = Path(path)
    header, rows = _read_tsv(path)
    _require_columns(path, header, TRIAD_COLUMNS)
    segments: List[TriadSegment] = []
    for n, row in enumerate(rows, start=2):
        seg_start = _parse_int(row["seg_start"], path, n, "seg_start")
        seg_end = _parse_int(row["seg_end"], path, n, "seg_end")
        if seg_end <= seg_start:
            raise Stage3SchemaError(
                f"{path}:{n}: seg_end ({seg_end}) must be > seg_start ({seg_start})"
            )
        segments.append(
            TriadSegment(
                triad_id=row["triad_id"],
                paternal_sample_id=row["paternal_sample_id"],
                maternal_sample_id=row["maternal_sample_id"],
                offspring_sample_id=row["offspring_sample_id"],
                chrom=row["chrom"],
                seg_start=seg_start,
                seg_end=seg_end,
                paternal_hap_inherited=_validate_enum(
                    row["paternal_hap_inherited"],
                    ALLOWED_PARENTAL_HAP,
                    path,
                    n,
                    "paternal_hap_inherited",
                ),
                maternal_hap_inherited=_validate_enum(
                    row["maternal_hap_inherited"],
                    ALLOWED_PARENTAL_HAP,
                    path,
                    n,
                    "maternal_hap_inherited",
                ),
                segment_confidence=_validate_enum(
                    row["segment_confidence"],
                    ALLOWED_CONFIDENCE,
                    path,
                    n,
                    "segment_confidence",
                ),
                recomb_event_left=_parse_bool(
                    row["recomb_event_left"], path, n, "recomb_event_left"
                ),
                recomb_event_right=_parse_bool(
                    row["recomb_event_right"], path, n, "recomb_event_right"
                ),
                n_informative_markers=_parse_int(
                    row["n_informative_markers"], path, n, "n_informative_markers"
                ),
                notes=row.get("notes", ""),
            )
        )
    return segments


def load_parent_phase(path: str | Path) -> Dict[PhaseKey, str]:
    path = Path(path)
    header, rows = _read_tsv(path)
    _require_columns(path, header, PHASE_COLUMNS)
    out: Dict[PhaseKey, str] = {}
    for n, row in enumerate(rows, start=2):
        pos = _parse_int(row["pos"], path, n, "pos")
        hap = _validate_enum(
            row["parental_hap"], ALLOWED_PHASED_HAP, path, n, "parental_hap"
        )
        variant_id = f"{row['chrom']}:{pos}:{row['ref']}:{row['alt']}"
        key: PhaseKey = (row["parent_sample_id"], variant_id)
        if key in out and out[key] != hap:
            raise Stage3SchemaError(
                f"{path}:{n}: contradictory parent_phase for {key}: "
                f"already {out[key]!r}, now {hap!r}"
            )
        out[key] = hap
    return out


def load_stage3(
    dyad_map_path: Optional[str | Path] = None,
    triad_map_path: Optional[str | Path] = None,
    parent_phase_path: Optional[str | Path] = None,
) -> Stage3Inputs:
    inputs = Stage3Inputs()
    if dyad_map_path is not None:
        inputs.dyad_segments = load_dyad_map(dyad_map_path)
    if triad_map_path is not None:
        inputs.triad_segments = load_triad_map(triad_map_path)
    if parent_phase_path is not None:
        inputs.parent_phase = load_parent_phase(parent_phase_path)
    return inputs
