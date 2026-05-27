"""
I/O adapter surface for HPP.

This module defines the *interface* every adapter must satisfy so that
upstream / downstream modules can be swapped without touching the
projection / gene-status / Mendelian / cross-check logic. Two concrete
adapters ship under different conditions:

  Stage 3 inheritance maps:
    placeholder loader  ──  src/hpp/stage3_placeholder.py     (MVP 1, present)
    real loader         ──  src/hpp/stage3_real.py            (MVP 6, stub)

  Variant annotation table:
    placeholder loader  ──  src/hpp/variant_master.py         (MVP 2, stub)

  KBC arrangement-assignment table:
    placeholder loader  ──  src/hpp/kbc_adapter.py            (MVP 5, stub)

  Joint VCF genotypes:
    fixture loader      ──  src/hpp/vcf_lite.py               (MVP 1, present)
    real loader         ──  to be supplied via cyvcf2/pysam   (MVP 2)

All adapters expose ``load_*`` callables returning typed records. Schema
validation lives inside each adapter, raising the adapter's own
SchemaError subclass on contract violation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Protocol

from .stage3_placeholder import (
    DyadSegment,
    PhaseKey,
    Stage3Inputs,
    TriadSegment,
)
from .vcf_lite import VariantRecord

# ----------------------------------------------------------------------
# Adapter protocols — what the swap targets must implement.
# ----------------------------------------------------------------------


class Stage3Adapter(Protocol):
    """Stage 3 inheritance-map source. Placeholder and real both implement."""

    def load_dyad_map(self, path: str | Path) -> Iterable[DyadSegment]: ...
    def load_triad_map(self, path: str | Path) -> Iterable[TriadSegment]: ...
    def load_parent_phase(self, path: str | Path) -> Dict[PhaseKey, str]: ...
    def load(
        self,
        dyad_map_path: Optional[str | Path] = None,
        triad_map_path: Optional[str | Path] = None,
        parent_phase_path: Optional[str | Path] = None,
    ) -> Stage3Inputs: ...


class VariantMasterAdapter(Protocol):
    """variant_master_scored.tsv consumer.

    HPP only reads the consequence-annotation columns it needs. The real
    adapter is implemented at MVP 2 against the MODULE_CONSERVATION
    STEP 16 output.
    """

    def is_damaging(self, variant_id: str, tier: str = "T1") -> bool: ...
    def lookup(self, variant_id: str) -> Optional[dict]: ...


class JointVcfAdapter(Protocol):
    """Joint multisample VCF consumer.

    MVP 1 uses vcf_lite (stdlib) for synthetic fixtures. MVP 2+ swap in
    a cyvcf2- or pysam-backed adapter against the real joint VCF.
    """

    def iter_variants(self, region: Optional[str] = None) -> Iterator[VariantRecord]: ...
    def sample_ids(self) -> list[str]: ...


class KbcAdapter(Protocol):
    """KBC table B (kbc_variant_arrangement_assignments.tsv) consumer.

    Used by the optional hpp_kbc_arrangement_crosscheck bloc. MVP 5 stub.
    """

    def lookup(self, variant_id: str) -> Optional[dict]: ...


# ----------------------------------------------------------------------
# Default adapter wiring for MVP 1.
# ----------------------------------------------------------------------


def default_stage3_adapter() -> Stage3Adapter:
    from . import stage3_placeholder

    return stage3_placeholder  # module satisfies the Protocol


def default_variant_master_adapter():
    from . import variant_master

    return variant_master.PlaceholderVariantMaster()


def default_kbc_adapter():
    from . import kbc_adapter

    return kbc_adapter.PlaceholderKbcAdapter()


# ----------------------------------------------------------------------
# Output-table writers — placeholder for MVP 2+ table A/B/C/D/E emission.
# ----------------------------------------------------------------------


def write_tsv(path: str | Path, columns: list[str], rows: Iterable[dict]) -> None:
    """Generic TSV writer used by table_A_writer etc."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write(
                "\t".join("" if row.get(c) is None else str(row[c]) for c in columns)
                + "\n"
            )
