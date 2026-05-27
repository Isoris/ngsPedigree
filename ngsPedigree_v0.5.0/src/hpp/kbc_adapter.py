"""
KBC table B adapter — MVP 5 stub.

Consumes ``kbc_variant_arrangement_assignments.tsv`` (KBC SPEC §3.B) and
exposes per-variant arrangement-background lookups for the optional
``hpp_kbc_arrangement_crosscheck`` bloc.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class KbcAssignment:
    variant_id: str
    inversion_id: str
    pod_segment: str  # L | M | R | whole
    arrangement_background: str  # A_private | B_private | shared | unassigned
    assignment_confidence: str  # high | medium | low | unassigned


class KbcNotReadyError(NotImplementedError):
    pass


class PlaceholderKbcAdapter:
    """In-memory placeholder for HPP cross-check unit tests."""

    def __init__(self, records: Optional[dict[str, KbcAssignment]] = None):
        self._records: dict[str, KbcAssignment] = records or {}

    def add(self, a: KbcAssignment) -> None:
        self._records[a.variant_id] = a

    def lookup(self, variant_id: str) -> Optional[KbcAssignment]:
        return self._records.get(variant_id)


def load_kbc_table_b(path) -> PlaceholderKbcAdapter:
    raise KbcNotReadyError(
        "kbc_adapter.load_kbc_table_b: real-table loader is MVP 5. "
        "Use PlaceholderKbcAdapter for unit tests until then."
    )
